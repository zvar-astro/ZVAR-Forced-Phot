#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "hdf5.h"

extern "C" {
#include "matchfile_pm.h"
}

/* C implementation of the analysis of variance algorithm for GPUs
 * with the data loaded directly from the ZTF matchfiles
 * 
 * Usage:
 * ./AOV_matchfile MATCHFILE PERIOD_MIN PERIOD_MAX
 * 
 * NOTES:
 * 
 * * USE cuda-memcheck --leak-check full to check for possible memory
 *   leaks (note that valgrind reports false-positives for a buch of
 *   CUDA stuff), see:
 *   https://stackoverflow.com/questions/20593450/valgrind-and-cuda-are-
 *   reported-leaks-real
 * * NOTE, as of CUDA 8, double-precision atomicAdd() is implemented
 *   in CUDA with hardware support in SM_6X (Pascal) GPUs. So make sure
 *   you compile the code for sm_60 (Tesla P100) or sm_61 (consumer
 *   Pascal) or sm_70 (Tesla V100, Titan V) target architecture:
 *   https://forums.developer.nvidia.com/t/why-does-atomicadd-not-work-
 *   with-doubles-as-input/56429
 * * USE nvidia-smi to check which GPU device is not busy
 * * USE nvprof for benchmarking:
 *   nvprof ./AOV_GPU 10853632043632_r.dat 0.02 20.0

 * 
 * COMPILATION:
 * nvcc -arch=sm_60 -o AOV AOV_matchfile.cu matchfile_pm.c 
 * -I /data/pmroz/bin/include/ -L /data/pmroz/bin/lib/ -lhdf5
 * You need to have the HDF5 C library installed.
 * You may need to export relevant LD_LIBRARY_PATH and PATH:
 * export PATH=$PATH:/usr/local/cuda-9.1/bin
 * export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-9.1/lib64
 * 
 * P. Mroz @ Caltech, 18 Jun 2020
 * 
 * 21 Dec 2020: handling new forced-phot ZTF matchfile format
 */

#define GPU_ID 1 // I am using GPU No 1 because No 0 was busy
#define BIN_NUMBER 10 // number of phase bins
#define MIN_GOOD_PTS 50 // minimal number of data points 
#define THREADS 16 // the actual number of threads is THREADS**2
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, 
    bool abort=true) {
    /* Check for GPU errors:
     * https://stackoverflow.com/questions/14038589/what-is-the-canonical-
     * way-to-check-for-errors-using-the-cuda-runtime-api/14038590#14038590
     */
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: %s %s %d\n",cudaGetErrorString(code),
            file,line);
        if (abort) exit(code);
    }
}

__global__
void calculate_phase (float *phase, double *time, double *freq, 
    int npts, int nfreq) {
    /* Kernel function to calculate orbital phases for all data points
     * and frequencies.
     * phase is a 1D array with npts*nfreq elements
     * phase[i*npts+j] points to the i-th frequency and the j-th data 
     * point (i=0..nfreq, j=0..npts)
     */
    int i,j;
    double aux;
    
    i = blockIdx.x*blockDim.x + threadIdx.x; // freqency
    j = blockIdx.y*blockDim.y + threadIdx.y; // data point
    if (i < nfreq && j < npts) {
        aux = (time[j]-time[0])*freq[i];
        phase[i*npts+j] = (float)(aux-floor(aux));
    }
}

__global__
void calculate_phase_bins (float *sum1, float *sum2, int *n, 
    float *phase, float *mag, int npts, int nfreq, int r) {
    /* Kernel function to calculate the sum of (squared) magnitudes
     * for all frequencies and phase bins.
     * sum1, sum2, n are 1D arrays with nfreq*r elements
     * sum1[i*r+j] points to the i-th frequency and the j-th phase bin
     * (i=0..nfreq, j=0..r)
     * 
     * NOTE: You need to use atomic functions to ensure that all
     * calculations are performed without interference from other
     * threads (no other thread can access a given memory address until
     * the operation is complete). For details, see:
     * https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
     */
    int i,j,idx,index;
    
    i = blockIdx.x*blockDim.x + threadIdx.x; // frequency
    j = blockIdx.y*blockDim.y + threadIdx.y; // data point
    if (i < nfreq && j < npts && (!isnan(mag[j]))) {
        idx = (int)(phase[npts*i+j]*r);
        index = r*i+idx;
        atomicAdd(&sum1[index],mag[j]);
        atomicAdd(&sum2[index],mag[j]*mag[j]);
        atomicAdd(&n[index],1);
    }
}

__global__
void calculate_phase_bins_2 (float *s1, float *s2, float *sum1, 
    float *sum2, int *n, float avg, int nfreq, int r) {
    /* Kernel function to calculate variances s1 and s2 for all
     * frequencies.
     * s1, s2 are 1D arrays with nfreq elements
     */
    int i,j,idx;
    double aux;
    
    i = blockIdx.x*blockDim.x + threadIdx.x; // frequency
    j = blockIdx.y*blockDim.y + threadIdx.y; // phase bin
    if (i < nfreq && j < r) {
        idx = r*i+j;
        if (n[idx] != 0) {
            aux = sum1[idx]/(float)n[idx];
            atomicAdd(&s1[i],n[idx]*(aux-avg)*(aux-avg));
            atomicAdd(&s2[i],sum2[idx]-n[idx]*aux*aux);
        }
    }
}

__global__
void AOV (float *power, float *s1, float *s2, int npts, int nfreq, 
    int r) {
    /* Kernel function to calculate AOV statistics for all frequencies.
     * power is a 1D array with nfreq elements
     */
     int i;
     
     i = blockDim.x * blockIdx.x + threadIdx.x;
     if (i < nfreq) {
         power[i] = s1[i]/s2[i];
         power[i] *= (float)(npts-r);
         power[i] /= (float)(r-1);
     }
}

float calculate_mean_magnitude (float *mag, int npts) {
    int i,ile;
    float mean_mag = 0.0;
    
    ile = 0;
    for (i=0; i<npts; i++) {
        if (!isnan(mag[i])) {mean_mag += mag[i]; ile +=1;}
    }
    mean_mag /= (float)ile;
    return mean_mag;
}

int main(int argc, char *argv[]) {
    
    int n,n_stars,n_obs,i,j,r,n_pts,*num;
    int nfreq;
    double *hjd_host,time_spend,*hjd,*freq,*freq_host;
    float *mag_host,*magerr_host;
    float *mag,*power,*power_host,*phase,*sum1,*sum2;
    float freq_min,freq_max,delta_freq,mean_mag,power_max,freq_best,*s1,*s2;
    float *mag_ref,*mag_ref_err,*ra,*dec,*flux,*ferr;
    unsigned long *gaia_id;
    unsigned short *objtype,*flag,*pid;
    cudaError_t err;
    dim3 threadsPerBlock,numBlocks;
    
    if (argc != 4) {
        fprintf(stderr,"Usage: ./read_matchfile MATCHFILE Period_min Period_max\n");
        return(0);
    }
    
    r = BIN_NUMBER; // number of phase bins
    
    err = cudaSetDevice(GPU_ID); // execute the device code on GPU_ID
    if (err != cudaSuccess) {
        fprintf(stderr,"Cannot set device ID=%d. Exiting\n",GPU_ID);
        return 1;
    }
    
    clock_t begin = clock();
    
    /* Reading data from the matchfile */
        
    read_data_from_matchfile(argv[1],&hjd_host,&pid,&gaia_id,&ra,&dec,&mag_ref,
        &mag_ref_err,&objtype,&flux,&ferr,&flag,&n_obs,&n_stars,&n);


    gpuErrchk(cudaMalloc(&hjd,n_obs*sizeof(double)));
    gpuErrchk(cudaMalloc(&mag,n_obs*sizeof(float)));
    
    cudaMemcpy(hjd,hjd_host,n_obs*sizeof(double),cudaMemcpyHostToDevice);
    
    /* Setting-up a frequency grid */

    freq_min = 1.0/(atof(argv[3]));
    freq_max = 1.0/(atof(argv[2]));
    delta_freq = 0.1/(hjd_host[n_obs-1]-hjd_host[0]);
    nfreq = (int)ceil((freq_max-freq_min)/delta_freq);

    freq_host = (double *)malloc(nfreq*sizeof(double));
    power_host = (float *)malloc(nfreq*sizeof(float));
    freq_host[0] = freq_min;
    for (i=1; i<nfreq; i++) freq_host[i] = freq_host[i-1] + delta_freq;

    gpuErrchk(cudaMalloc(&freq,nfreq*sizeof(double)));
    gpuErrchk(cudaMalloc(&power,nfreq*sizeof(float)));
    gpuErrchk(cudaMemcpy(freq,freq_host,nfreq*sizeof(double),
                         cudaMemcpyHostToDevice));
    
    /* Allocating device memory for the calculations */

    gpuErrchk(cudaMalloc(&phase,nfreq*n_obs*sizeof(float)));
    gpuErrchk(cudaMalloc(&sum1,nfreq*r*sizeof(float)));
    gpuErrchk(cudaMalloc(&sum2,nfreq*r*sizeof(float)));
    gpuErrchk(cudaMalloc(&num,nfreq*r*sizeof(int)));
    gpuErrchk(cudaMalloc(&s1,nfreq*sizeof(float)));
    gpuErrchk(cudaMalloc(&s2,nfreq*sizeof(float)));

    /* Calculating orbital phases: you need to calculate them only once */

    threadsPerBlock.x = THREADS;
    threadsPerBlock.y = THREADS;
    numBlocks.x = (nfreq+THREADS-1)/threadsPerBlock.x;
    numBlocks.y = (n_obs+THREADS-1)/threadsPerBlock.y;

    calculate_phase<<<numBlocks, threadsPerBlock>>>(phase,hjd,freq,n_obs,
        nfreq);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk(cudaDeviceSynchronize());
    
    /* Loop over all stars in the matchfile */

    mag_host = (float *) malloc(n_obs*sizeof(float));
    magerr_host = (float *) malloc(n_obs*sizeof(float));
    
    for (i=0; i<n_stars; i++) {
        n_pts = getobj_single(i,mag_host,magerr_host,flux,ferr,flag,n_obs);
        n_pts = 0;
        for (j=0; j<n_obs; j++) if (!isnan(mag_host[j])) n_pts += 1;
        if (n_pts < MIN_GOOD_PTS) continue;
        
        // copy data from CPU to GPU
        gpuErrchk( cudaMemcpy(mag,mag_host,n_obs*sizeof(float),
                              cudaMemcpyHostToDevice));
        
        // calculate mean magnitude
        mean_mag = calculate_mean_magnitude(mag_host,n_obs);
        
        // fill these arrays with zeros
        gpuErrchk(cudaMemset(sum1,0,nfreq*r*sizeof(float)));
        gpuErrchk(cudaMemset(sum2,0,nfreq*r*sizeof(float)));
        gpuErrchk(cudaMemset(num,0,nfreq*r*sizeof(int)));
        gpuErrchk(cudaMemset(s1,0,nfreq*sizeof(float)));
        gpuErrchk(cudaMemset(s2,0,nfreq*sizeof(float)));
        
        /* Calculating AOV statistics */

        numBlocks.y = (n_obs+THREADS-1)/threadsPerBlock.y;

        calculate_phase_bins<<<numBlocks, threadsPerBlock>>>(sum1,sum2,
            num,phase,mag,n_obs,nfreq,r);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        numBlocks.y = (r+THREADS-1)/threadsPerBlock.y;
        
        calculate_phase_bins_2<<<numBlocks, threadsPerBlock>>>(s1,s2,
            sum1,sum2,num,mean_mag,nfreq,r);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        threadsPerBlock.x = THREADS*THREADS;
        threadsPerBlock.y = 1;
        numBlocks.x = (nfreq+threadsPerBlock.x-1)/threadsPerBlock.x;
        numBlocks.y = 1;
        
        AOV<<<numBlocks, threadsPerBlock>>>(power,s1,s2,n_obs,nfreq,r);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        // copy AOV stat from GPU to CPU
        gpuErrchk( cudaMemcpy(power_host,power,nfreq*sizeof(float),
                              cudaMemcpyDeviceToHost) );

        /* Finding the strongest frequency */

        power_max = 0.0;
        for (j=0; j<nfreq; j++) {
        if (power_host[j] > power_max) {
                power_max = power_host[j];
                freq_best = freq_host[j];
            }
        }
        printf("%-20ld %10.6f %10.6f %7.3f %5.3f %1d %3d %10.6f %8.3f\n",
            gaia_id[i],ra[i],dec[i],mag_ref[i],mag_ref_err[i],objtype[i],
            n_pts,1.0/freq_best,power_max);
        
    }

    clock_t end = clock();
    time_spend = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("# Execution time: %.3f\n",time_spend);
    printf("# Number of stars: %d\n",n_stars);
    printf("# Number of epochs: %d\n",n_obs);

    /* Deallocate memory */

    free(hjd_host);
    free(gaia_id);
    free(mag_ref);
    free(mag_ref_err);
    free(objtype);
    free(flux);
    free(ferr);
    free(flag);
    free(ra);
    free(dec);
    free(magerr_host);
    cudaFree(hjd);
    cudaFree(mag);
    cudaFree(freq);
    cudaFree(phase);
    cudaFree(sum1);
    cudaFree(sum2);
    cudaFree(num);
    cudaFree(power);
    cudaFree(s1);
    cudaFree(s2);
    free(power_host);
    free(freq_host);
    free(mag_host);
    
    return(0);
    
}

