#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

#include "sbls.h"

/* C implementation of the sparse BLS algorithm
 * 
 * gcc -Wall -o period_find period_find.c sbls.c -lm
 * 
 * Usage: ./period_find FILENAME PERIOD_MIN PERIOD_MAX
 * 
 * P. Mroz @ Caltech, Mar 18, 2021
 */

#define NMAX 2000

int read_photometry (char *filename, double *time, double *flux, 
    double *ferr) {
    
    /* Reading data from a text file. This function returns the number
     * of data points */
    
    FILE *fp;
    int i,j;
    
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    i = 0;
    while (fscanf(fp,"%lf %lf %lf%*[^\n]\n",
        &time[i],&flux[i],&ferr[i]) != EOF) i += 1;
    fclose(fp);
    
    // Transform magnitudes to flux:
    for (j=0; j<i; j++) {
        flux[i] = pow(10.0,0.4*(18.0-flux[i]));
        ferr[i] = 0.4*M_LN10*flux[i]*ferr[i];
    }
    
    return(i);
}

double get_baseline (double *time, int n) {
    
    double min,max;
    int i;
    
    min = 1.0e10;
    max = -1.0e10;
    
    for (i=0; i<n; i++) {
        if (time[i] < min) min = time[i];
        if (time[i] > max) max = time[i];
    }
    
    return (max-min);
    
}

int main (int argc, char *argv[]) {
    
    double time[NMAX],flux[NMAX],ferr[NMAX],freq_min,freq_max,delta_freq;
    double *periods,*bls,freq,baseline,best_sr,best_period;
    int opt,verbose,n,i,n_freq;
    
    verbose = 0;
    while ((opt = getopt(argc,argv,"v")) != -1) {
        switch (opt) {
            case 'v': verbose = 1; break;
            default: verbose = 0; break;
        }
    }

    if (argc != 4 && argc != 5) {
        fprintf(stderr,"Usage: ./period_find [-v] FILENAME PERIOD_MIN PERIOD_MAX\n");
        return(0);
    }
    if (argc == 5 && verbose == 0) {
        fprintf(stderr,"Usage: ./period_find [-v] FILENAME PERIOD_MIN PERIOD_MAX\n");
        return(0);
    }
    
    n = read_photometry(argv[optind],time,flux,ferr);
    if (n < MIN_NUMBER_OF_POINTS) return(0);
    
    freq_min = 1.0/(atof(argv[optind+2]));
    freq_max = 1.0/(atof(argv[optind+1]));
    baseline = get_baseline(time,n);
    delta_freq = 0.1/baseline;
    n_freq = (int)((freq_max-freq_min)/delta_freq);

    periods = (double *)malloc(n_freq*sizeof(double));
    bls = (double *)malloc(n_freq*sizeof(double));
    freq = freq_min;
    for (i=0; i<n_freq; i++) {
        periods[i] = 1.0/freq;
        freq += delta_freq;
    }
    
    sparse_bls(periods,bls,n_freq,time,flux,ferr,n);

    best_sr = 0.0;
    best_period = 0.0;
    for (i=0; i<n_freq; i++) {
        if (bls[i] > best_sr) {best_sr = bls[i]; best_period = periods[i];}
    }
    
    if (verbose == 0) {
        printf("%-20s %12.8f %.4f\n",argv[1],best_period,best_sr);
    }
    else if (verbose == 1) {
        for (i=0; i<n_freq; i++) {
            printf("%12.8f %9.4f\n",periods[i],bls[i]);
        }
    }
    free(periods);
    free(bls);
    
    return(0);
    
}
