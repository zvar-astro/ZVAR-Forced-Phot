#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "sbls.h"

/* C implementation of the sparse BLS algorithm
 * 
 * See:
 * https://arxiv.org/abs/2103.06193
 * https://github.com/aviadpi/SparseBLS
 * 
 * P. Mroz @ Caltech, Mar 18, 2021
 */

int prepare_data (double *x_norm, double *w_norm, double *x, 
    double *sigma, int n) {
    
    /* Prepare data. Given a time series of values x_j and errors 
     * sigma_j, we first convert it to a zero-weighted-mean dataset
     * with normalized weights (x_norm, w_norm)
     * w_j = 1.0/sigma_j^2
     * W = sum w_j
     * w_norm_j = w_j / W
     * mu = sum w_norm_j * x_j = (sum w_j*x_j) / W
     * x_norm_j = x_j - mu
     */
    
    int i;
    double W,mu;
    
    W = 0.0;
    mu = 0.0;
    for (i=0; i<n; i++) {
        w_norm[i] = 1.0/(sigma[i]*sigma[i]);
        W += w_norm[i];
        mu += w_norm[i]*x[i];
    }
    mu = mu / W;
    for (i=0; i<n; i++) {
        w_norm[i] = w_norm[i] / W;
        x_norm[i] = x[i] - mu;
    }
    
    return 0;
    
}

void heap_sort_pm (double *arr, int *arr1, int N)  {
    
    /* Heap sort */
    
    if (N==0) return;

    double t;
    int t1; 
    int n = N, parent = N/2, index, child; 
    
    while (1) { 
        if (parent > 0) { 
            /* first stage - Sorting the heap */
            parent--;
            t = arr[parent];  /* save old value to t */
            t1 = arr1[parent];
        } else {
            /* second stage - Extracting elements in-place */
            n--;                /* make the heap smaller */
            if (n == 0) {
                return; /* When the heap is empty, we are done */
            }
            t = arr[n];         /* save lost heap entry to temporary */
            t1 = arr1[n];
            arr[n] = arr[0];    /* save root entry beyond heap */
            arr1[n] = arr1[0];
        }
        /* insert operation - pushing t down the heap to replace the parent */
        index = parent; /* start at the parent index */
        child = index * 2 + 1; /* get its left child index */
        while (child < n) {
            /* choose the largest child */
            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
                child++; /* right child exists and is bigger */
            }
            /* is the largest child larger than the entry? */
            if (arr[child] > t) {
                arr[index] = arr[child]; /* overwrite entry with child */
                arr1[index] = arr1[child];
                index = child; /* move index to the child */
                child = index * 2 + 1; /* get the left child and go around again */
            } else {
                break; /* t's place is found */
            }
        }
        /* store the temporary value at its new location */
        arr[index] = t; 
        arr1[index] = t1;
    }
}

void insertion_sort_pm (double *arr1, int *arr2, int n) {
    
    /* Insertion sort */
    
    double key1;
    int i,j,key2;
    
    for (i=1; i<n; i++) {
        key1 = arr1[i];
        key2 = arr2[i];
        j = i-1;
        
        while (j >= 0 && arr1[j] > key1) {
            arr1[j+1] = arr1[j];
            arr2[j+1] = arr2[j];
            j = j-1;
        }
        arr1[j+1] = key1;
        arr2[j+1] = key2;
    }
    
}

void fold_and_sort_phases (double *phase, int *phase_index, 
    double *time, double period, int n) {
    
    int i,idx;
    
    for (i=0; i<n; i++) {
        idx = phase_index[i];
        phase[i] = fmod(time[idx],period)/period;
    }
    
    //heap_sort_pm(phase,phase_index,n);
    insertion_sort_pm(phase,phase_index,n);
       
}

void sparse_bls (double *period, double *bls, int n_periods, double *jd, 
    double *flux, double *ferr, int n) {
    
    int i,k,*phase_index,i2MinPrev,i1,i2,j2,i1ln,j1ln,i2rn;
    double *x,*w,*phase,SR,SRmax,s,r,egress,ingress,sMinPrev,rMinPrev;
    double transitDurationInPhase;
    
    x = (double *)malloc(n*sizeof(double));
    w = (double *)malloc(n*sizeof(double));
    phase = (double *)malloc(n*sizeof(double));
    phase_index = (int *)malloc(n*sizeof(int));
    
    // Prepare data and initialize phase_index array
    
    prepare_data(x,w,flux,ferr,n);
    for (i=0; i<n; i++) phase_index[i] = i;
    
    // Run the main loop over all periods
    
    for (i=0; i<n_periods; i++) {
        
        SRmax = 0.0;
            
        fold_and_sort_phases(phase,phase_index,jd,period[i],n);
        
        i2MinPrev = 0;
        s = r = sMinPrev = rMinPrev = 0.0;
        
        for (i1=0; i1<n; i1++) {
                        
            i2 = i2MinPrev;
            j2 = phase_index[i2];
            s = sMinPrev;
            r = rMinPrev;
            i1ln = (i1 > 0) ? i1 - 1 : n - 1;  
            j1ln = phase_index[i1ln];
            
            if (i1 == 0) ingress = 0.5*(phase[i1ln]-1.0 + phase[i1]);
            else {
                s -= x[j1ln] * w[j1ln];
                r -= w[j1ln];
                s -= x[j2] * w[j2];
                r -= w[j2];
                ingress = 0.5*(phase[i1ln] + phase[i1]);
            }
            
            for (k=0; k<n; k++) {
                
                i2 = (i1 + k) % n;
                j2 = phase_index[i2];
                if (k < MIN_POINTS_IN_TRANSIT-1) {
                    s += x[j2] * w[j2];
                    r += w[j2];
                    i2MinPrev = i2;
                    sMinPrev = s;
                    rMinPrev = r;
                    continue;
                }
                                
                i2rn = (i2 + 1) % n;
                egress = 0.5*(phase[i2] + phase[i2rn]);
                if (i2rn < i2) egress += 0.5;
                
                transitDurationInPhase = egress-ingress;
                if (transitDurationInPhase < 0) transitDurationInPhase += 1;
                
                if (transitDurationInPhase > MAX_TRANSIT_DURATION) break;
                
                s += x[j2] * w[j2];
                r += w[j2];

                SR = s * s / (r * (1.0 - r));
                
                if (SR > SRmax) SRmax = SR; 
            } // end of k / i2
            
        } // end of i1

        bls[i] = SRmax;
        
    } // end of i (periods)

    free(x);
    free(w);
    free(phase);
    free(phase_index);
}


