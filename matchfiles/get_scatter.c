#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf5.h"
#include "matchfile_pm.h"

/* Example code for reading data from forced-phot ZTF matchfiles
 * 
 * Compilation:
 * gcc -o get_scatter get_scatter.c matchfile_pm.c -I /usr/local/hdf5/include/ -L /usr/local/hdf5/lib/ -lhdf5 -lm
 * 
 * P. Mroz @ Caltech, 21 Dec 2020
 */

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))
#define MIN_NUM_PTS 10
#define MAG_ZP 26.275

double kth_smallest(double *a, int n, int k) {
    
    int i,j,l,m ;
    double x, tmp ;
    double *copy;
    
    copy = (double *) malloc(sizeof(double)*n);
    for (i=0; i<n; i++) copy[i] = a[i];

    l=0 ; m=n-1 ;
    while (l<m) {
        x=copy[k] ;
        i=l ;
        j=m ;
        do {
            while (copy[i]<x) i++ ;
            while (x<copy[j]) j-- ;
            if (i<=j) {
                tmp = copy[i];
                copy[i] = copy[j];
                copy[j] = tmp;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    tmp = copy[k];
    
    free(copy);
    return tmp;
}

int calc_stats (double *flux_mean, double *flux_rms, int *n, 
    double *flux, double *ferr, int n_pts) {
    
    int i,ng;
    double med_ferr,med_flux,mad_flux,sum1,sum2,*diff;
    
    *n = 0;
    *flux_mean = 0.0;
    *flux_rms = 0.0;
    
    // Rejecting NAN values
    ng = 0;
    for (i=0; i<n_pts; i++) {
        if (isnan(flux[i]) || isnan(ferr[i])) continue;
        flux[ng] = flux[i];
        ferr[ng] = ferr[i];
        ng += 1;
    }
    if (ng < MIN_NUM_PTS) return -1;
    
    // Calculate median brightness, median error bars and MAD
    
    med_ferr = median(ferr,ng);
    med_flux = median(flux,ng);
    
    diff = (double *) malloc(sizeof(double)*ng);
    for (i=0; i<ng; i++) diff[i] = fabs(flux[i]-med_flux);
    mad_flux = median(diff,ng);
      
    // Select good measurements
    
    mad_flux *= 7.5;
    med_ferr *= 7.5;
    sum1 = 0.0;
    sum2 = 0.0;
    *n = 0;
    
    for (i=0; i<ng; i++) {
        if (ferr[i] >= med_ferr || diff[i] >= mad_flux) continue;
        sum1 += flux[i];
        sum2 += flux[i]*flux[i];
        *n += 1;
    }
    free(diff);
    
    if (*n < MIN_NUM_PTS) return -1;
    
    *flux_mean = sum1/(double)(*n);
    *flux_rms = sqrt(sum2/(double)(*n)-(*flux_mean)*(*flux_mean));
    
    return 1;
}

int main(int argc, char *argv[]) {
    
    int i,n_pts,n,n_stars,n_obs,status;
    double *bjd,*flux,*ferr,flux_mean,flux_rms,flux_ref,mag,err;
    float *ra,*dec,*mag_ref,*mag_ref_err,*flux_all,*ferr_all;
    unsigned short *objtype,*flag,*pid;
    unsigned long *gaia_id;
    
    if (argc != 2) {
        fprintf(stderr,"Usage: ./read_matchfile filename\n");
        return(0);
    }

    /* Reading the matchfile content into memory */
    read_data_from_matchfile(argv[1],&bjd,&pid,&gaia_id,&ra,&dec,&mag_ref,
        &mag_ref_err,&objtype,&flux_all,&ferr_all,&flag,&n_obs,&n_stars,&n);
    
    /* Allocating memory to store the light curve data */
    flux = (double *) malloc(sizeof(double)*n_obs);
    ferr = (double *) malloc(sizeof(double)*n_obs);
    
    /* Loop over all stars */
    for (i=0; i<n_stars; i++) {
        if (objtype[i] != 1) continue;
        n_pts = getobj_double(i,flux,ferr,flux_all,ferr_all,flag,n_obs);
        status = calc_stats(&flux_mean,&flux_rms,&n_pts,flux,ferr,n_obs);
        if (n_pts < MIN_NUM_PTS) continue;
        flux_ref = pow(10.0,-0.4*mag_ref[i]);
        flux_ref += flux_mean;
        if (flux_ref < 0.0) continue;
        mag = MAG_ZP - 2.5*log10(flux_ref);
        err = (2.5/M_LN10)*flux_rms/flux_ref;
        printf("%10.6f %10.6f %6.3f %8.6f %4d\n",ra[i],dec[i],mag,err,n_pts);
    }
    
    /* Remember to free the memory */
    free(bjd);
    free(gaia_id);
    free(ra);
    free(dec);
    free(mag_ref);
    free(mag_ref_err);
    free(objtype);
    free(flux);
    free(ferr);
    free(flux_all);
    free(ferr_all);
    free(flag);
    free(pid);
    
    return(0);
}


