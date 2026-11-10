#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf5.h"
#include "matchfile_pm.h"

/* Example code for reading data from forced-phot ZTF matchfiles
 * 
 * Compilation:
 * gcc -o read_data read_data.c matchfile_pm.c -I /usr/local/hdf5/include/ -L /usr/local/hdf5/lib/ -lhdf5 -lm
 * 
 * P. Mroz @ Caltech, 21 Dec 2020
 */

int main(int argc, char *argv[]) {
    
    int i,n_pts,n,n_stars,n_obs;
    double *bjd,*mag,*magerr;
    float *ra,*dec,*mag_ref,*mag_ref_err,*flux,*ferr;
    unsigned short *objtype,*flag,*pid;
    unsigned long *gaia_id;
    
    if (argc != 2) {
        fprintf(stderr,"Usage: ./read_matchfile filename\n");
        return(0);
    }

    /* Reading the matchfile content into memory */
    read_data_from_matchfile(argv[1],&bjd,&pid,&gaia_id,&ra,&dec,&mag_ref,
        &mag_ref_err,&objtype,&flux,&ferr,&flag,&n_obs,&n_stars,&n);
    
    /* Allocating memory to store the light curve data */
    mag = (double *) malloc(sizeof(double)*n_obs);
    magerr = (double *) malloc(sizeof(double)*n_obs);
    
    /* Loop over all stars */
    for (i=0; i<10; i++) {
        // grab data for i-th star
        n_pts = getobj_double(i,mag,magerr,flux,ferr,flag,n_obs);
        //if (n_pts > 50) for (j=0; j<n_obs; j++) printf("%f %f\n",mag[j],magerr[j]);
        printf("%d %d %f\n",i,n_pts,mag_ref[i]);
    }
    
    /* Remember to free the memory */
    free(mag);
    free(magerr);
    free(bjd);
    free(gaia_id);
    free(ra);
    free(dec);
    free(mag_ref);
    free(mag_ref_err);
    free(objtype);
    free(flux);
    free(ferr);
    free(flag);
    free(pid);
    
    return(0);
}


