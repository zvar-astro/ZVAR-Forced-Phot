#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "hdf5.h"
#include "matchfile_pm.h"

/* Example code for reading data from forced-phot ZTF matchfiles
 * 
 * Compilation:
 * gcc -Wall -o getobj getobj.c matchfile_pm.c -I /usr/local/hdf5/include/ -L /usr/local/hdf5/lib/ -lhdf5 -lm
 * P. Mroz @ Caltech, 21 Dec 2020
 * Updated: 8 Feb 2021 (transforming flux to magnitudes)
 *
 * Modified by Zach Vanderbosch @ Caltech, 2022-2024
 */

int main(int argc, char *argv[]) {
    
    int opt,print_mags,print_metadata,print_cmds;
    int i,n_obs,success;
    double *bjd,*mzpsci,*mzpsciunc,*mzpscirms;
    double *clrco,*clrcounc;
    float *flux,*ferr,*exptime,*seeing,*maglim,*airmass;
    float mag_ref,mag_zp,flux_ref,fi,fei;
    unsigned short *flag,*pid,*filterid,*imstat,typ_ref;
    unsigned short *field,*ccd,*quad;
    unsigned long *infobits,*nps1matches,target_id;
    char *ptr;
    
    print_mags = 0;
    print_metadata = 0;
    while ((opt = getopt(argc,argv,"md")) != -1) {
        switch (opt) {
            case 'm': print_mags = 1; break;
            case 'd': print_metadata = 1; break;
            default: 
                print_mags = 0; 
                print_metadata = 0;
                break;
        }
    }
    print_cmds = print_mags + print_metadata;

    if (argc != 3 && argc != 4 && argc != 5) {
        fprintf(stderr,"Usage: ./getobj [-m] [-d] filename ident\n");
        return(0);
    }
    if (argc == 4 && print_cmds != 1) {
        fprintf(stderr,"Usage: ./getobj [-m] [-d] filename ident\n");
        return(0);
    }
    if (argc == 5 && print_cmds != 2) {
        fprintf(stderr,"Usage: ./getobj [-m] [-d] filename ident\n");
        return(0);
    }

    /* Check if the matchfile exists */
    if (access(argv[optind],F_OK) != 0) {
        fprintf(stderr,"File %s does not exist\n",argv[optind]);
        return(0);
    }

    target_id = strtoul(argv[optind+1],&ptr,10);

    /*success = getobj(argv[optind],target_id,&bjd,&pid,&flux,&ferr,&flag,
        &n_obs,&mag_ref,&typ_ref,&mag_zp);*/
    success = getobj(argv[optind],target_id,&bjd,&filterid,&exptime,&pid,&field,&ccd,&quad,&imstat,
        &infobits,&seeing,&mzpsci,&mzpsciunc,&mzpscirms,&clrco,&clrcounc,&maglim,
        &airmass,&nps1matches,&flux,&ferr,&flag,&n_obs,&mag_ref,
        &typ_ref,&mag_zp);
    
    if (success > 0) {

        printf("#\n");
        printf("# Instrumental magnitude: %.3f mag\n",mag_ref);
        printf("# DoPHOT object type: %d\n",typ_ref);
        printf("#\n");
        printf("# Column Definitions:\n");
        printf("# -------------------\n");
        printf("# bjd = Mid-exposure barycentric julian date, corrected using reference image central RA and Dec [days]\n");

        if (print_mags == 0) {

            printf("# diffimflux = Forced difference image PSF flux [DN]\n");
            printf("# diffimfluxunc = 1-sigma uncertainty on diffimflux [DN]\n");
            printf("# flag = 0 (if diffimflux is not NaN), or 1 (if diffimflux is NaN)\n");
            printf("# pid = ZTF Program ID (1, 2, or 3)\n");

            if (print_metadata == 1) {

                printf("# filterid = 1 for g-band, 2 for r-band, or 3 for i-band\n");
                printf("# exptime = Exposure time of the image [s]\n");
                printf("# field = ZTF Field ID\n");
                printf("# ccd = ZTF CCD ID (1-16)\n");
                printf("# quad = ZTF Quadrant ID (1-4)\n");
                printf("# status = Image status (0 = bad image, 1 = good image)\n");
                printf("# infobits = Processing summary/QA bits for science image\n");
                printf("# seeing = Effective FWHM of sci image [pixels]\n");
                printf("# zpmagsci = Photometric zeropoint for science image [mag]\n");
                printf("# zpmagsciunc = 1-sigma uncertainty on zpmagsci [mag]\n");
                printf("# zpmagscirms = RMS deviation in difference between instrumental mags and PS1 calibrators [mag]\n");
                printf("# clrcoeff = Linear color coefficient from calibration\n");
                printf("# clrcoeffunc = 1-sigma uncertainty on clrcoeff\n");
                printf("# maglim = Magnitude limit of the science image [mag]\n");
                printf("# airmass = Airmass of science image\n");
                printf("# nps1matches = Number of PS1 calibrators used in initial calibration of sci image\n");
                printf("# -------------------\n");
                printf("#\n");
                printf("# Column names in order:\n");
                printf("bjd diffimflux diffimfluxunc flag pid filterid exptime field ccd quad status infobits seeing zpmagsci zpmagsciunc zpmagscirms clrcoeff clrcoeffunc maglim airmass nps1matches\n");
                for (i=0; i<n_obs; i++) {
                    printf("%14.6f %10.3f %10.3f %1d %1d %1d %5.2f %4d %2d %1d %1d %9lu %6.3f %9.6f %12.10f %8.6f %9.6f %12.10f %5.2f %5.3f %6lu\n",
                        bjd[i],flux[i],ferr[i],flag[i],pid[i],filterid[i],exptime[i],field[i],ccd[i],quad[i],imstat[i],
                        infobits[i],seeing[i],mzpsci[i],mzpsciunc[i],mzpscirms[i],clrco[i],clrcounc[i],maglim[i],
                        airmass[i],nps1matches[i]);
                }
            }
            else {
                printf("# -------------------\n");
                printf("#\n");
                printf("# Column names in order:\n");
                printf("bjd diffimflux diffimfluxunc flag pid\n");

                for (i=0; i<n_obs; i++) {
                    printf("%14.6f %10.3f %10.3f %1d %1d\n",
                        bjd[i],flux[i],ferr[i],flag[i],pid[i]);
                }
            }
        }
        if (print_mags == 1) {
            flux_ref = pow(10.0,-0.4*mag_ref);
            
            printf("# mag = Source magnitude using above intrumental magnitude to scale the difference image flux [mag]\n");
            printf("# magerr = Uncertainty on source magnitude [mag]\n");
            printf("# diffimflux = Forced difference image PSF flux [DN]\n");
            printf("# diffimfluxunc = 1-sigma uncertainty on diffimflux [DN]\n");
            printf("# flag = 0 (if diffimflux is not NaN), or 1 (if diffimflux is NaN)\n");
            printf("# pid = ZTF Program ID (1, 2, or 3)\n");

            if (print_metadata == 1) {

                printf("# filterid = 1 for g-band, 2 for r-band, or 3 for i-band\n");
                printf("# exptime = Exposure time of the image [s]\n");
                printf("# field = ZTF Field ID\n");
                printf("# ccd = ZTF CCD ID (1-16)\n");
                printf("# quad = ZTF Quadrant ID (1-4)\n");
                printf("# status = Image status (0 = bad image, 1 = good image)\n");
                printf("# infobits = Processing summary/QA bits for science image\n");
                printf("# seeing = Effective FWHM of sci image [pixels]\n");
                printf("# zpmagsci = Photometric zeropoint for science image [mag]\n");
                printf("# zpmagsciunc = 1-sigma uncertainty on zpmagsci [mag]\n");
                printf("# zpmagscirms = RMS deviation in difference between instrumental mags and PS1 calibrators [mag]\n");
                printf("# clrcoeff = Linear color coefficient from calibration\n");
                printf("# clrcoeffunc = 1-sigma uncertainty on clrcoeff\n");
                printf("# maglim = Magnitude limit of the science image [mag]\n");
                printf("# airmass = Airmass of science image\n");
                printf("# nps1matches = Number of PS1 calibrators used in initial calibration of sci image\n");
                printf("# -------------------\n");
                printf("#\n");
                printf("# Column names in order:\n");
                printf("bjd mag magerr diffimflux diffimfluxunc flag pid filterid exptime field ccd quad status infobits seeing zpmagsci zpmagsciunc zpmagscirms clrcoeff clrcoeffunc maglim airmass nps1matches\n");
                for (i=0; i<n_obs; i++) {
                    fi = flux_ref + flux[i];
                    if (fi > 0.0) {
                        fei = (2.5/M_LN10) * ferr[i]/fi;
                        fi = mag_zp - 2.5*log10(fi);
                    } else {
                        fei = NAN;
                        fi = NAN;
                    }
                    printf("%14.6f %6.3f %5.3f %10.3f %10.3f %1d %1d %1d %5.2f %4d %2d %1d %1d %9lu %6.3f %9.6f %12.10f %8.6f %9.6f %12.10f %5.2f %5.3f %6lu\n",
                        bjd[i],fi,fei,flux[i],ferr[i],flag[i],pid[i],filterid[i],exptime[i],field[i],ccd[i],quad[i],imstat[i],
                        infobits[i],seeing[i],mzpsci[i],mzpsciunc[i],mzpscirms[i],clrco[i],clrcounc[i],maglim[i],
                        airmass[i],nps1matches[i]);
                }
            }
            else {
                printf("# -------------------\n");
                printf("#\n");
                printf("# Column names in order:\n");
                printf("bjd mag magerr diffimflux diffimfluxunc flag pid\n");
                for (i=0; i<n_obs; i++) {
                    fi = flux_ref + flux[i];
                    if (fi > 0.0) {
                        fei = (2.5/M_LN10) * ferr[i]/fi;
                        fi = mag_zp - 2.5*log10(fi);
                    } else {
                        fei = NAN;
                        fi = NAN;
                    }
                    printf("%14.6f %6.3f %5.3f %10.3f %10.3f %1d %1d\n",
                        bjd[i],fi,fei,flux[i],ferr[i],flag[i],pid[i]);
                }
            }
        }
        free(flux);
        free(ferr);
        free(flag);
    }
    
    /* Remember to free the memory */
    
    free(bjd);
    free(filterid);
    free(exptime);
    free(pid);
    free(field);
    free(ccd);
    free(quad);
    free(imstat);
    free(infobits);
    free(seeing);
    free(mzpsci);
    free(mzpsciunc);
    free(mzpscirms);
    free(clrco);
    free(clrcounc);
    free(maglim);
    free(airmass);
    free(nps1matches);
    
    return(0);
}


