#include <stdio.h>
#include <math.h>
#include <time.h>
#include "fitsio.h"

int read_fits_file (char *filename, double **img, int *nx, int *ny) {
    
    fitsfile *fp;
    int status,hdunum,hdutype;
    long naxes[2],fpixel[2],nelements;
    double nan_val = nan("");
    
    status = 0; /* CFITSIO status value MUST be initialized to zero */
    fpixel[0] = fpixel[1] = 1;
    
    fits_open_file(&fp,filename,READONLY,&status);
    
    fits_get_num_hdus(fp,&hdunum,&status);
    if (hdunum == 2) fits_movabs_hdu(fp,2,&hdutype,&status);
    
    fits_get_img_size(fp,2,naxes,&status);
    *nx = (int)naxes[0];
    *ny = (int)naxes[1];
    nelements = naxes[0]*naxes[1];
    
    (*img) = (double *) malloc(nelements*sizeof(double));
    
    fits_read_pix(fp,TDOUBLE,fpixel,nelements,&nan_val,*img,0,&status);
    
    fits_close_file(fp,&status);
    
    return status;
}

double interp2d (double x, double y, double x0, double y0, double *val, 
    int nx, int ny) {
    
    int ix,iy;
    double res,dx,dy;
    
    dx = x-x0+0.5+nx/2;
    dy = y-y0+0.5+ny/2;
    ix = (int)(dx);
    iy = (int)(dy);
    
    if (ix < 0 || ix >= nx) return 0.0;
    if (iy < 0 || iy >= ny) return 0.0;
    
    res = val[iy*nx+ix]*(ix+1-dx)*(iy+1-dy);
    res += val[iy*nx+ix+1]*(dx-ix)*(iy+1-dy);
    res += val[(iy+1)*nx+ix]*(ix+1-dx)*(dy-iy);
    res += val[(iy+1)*nx+ix+1]*(dx-ix)*(dy-iy);
    
    return res;
}

void phot (double *flux, double *ferr, double *img_diff, double *img_sci,
    double *psf, double *x, double *y, double gain, int nstars, int npsf, 
    int nx, int ny, int stamp) {
    
    int i,n,ngmin,kx,ky,idx,xc,yc;
    double T_ki,sum1,sum2,sigma2;
    
    ngmin = (2*stamp+1)*(2*stamp+1)/2;
    
    for (i=0; i<nstars; i++) {
        
        xc = (int)(x[i]);
        yc = (int)(y[i]);
        
        sum1 = 0.0;
        sum2 = 0.0;
        n = 0;
        
        for (kx = xc-stamp; kx <= (xc+stamp); kx++) {
            for (ky = yc-stamp; ky <= (yc+stamp); ky++) {
                
                if (kx < 0 || kx >= nx) continue;
                if (ky < 0 || ky >= ny) continue;
                
                idx = kx + ky*nx;
                if (isnan(img_diff[idx])) continue;
                
                T_ki = interp2d(kx,ky,x[i],y[i],psf,npsf,npsf);
                
                sigma2 = img_sci[idx]/gain;
                sum1 += T_ki*T_ki/sigma2;
                sum2 += img_diff[idx]*T_ki/sigma2;
                n += 1;
                
            }
        }
        
        if (sum1 > 0.0 && n > ngmin) {
            flux[i] = sum2/sum1;
            ferr[i] = 1.0/sqrt(sum1);
        } else {        
            flux[i] = NAN;
            ferr[i] = NAN;
        }
        
    }
        
}

void phot2 (double *flux, double *ferr, char *diff_imagename, 
    char *sci_imagename, char *psf_imagename, double *x, double *y, 
    double gain, int nstars, int stamp) {
    
    double *img_diff,*img_psf,*img_sci;
    int nx_psf,ny_psf,nx_diff,ny_diff,nx_sci,ny_sci;
    
    /* Read data from FITS files */
    read_fits_file(psf_imagename,&img_psf,&nx_psf,&ny_psf);
    read_fits_file(diff_imagename,&img_diff,&nx_diff,&ny_diff);
    read_fits_file(sci_imagename,&img_sci,&nx_sci,&ny_sci);
    
    /* Calculate forced photometry */
    phot(flux,ferr,img_diff,img_sci,img_psf,x,y,gain,nstars,nx_psf,
        nx_diff,ny_diff,stamp);

    free(img_diff);
    free(img_psf);
    free(img_sci);

}

