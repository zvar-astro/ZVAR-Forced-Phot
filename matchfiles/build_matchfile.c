#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "hdf5.h"
#include "matchfile_pm.h"

/* Writing forced photometry to HDF5 files
 * 
 * Compilation:
 * gcc -o build_matchfile build_matchfile.c matchfile_pm.c -I /usr/local/hdf5/include/ -L /usr/local/hdf5/lib/ -lhdf5 -lm
 * 
 * Usage:
 * build_matchfile starlist_file dophot_output_file exposure_list 
 *  mag_zp filename
 * 
 * P. Mroz @ Caltech, 13 Jan 2021
 *
 * Modified by Zach Vanderbosch @ Caltech, 2022-2024
 */

#define NMAX 1000000
#define NMAXEXP 2000
#define STRING_LEN 128
/* 0 - no compression, 1 - best compression speed, least compression, 
 * 9 - best compression ratio, slowest speed
 */
#define COMPRESSION_LEVEL 1 
 
int read_starlist (char *filename, unsigned long *gaia_id, float *ra,
    float *dec) {
    
    FILE *fp;
    int i;
    
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    i = 0;
    while (fscanf(fp,"%lu %f %f %*[^\n]\n",&gaia_id[i],&ra[i],&dec[i]) != EOF) i += 1;
    fclose(fp);
    
    return(i);
    
}

int read_dophot_output (char *filename, float *mag_ref, float *mag_ref_err,
    unsigned short *objtype, int n_stars) {
        
    
    FILE *fp;
    int i;
    
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    for (i=0; i<n_stars; i++) {
        fscanf(fp,"%*d %hu %*f %*f %f %f %*[^\n]\n",&objtype[i],&mag_ref[i],
        &mag_ref_err[i]);
    }
    fclose(fp);
    
    return(i);
    
}

int read_exposure_table (char *filename, double *jd, double *bjd, 
    unsigned short *filterid, float *exptime, unsigned short *pid, unsigned 
    short *field, unsigned short *ccd, unsigned short *quad, unsigned
    short *imstat, unsigned long *infobits, float *seeing, 
    double *mzpsci, double *mzpsciunc, double *mzpscirms,
    double *clrco, double *clrcounc, float *maglim, float *airmass,
    unsigned long *nps1matches, char **qcomment, char **files) {

    FILE *fp;
    char tmp1[STRING_LEN];
    char tmp2[STRING_LEN];
    int i;
    
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    i = 0;
    while (fscanf(fp,"%lf %lf %hu %f %hu %hu %hu %hu %hu %lu %f %lf %lf %lf %lf %lf %f %f %lu %s %s",
        &jd[i],&bjd[i],&filterid[i],&exptime[i],&pid[i],&field[i],
        &ccd[i],&quad[i],&imstat[i],&infobits[i],&seeing[i],&mzpsci[i],
        &mzpsciunc[i],&mzpscirms[i],&clrco[i],&clrcounc[i],&maglim[i],
        &airmass[i],&nps1matches[i],tmp1,tmp2) != EOF) {
        strcpy(qcomment[i],tmp1);
        strcpy(files[i],tmp2);
        i += 1;
    }
    fclose(fp);
    
    return(i);
    
}

int read_metrics_table (char *filename, unsigned long *ps_id,
    double *medianmag, double *period, double *psig, double *stdev, 
    double *ppavg, double *mssd, double *skew, double *skewp, 
    double *kurtosis, double *normev, double *vonn,double *stetsonj, 
    double *stetsonk, int n_stars) {

    FILE *fp;
    char tmp[STRING_LEN];
    int i;
    
    fp = fopen(filename,"r");
    if (fp == NULL) {
        fprintf(stderr,"Error while opening file %s\n",filename);
        return(0);
    }
    i = 0;
    while (fscanf(fp,"%lu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &ps_id[i],&medianmag[i],&period[i],&psig[i],&stdev[i],&ppavg[i],
        &mssd[i],&skew[i],&skewp[i],&kurtosis[i],&normev[i],&vonn[i],
        &stetsonj[i],&stetsonk[i]) != EOF) {
        i += 1;
    }
    fclose(fp);
    
    return(i);
    
}

hid_t save_sources_table (hid_t group, char *file_starlist, 
    char *file_dophot_output, char *file_lc_metrics, int *n_stars) {

    unsigned long *gaia_id,*ps_id;
    unsigned short *objtype;
    float *ra,*dec,*mag_ref,*mag_ref_err;
    double *medianmag,*period,*psig,*stdev;
    double *ppavg,*mssd,*skew,*skewp,*kurtosis;
    double *normev,*vonn,*stetsonj,*stetsonk;
    int i;
    sources *data;
    hid_t space,dset,dcpl;
    herr_t status;
    
    /* Allocate memory and read the starlist file */
    
    gaia_id = (unsigned long *) malloc(sizeof(unsigned long)*NMAX);
    ra = (float *) malloc(sizeof(float)*NMAX);
    dec = (float *) malloc(sizeof(float)*NMAX);
    
    *n_stars = read_starlist(file_starlist,gaia_id,ra,dec);
    
    /* Allocate memory and read DoPhot output */
    
    objtype = (unsigned short *) malloc(sizeof(unsigned short)*(*n_stars));
    mag_ref = (float *) malloc(sizeof(float)*(*n_stars));
    mag_ref_err = (float *) malloc(sizeof(float)*(*n_stars));
    
    *n_stars = read_dophot_output(file_dophot_output,mag_ref,mag_ref_err,
        objtype,*n_stars);

    /* Allocate memory and read LC Metrics file */
    
    // ps_id = (unsigned long *) malloc(sizeof(unsigned long)*(*n_stars));
    // medianmag = (double *) malloc(sizeof(double)*(*n_stars));
    // period = (double *) malloc(sizeof(double)*(*n_stars));
    // psig = (double *) malloc(sizeof(double)*(*n_stars));
    // stdev = (double *) malloc(sizeof(double)*(*n_stars));
    // ppavg = (double *) malloc(sizeof(double)*(*n_stars));
    // mssd = (double *) malloc(sizeof(double)*(*n_stars));
    // skew = (double *) malloc(sizeof(double)*(*n_stars));
    // skewp = (double *) malloc(sizeof(double)*(*n_stars));
    // kurtosis = (double *) malloc(sizeof(double)*(*n_stars));
    // normev = (double *) malloc(sizeof(double)*(*n_stars));
    // vonn = (double *) malloc(sizeof(double)*(*n_stars));
    // stetsonj = (double *) malloc(sizeof(double)*(*n_stars));
    // stetsonk = (double *) malloc(sizeof(double)*(*n_stars));
    
    // *n_stars = read_metrics_table(file_lc_metrics,ps_id,medianmag,
    //     period,psig,stdev,ppavg,mssd,skew,skewp,kurtosis,normev,
    //     vonn,stetsonj,stetsonk,*n_stars);
    
    /* Create compound data type, allocate memory and copy data */
    
    hid_t sources_tid = H5Tcreate(H5T_COMPOUND,sizeof(sources));
    H5Tinsert(sources_tid,"gaia_id",HOFFSET(sources,gaia_id),H5T_NATIVE_ULONG);
    H5Tinsert(sources_tid,"ra",HOFFSET(sources,ra),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"decl",HOFFSET(sources,dec),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"mag_ref",HOFFSET(sources,mag_ref),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"mag_err_ref",HOFFSET(sources,mag_err_ref),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"objtype",HOFFSET(sources,objtype),H5T_NATIVE_USHORT);
    // H5Tinsert(sources_tid,"medianmag",HOFFSET(sources,medianmag),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"period",HOFFSET(sources,period),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"period_sig",HOFFSET(sources,psig),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"stdev",HOFFSET(sources,stdev),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"ppavg",HOFFSET(sources,ppavg),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"mssd",HOFFSET(sources,mssd),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"skew",HOFFSET(sources,skew),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"skewp",HOFFSET(sources,skewp),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"kurtosis",HOFFSET(sources,kurtosis),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"normev",HOFFSET(sources,normev),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"vonn",HOFFSET(sources,vonn),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"stetsonj",HOFFSET(sources,stetsonj),H5T_NATIVE_DOUBLE);
    // H5Tinsert(sources_tid,"stetsonk",HOFFSET(sources,stetsonk),H5T_NATIVE_DOUBLE);

    data = (sources *) malloc(sizeof(sources)*(*n_stars));
    for (i=0; i<(*n_stars); i++) {
        data[i].gaia_id = gaia_id[i];
        data[i].ra = ra[i];
        data[i].dec = dec[i];
        data[i].mag_ref = mag_ref[i];
        data[i].mag_err_ref = mag_ref_err[i];
        data[i].objtype = objtype[i];
    }
    // data[i].medianmag = medianmag[i];
    // data[i].period = period[i];
    // data[i].psig = psig[i];
    // data[i].stdev = stdev[i];
    // data[i].ppavg = ppavg[i];
    // data[i].mssd = mssd[i];
    // data[i].skew = skew[i];
    // data[i].skewp = skewp[i];
    // data[i].kurtosis = kurtosis[i];
    // data[i].normev = normev[i];
    // data[i].vonn = vonn[i];
    // data[i].stetsonj = stetsonj[i];
    // data[i].stetsonk = stetsonk[i];
    
    hsize_t dims[1] = {*n_stars};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1024};
    
    /* Create dataspace. Setting maximum size to NULL sets the maximum 
     * size to be the current size.*/
    space = H5Screate_simple(1,dims,maxdims);
    /* Create the dataset creation property list, add the gzip compression
     * filter and set the chunk size */
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_deflate(dcpl,COMPRESSION_LEVEL);
    status = H5Pset_chunk(dcpl,1,chunk);
    
    /* Create the dataset */
    dset = H5Dcreate(group,"/data/sources",sources_tid,space,H5P_DEFAULT,
        dcpl,H5P_DEFAULT);
    
    /* Write the data to the dataset */
    status = H5Dwrite(dset,sources_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    
    /*Close and release resources, free memory */
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    
    free(data);
    free(gaia_id);
    free(ra);
    free(dec);
    free(objtype);
    free(mag_ref);
    // free(medianmag);
    // free(period);
    // free(psig);
    // free(stdev);
    // free(ppavg);
    // free(mssd);
    // free(skew);
    // free(skewp);
    // free(kurtosis);
    // free(normev);
    // free(vonn);
    // free(stetsonj);
    // free(stetsonk);
    
    return status;
    
}
    
hid_t save_exposures_table (hid_t group, double *jd, double *bjd, 
    unsigned short *filterid, float *exptime, unsigned short *pid, unsigned 
    short *field, unsigned short *ccd, unsigned short *quad, unsigned
    short *imstat, unsigned long *infobits, float *seeing, 
    double *mzpsci, double *mzpsciunc, double *mzpscirms,
    double *clrco, double *clrcounc, float *maglim, float *airmass,
    unsigned long *nps1matches, float mag_zp, int n_exp) {

    int i;
    exposures *data;
    hid_t space,dset,dcpl,dspace,attr,filtname_type;
    herr_t status;
    
    /* Create compound data type, allocate memory and copy data */

    hid_t exposures_tid = H5Tcreate(H5T_COMPOUND,sizeof(exposures));
    H5Tinsert(exposures_tid,"jd",HOFFSET(exposures,jd),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"bjd",HOFFSET(exposures,bjd),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"filterid",HOFFSET(exposures,filterid),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"exptime",HOFFSET(exposures,exptime),H5T_NATIVE_FLOAT);
    H5Tinsert(exposures_tid,"pid",HOFFSET(exposures,pid),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"field",HOFFSET(exposures,field),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"ccd",HOFFSET(exposures,ccd),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"quad",HOFFSET(exposures,quad),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"imstat",HOFFSET(exposures,imstat),H5T_NATIVE_USHORT);
    H5Tinsert(exposures_tid,"infobits",HOFFSET(exposures,infobits),H5T_NATIVE_ULONG);
    H5Tinsert(exposures_tid,"seeing",HOFFSET(exposures,seeing),H5T_NATIVE_FLOAT);
    H5Tinsert(exposures_tid,"mzpsci",HOFFSET(exposures,mzpsci),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"mzpsciunc",HOFFSET(exposures,mzpsciunc),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"mzpscirms",HOFFSET(exposures,mzpscirms),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"clrco",HOFFSET(exposures,clrco),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"clrcounc",HOFFSET(exposures,clrcounc),H5T_NATIVE_DOUBLE);
    H5Tinsert(exposures_tid,"maglim",HOFFSET(exposures,maglim),H5T_NATIVE_FLOAT);
    H5Tinsert(exposures_tid,"airmass",HOFFSET(exposures,airmass),H5T_NATIVE_FLOAT);
    H5Tinsert(exposures_tid,"nps1matches",HOFFSET(exposures,nps1matches),H5T_NATIVE_ULONG);

    data = (exposures *) malloc(sizeof(exposures)*n_exp);
    for (i=0; i<n_exp; i++) {
        data[i].jd = jd[i];
        data[i].bjd = bjd[i];
        data[i].filterid = filterid[i];
        data[i].exptime = exptime[i];
        data[i].pid = pid[i];
        data[i].field = field[i];
        data[i].ccd = ccd[i];
        data[i].quad = quad[i];
        data[i].imstat = imstat[i];
        data[i].infobits = infobits[i];
        data[i].seeing = seeing[i];
        data[i].mzpsci = mzpsci[i];
        data[i].mzpsciunc = mzpsciunc[i];
        data[i].mzpscirms = mzpscirms[i];
        data[i].clrco = clrco[i];
        data[i].clrcounc = clrcounc[i];
        data[i].maglim = maglim[i];
        data[i].airmass = airmass[i];
        data[i].nps1matches = nps1matches[i];
    }
    
    hsize_t dims[1] = {n_exp};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1024};
    
    /* Add attribute */
    dspace = H5Screate(H5S_SCALAR);
    attr = H5Acreate(group,"magzp",H5T_NATIVE_FLOAT,dspace,H5P_DEFAULT,H5P_DEFAULT);
    status = H5Awrite(attr,H5T_NATIVE_FLOAT,&mag_zp);
    H5Aclose(attr);
    H5Sclose(dspace);
    
    /* Create dataspace. Setting maximum size to NULL sets the maximum 
     * size to be the current size.*/
    space = H5Screate_simple(1,dims,maxdims);
    /* Create the dataset creation property list, add the gzip compression
     * filter and set the chunk size */
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_deflate(dcpl,COMPRESSION_LEVEL);
    status = H5Pset_chunk(dcpl,1,chunk);

    /* Create the dataset */
    dset = H5Dcreate(group,"/data/exposures",exposures_tid,space,
        H5P_DEFAULT,dcpl,H5P_DEFAULT);
    
    /* Write the data to the dataset */
    status = H5Dwrite(dset,exposures_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

    /*Close and release resources, free memory */
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    
    free(data);
        
    return status;
    
}

hid_t save_sourcedata_table (hid_t group, char **files, int n_stars, int n_exp) {

    FILE *fp;
    int i,j;
    sourcedata *data;
    hid_t space,dset,dcpl;
    herr_t status;
    
    /* Create compound data type and allocate memory */
    
    hid_t sourcedata_tid = H5Tcreate(H5T_COMPOUND,sizeof(sourcedata));
    H5Tinsert(sourcedata_tid,"flux",HOFFSET(sourcedata,flux),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flux_err",HOFFSET(sourcedata,flux_err),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flag",HOFFSET(sourcedata,flag),H5T_NATIVE_USHORT);

    data = (sourcedata *) malloc(sizeof(sourcedata)*n_exp*n_stars);
    
    /* Reading data directly from cat files */
    
    for (i=0; i<n_exp; i++) {
        fp = fopen(files[i],"rb");
        for (j=0; j<n_stars; j++) fread(&(data[j*n_exp+i].flux),sizeof(float),1,fp);
        for (j=0; j<n_stars; j++) fread(&(data[j*n_exp+i].flux_err),sizeof(float),1,fp);
        for (j=0; j<n_stars; j++) {
            if (isnan(data[j*n_exp+i].flux) || isnan(data[j*n_exp+i].flux_err)) data[j*n_exp+i].flag = 1;
            else data[j*n_exp+i].flag = 0;
        }
        fclose(fp);
    }
    
    /* Saving file */
    
    hsize_t dims[1] = {n_exp*n_stars};
    hsize_t maxdims[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {4*n_exp};
    
    /* Create dataspace. Setting maximum size to NULL sets the maximum 
     * size to be the current size.*/
    space = H5Screate_simple(1,dims,maxdims);
    /* Create the dataset creation property list, add the gzip compression
     * filter and set the chunk size */
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_deflate(dcpl,COMPRESSION_LEVEL);
    status = H5Pset_chunk(dcpl,1,chunk);

    /* Create the dataset */
    dset = H5Dcreate(group,"/data/sourcedata",sourcedata_tid,space,
        H5P_DEFAULT,dcpl,H5P_DEFAULT);
    
    /* Write the data to the dataset */
    status = H5Dwrite(dset,sourcedata_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

    /*Close and release resources, free memory */
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    
    free(data);
    
    return status;
    
}

int main (int argc, char *argv[]) {
    
    hid_t file,group;
    herr_t status;
    
    char *catfiles[NMAXEXP],*qcomment[NMAXEXP];
    double jd[NMAXEXP],bjd[NMAXEXP];
    double mzpsci[NMAXEXP],mzpsciunc[NMAXEXP],mzpscirms[NMAXEXP];
    double clrco[NMAXEXP],clrcounc[NMAXEXP];
    float seeing[NMAXEXP],maglim[NMAXEXP],exptime[NMAXEXP],airmass[NMAXEXP];
    unsigned short pid[NMAXEXP],filterid[NMAXEXP],imstat[NMAXEXP];
    unsigned short field[NMAXEXP],ccd[NMAXEXP],quad[NMAXEXP];
    unsigned long infobits[NMAXEXP],nps1matches[NMAXEXP];
    float mag_zp;
    int i,n_exp,n_stars;
    
    if (argc != 7) {
        printf("Usage: ./build_matchfile STARLIST_FILE DOPHOT_FILE METRICS_FILE EXPOSURES MAG_ZP NAME\n");
        return(0);
    }
    
    mag_zp = atof(argv[5]);
    
    for (i=0; i<NMAXEXP; i++) catfiles[i] = malloc(STRING_LEN*sizeof(char));
    for (i=0; i<NMAXEXP; i++) qcomment[i] = malloc(STRING_LEN*sizeof(char));
    n_exp = read_exposure_table(argv[4],jd,bjd,filterid,exptime,pid,
        field,ccd,quad,imstat,infobits,seeing,mzpsci,mzpsciunc,mzpscirms,
        clrco,clrcounc,maglim,airmass,nps1matches,qcomment,catfiles);
    
    file = H5Fcreate(argv[6],H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    group = H5Gcreate(file,"/data/",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    
    status = save_exposures_table(group,jd,bjd,filterid,exptime,pid,
        field,ccd,quad,imstat,infobits,seeing,mzpsci,mzpsciunc,mzpscirms,
        clrco,clrcounc,maglim,airmass,nps1matches,mag_zp,n_exp);
    status = save_sources_table(group,argv[1],argv[2],argv[3],&n_stars);
    status = save_sourcedata_table(group,catfiles,n_stars,n_exp);
    
    status = H5Gclose(group);
    status = H5Fclose(file);
    
    for (i=0; i<NMAXEXP; i++) free(catfiles[i]);
    for (i=0; i<NMAXEXP; i++) free(qcomment[i]);
        
    return(0);
}
