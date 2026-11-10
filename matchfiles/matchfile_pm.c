#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf5.h"
#include "matchfile_pm.h"

/* Code for reading data from forced-phot ZTF matchfiles. 
 * 
 * For reading data from HDF5 files, see:
 * http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html
 * 
 * read_data_from_matchfile - reads data from a matchfile into memory
 * getobj - grabs photometry for a given object
 * 
 * P. Mroz @ Caltech, 21 Dec 2020
 *
 * Modified by Zach Vanderbosch @ Caltech, 2022-2024
 */
 
int get_number_of_rows(hid_t dataset) {
    
    /* Calculate number of rows in the dataset */
    
    hid_t dspace;
    int ndims;
    dspace = H5Dget_space(dataset);
    ndims =  H5Sget_simple_extent_ndims(dspace);
    if (ndims != 1) return(0);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace,dims,NULL);
    
    return(dims[0]);
}

int read_exposures_table(hid_t dataset, double *jd, double *bjd, 
    unsigned short *filterid, float *exptime, unsigned short *pid, 
    unsigned short *field, unsigned short *ccd, unsigned short *quad, 
    unsigned short *imstat, unsigned long *infobits, float *seeing, 
    double *mzpsci, double *mzpsciunc, double *mzpscirms,
    double *clrco, double *clrcounc, float *maglim, float *airmass,
    unsigned long *nps1matches, int n) {
    
    /* Read the entire exposures table */
    
    int i,status;
    exposures *data;
    
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
    
    data = (exposures *) malloc(sizeof(exposures)*n);
    status = H5Dread(dataset,exposures_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    
    for (i=0; i<n; i++) {
        jd[i] = data[i].jd;
        bjd[i] = data[i].bjd;
        filterid[i] = data[i].filterid;
        exptime[i] = data[i].exptime;
        pid[i] = data[i].pid;
        field[i] = data[i].field;
        ccd[i] = data[i].ccd;
        quad[i] = data[i].quad;
        imstat[i] = data[i].imstat;
        infobits[i] = data[i].infobits;
        seeing[i] = data[i].seeing;
        mzpsci[i] = data[i].mzpsci;
        mzpsciunc[i] = data[i].mzpsciunc;
        mzpscirms[i] = data[i].mzpscirms;
        clrco[i] = data[i].clrco;
        clrcounc[i] = data[i].clrcounc;
        maglim[i] = data[i].maglim;
        airmass[i] = data[i].airmass;
        nps1matches[i] = data[i].nps1matches;
    }

    free(data);
    
    return(status);
}

int read_sources_table(hid_t dataset, unsigned long *gaia_id, float *ra, 
    float *dec, float *mag_ref, float *mag_err_ref, unsigned short 
    *objtype, int n) {
    
    /* Read the entire sources table */
    
    int i,status;
    sources *data;
    
    hid_t sources_tid = H5Tcreate(H5T_COMPOUND,sizeof(sources));
    H5Tinsert(sources_tid,"gaia_id",HOFFSET(sources,gaia_id),H5T_NATIVE_ULONG);
    H5Tinsert(sources_tid,"ra",HOFFSET(sources,ra),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"decl",HOFFSET(sources,dec),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"mag_ref",HOFFSET(sources,mag_ref),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"mag_err_ref",HOFFSET(sources,mag_err_ref),H5T_NATIVE_FLOAT);
    H5Tinsert(sources_tid,"objtype",HOFFSET(sources,objtype),H5T_NATIVE_USHORT);
    
    data = (sources *) malloc(sizeof(sources)*n);
    status = H5Dread(dataset,sources_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    
    for (i=0; i<n; i++) {
        gaia_id[i] = data[i].gaia_id;
        ra[i] = data[i].ra;
        dec[i] = data[i].dec;
        mag_ref[i] = data[i].mag_ref;
        mag_err_ref[i] = data[i].mag_err_ref;
        objtype[i] = data[i].objtype;
    }
    
    free(data);
    
    return(status);
}

int read_sourcedata_table(hid_t dataset, float *flux, float *flux_err, 
    unsigned short *flag, int n) {
    
    /* Read the entire sourcedata table */
    
    int i,status;
    sourcedata *data;
    
    hid_t sourcedata_tid = H5Tcreate(H5T_COMPOUND,sizeof(sourcedata));
    H5Tinsert(sourcedata_tid,"flux",HOFFSET(sourcedata,flux),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flux_err",HOFFSET(sourcedata,flux_err),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flag",HOFFSET(sourcedata,flag),H5T_NATIVE_USHORT);
    
    data = (sourcedata *) malloc(sizeof(sourcedata)*n);
    status = H5Dread(dataset,sourcedata_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
    
    for (i=0; i<n; i++) {
        flux[i] = data[i].flux;
        flux_err[i] = data[i].flux_err;
        flag[i] = data[i].flag;
    }
    
    free(data);
    
    return(status);
}

int read_sourcedata_table_obj(int objno, hid_t dataset, float *flux, 
    float *flux_err, unsigned short *flag, int n_obs) {
    
    /* Read only part of the sourcedata table 
     * See: http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/12_Dataspaces.html
     * (Section 7.4.2)
     */
    
    hid_t dataspace,memspace;
    int i,status;
    sourcedata *data;
    
    // Define hyperslab in the dataset:
    dataspace = H5Dget_space(dataset);
    hsize_t offset[] = {objno*n_obs};
    hsize_t count[] = {n_obs};
    status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,offset,NULL,count,NULL);
    
    // Define memory hyperslab:
    hsize_t dim[] = {n_obs};
    memspace = H5Screate_simple(1,dim,NULL);
    hsize_t offset_out[] = {0};
    hsize_t count_out[] = {n_obs};
    status = H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset_out,NULL,count_out,NULL);
    
    hid_t sourcedata_tid = H5Tcreate(H5T_COMPOUND,sizeof(sourcedata));
    H5Tinsert(sourcedata_tid,"flux",HOFFSET(sourcedata,flux),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flux_err",HOFFSET(sourcedata,flux_err),H5T_NATIVE_FLOAT);
    H5Tinsert(sourcedata_tid,"flag",HOFFSET(sourcedata,flag),H5T_NATIVE_USHORT);
    
    data = (sourcedata *) malloc(sizeof(sourcedata)*n_obs);
    status = H5Dread(dataset,sourcedata_tid,memspace,dataspace,H5P_DEFAULT,data);
    
    for (i=0; i<n_obs; i++) {
        flux[i] = data[i].flux;
        flux_err[i] = data[i].flux_err;
        flag[i] = data[i].flag;
    }
    
    free(data);
    status = H5Sclose(memspace);
    
    return(status);
}


int get_exposures (hid_t group, double **jd, double **bjd, 
    unsigned short **filterid, float **exptime, unsigned short **pid, 
    unsigned short **field, unsigned short **ccd, unsigned short **quad, 
    unsigned short **imstat, unsigned long **infobits, float **seeing, 
    double **mzpsci, double **mzpsciunc, double **mzpscirms,
    double **clrco, double **clrcounc, float **maglim, float **airmass,
    unsigned long **nps1matches, int *n_obs) {
    
    /* Allocate memory and read data from exposures table */
    
    hid_t dataset;
    int status;
    
    dataset = H5Dopen2(group,"/data/exposures",H5P_DEFAULT);
    
    // get number of observations:
    *n_obs = get_number_of_rows(dataset);
    
    *jd = (double *) malloc(sizeof(double)*(*n_obs));
    *bjd = (double *) malloc(sizeof(double)*(*n_obs));
    *filterid = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *exptime = (float *) malloc(sizeof(float)*(*n_obs));
    *pid = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *field = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *ccd = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *quad = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *imstat = (unsigned short *) malloc(sizeof(unsigned short)*(*n_obs));
    *infobits = (unsigned long *) malloc(sizeof(unsigned long)*(*n_obs));
    *seeing = (float *) malloc(sizeof(float)*(*n_obs));
    *mzpsci = (double *) malloc(sizeof(double)*(*n_obs));
    *mzpsciunc = (double *) malloc(sizeof(double)*(*n_obs));
    *mzpscirms = (double *) malloc(sizeof(double)*(*n_obs));
    *clrco = (double *) malloc(sizeof(double)*(*n_obs));
    *clrcounc = (double *) malloc(sizeof(double)*(*n_obs));
    *maglim = (float *) malloc(sizeof(float)*(*n_obs));
    *airmass = (float *) malloc(sizeof(float)*(*n_obs));
    *nps1matches = (unsigned long *) malloc(sizeof(unsigned long)*(*n_obs));
    
    status = read_exposures_table(dataset,*jd,*bjd,*filterid,*exptime,*pid,
        *field,*ccd,*quad,*imstat,*infobits,*seeing,*mzpsci,*mzpsciunc,*mzpscirms,
        *clrco,*clrcounc,*maglim,*airmass,*nps1matches,*n_obs);
    status = H5Dclose(dataset);
    return status;
}

int get_sources (hid_t group, unsigned long **gaia_id, float **ra, 
    float **dec, float **mag_ref, float **mag_ref_err, unsigned short 
    **objtype, int *n_stars) {
    
    /* Allocate memory and read data from sources table */
    
    hid_t dataset;
    int status;
    
    dataset = H5Dopen2(group,"/data/sources",H5P_DEFAULT);
    
    // get number of stars:
    *n_stars = get_number_of_rows(dataset);
    
    *gaia_id = (unsigned long *) malloc(sizeof(unsigned long)*(*n_stars));
    *ra = (float *) malloc(sizeof(float)*(*n_stars));
    *dec = (float *) malloc(sizeof(float)*(*n_stars));
    *mag_ref = (float *) malloc(sizeof(float)*(*n_stars));
    *mag_ref_err = (float *) malloc(sizeof(float)*(*n_stars));
    *objtype = (unsigned short *) malloc(sizeof(unsigned short)*(*n_stars));
    
    status = read_sources_table(dataset,*gaia_id,*ra,*dec,*mag_ref,
        *mag_ref_err,*objtype,*n_stars);
    
    status = H5Dclose(dataset);
    return status;
    
}

int get_sourcedata (hid_t group, float **flux, float **ferr, unsigned 
    short **flag, int *n) {
    
    /* Allocate memory and read data from sourcedata table */
    
    hid_t dataset;
    int status;
    
    dataset = H5Dopen2(group,"/data/sourcedata",H5P_DEFAULT);
    
    // get number of data points:
    *n = get_number_of_rows(dataset);
    
    *flux = (float *) malloc(sizeof(float)*(*n));
    *ferr = (float *) malloc(sizeof(float)*(*n));
    *flag = (unsigned short *) malloc(sizeof(unsigned short)*(*n));
    
    status = read_sourcedata_table(dataset,*flux,*ferr,*flag,*n);
    
    status = H5Dclose(dataset);
    return status;
}

int get_sourcedata_obj (int objno, hid_t group, float *flux, float *ferr, unsigned 
    short *flag, int n_obs) {
    
    /* Allocate memory and read data from sourcedata table */
    
    hid_t dataset;
    int status;
    
    dataset = H5Dopen2(group,"/data/sourcedata",H5P_DEFAULT);
            
    status = read_sourcedata_table_obj(objno,dataset,flux,ferr,flag,n_obs);
    
    status = H5Dclose(dataset);
    return status;
}


int getobj (char *filename, unsigned long objno, double **bjd, 
    unsigned short **filterid, float **exptime, unsigned short **pid, 
    unsigned short **field, unsigned short **ccd, unsigned short **quad, 
    unsigned short **imstat, unsigned long **infobits, float **seeing, 
    double **mzpsci, double **mzpsciunc, double **mzpscirms,
    double **clrco, double **clrcounc, float **maglim, float **airmass,
    unsigned long **nps1matches, float **flux, float **ferr, 
    unsigned short **flag, int *n_obs, float *mag, unsigned short *typ, 
    float *magzp) {
    
    hid_t file,group,attr;
    int i,success,n_stars;
    double *jd;
    float *ra,*dec,*mag_ref,*mag_ref_err;
    unsigned short *objtype;
    unsigned long *gaia_id;
    
    file = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr,"Cannot read %s\n",filename); 
        return(0);
    }
    
    group = H5Gopen2(file,"/data",H5P_DEFAULT);
    
    get_exposures(group,&jd,bjd,filterid,exptime,pid,field,ccd,quad,imstat,
        infobits,seeing,mzpsci,mzpsciunc,mzpscirms,clrco,clrcounc,maglim,
        airmass,nps1matches,n_obs);
    get_sources(group,&gaia_id,&ra,&dec,&mag_ref,&mag_ref_err,&objtype,
        &n_stars);
    
    // get magnitude zero-point
    attr = H5Aopen(group,"magzp",H5P_DEFAULT);
    H5Aread(attr,H5T_NATIVE_FLOAT,magzp);
    H5Aclose(attr);
    
    success = -1;
    for (i=0; i<n_stars; i++) {
        if (gaia_id[i] == objno) {success = 1; break;}
    }
    
    if (success == 1) {
        
        *mag = mag_ref[i];
        *typ = objtype[i];
        
        *flux = (float *) malloc(*n_obs*sizeof(float));
        *ferr = (float *) malloc(*n_obs*sizeof(float));
        *flag = (unsigned short *) malloc(*n_obs*sizeof(unsigned short));
        
        get_sourcedata_obj(i,group,*flux,*ferr,*flag,*n_obs);
    
    } 
    
    H5Gclose(group);
    H5Fclose(file);
    
    free(jd);
    free(ra);
    free(dec);
    free(mag_ref);
    free(mag_ref_err);
    free(objtype);
    free(gaia_id);

    return success;
    
}

int getobj_single (int objno, float *flux, float *ferr, float *flux_all, 
    float *ferr_all, unsigned short *flag, int n_obs) {
    
    int i,idx;
    
    idx = objno*n_obs;
    for (i=0; i<n_obs; i++) {
        flux[i] = flux_all[idx];
        ferr[i] = ferr_all[idx];
        idx += 1;
    }
    
    return n_obs;
}

int getobj_double (int objno, double *flux, double *ferr, float *flux_all, 
    float *ferr_all, unsigned short *flag, int n_obs) {
    
    int i,idx;
    
    idx = objno*n_obs;
    for (i=0; i<n_obs; i++) {
        flux[i] = (double)flux_all[idx];
        ferr[i] = (double)ferr_all[idx];
        idx += 1;
    }
    
    return n_obs;
}
