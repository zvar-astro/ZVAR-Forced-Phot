/* Codes for reading data from forced-phot ZTF matchfiles. 
 * 
 * For reading data from HDF5 files, see:
 * http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html
 * 
 * P. Mroz @ Caltech, 21 Dec 2020
 *
 * Modified by Zach Vanderbosch @ Caltech, 2022-2024
 */
 
#ifndef MATCHFILE_PM_H_
#define MATCHFILE_PM_H_
 
typedef struct sources {
    unsigned long gaia_id;
    float ra;
    float dec;
    float mag_ref;
    float mag_err_ref;
    unsigned short objtype;
} sources;
    
typedef struct exposures {
    double jd;
    double bjd;
    unsigned short filterid;
    float exptime;
    unsigned short pid;
    unsigned short field;
    unsigned short ccd;
    unsigned short quad;
    unsigned short imstat;
    unsigned long infobits;
    float seeing;
    double mzpsci;
    double mzpsciunc;
    double mzpscirms;
    double clrco;
    double clrcounc;
    float maglim;
    float airmass;
    unsigned long nps1matches;
} exposures;

typedef struct sourcedata {
    float flux;
    float flux_err;
    unsigned short flag;
} sourcedata;

int get_number_of_rows(hid_t);

int read_exposures_table(hid_t, double*, double*, unsigned short*, float*, 
    unsigned short*, unsigned short*, unsigned short*, unsigned short*, 
    unsigned short*, unsigned long*, float*, double*, double*, double*,
    double*, double*, float*, float*, unsigned long*, int);
int read_sources_table(hid_t, unsigned long*, float*, float*, float*, 
    float*, unsigned short*, int);
int read_sourcedata_table(hid_t, float*, float*, unsigned short*, int);

int get_exposures (hid_t, double**, double**, unsigned short**, float**, 
    unsigned short**, unsigned short**, unsigned short**, unsigned short**, 
    unsigned short**, unsigned long**, float**, double**, double**, double**,
    double**, double**, float**, float**, unsigned long**, int*);
int get_sources (hid_t, unsigned long**, float**, float**, float**, 
    float**, unsigned short**, int*);
int get_sourcedata (hid_t, float**, float**, unsigned short**, int*);
    
int getobj_single (int, float*, float*, float*, float*, unsigned short*,
    int);
int getobj_double (int, double*, double*, float*, float*, unsigned short*,
    int);

int getobj (char*, unsigned long, double**, unsigned short**, float**, 
    unsigned short**, unsigned short**, unsigned short**, unsigned short**, 
    unsigned short**, unsigned long**, float**, double**, double**, double**,
    double**, double**, float**, float**, unsigned long**, float**, 
    float**, unsigned short**, int*, float*, unsigned short*, float*);

#endif
