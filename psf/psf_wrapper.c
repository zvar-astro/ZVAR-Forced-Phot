#include <Python.h>
#include <numpy/arrayobject.h>

void phot2 (double *flux, double *ferr, char *diff_imagename, 
    char *sci_imagename, char *psf_imagename, double *x, double *y, 
    double gain, int nstars, int stamp);
    
static PyObject * 
get_phot_wrapper (PyObject *self, PyObject *args) {
    
    int n,stamp;
    char *str1,*str2,*str3;
    double *x,*y,*flux,*ferr,gain;
    PyObject *x_obj,*y_obj,*flux_obj,*ferr_obj;
    
    if (!PyArg_ParseTuple(args,"sssOOdii",\
        &str1,&str2,&str3,&x_obj,&y_obj,&gain,&n,&stamp)) return NULL;
        
    PyObject *x_obj_ = PyArray_FROM_OTF(x_obj,NPY_DOUBLE,NPY_IN_ARRAY);
    PyObject *y_obj_ = PyArray_FROM_OTF(y_obj,NPY_DOUBLE,NPY_IN_ARRAY);
    
    if (x_obj_ == NULL || y_obj_ == NULL) {
        Py_XDECREF(x_obj_);
        Py_XDECREF(y_obj_);
        return NULL;
    }
    
    x = (double *) PyArray_DATA(x_obj_);
    y = (double *) PyArray_DATA(y_obj_);
    
    flux = (double *) PyMem_Malloc(n*sizeof(double));
    ferr = (double *) PyMem_Malloc(n*sizeof(double));

    phot2(flux,ferr,str1,str2,str3,x,y,gain,n,stamp);
    
    Py_DECREF(x_obj_);
    Py_DECREF(y_obj_);
    
    npy_intp dims = {n};
    
    flux_obj = PyArray_SimpleNew(1,&dims,NPY_DOUBLE);
    memcpy(PyArray_DATA(flux_obj),flux,n*sizeof(double));
    ferr_obj = PyArray_SimpleNew(1,&dims,NPY_DOUBLE);
    memcpy(PyArray_DATA(ferr_obj),ferr,n*sizeof(double));
    
    PyMem_Free(flux);
    PyMem_Free(ferr);
    
    return Py_BuildValue("NN",flux_obj,ferr_obj);
}

static PyMethodDef PSF_Phot_PMMethods[] = {
    {"get_phot", (PyCFunction) get_phot_wrapper, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef PSF_Phot_PM = {
    PyModuleDef_HEAD_INIT,
    "PSF_Phot_PM",
    "usage",
    -1,
    PSF_Phot_PMMethods
};

PyMODINIT_FUNC PyInit_PSF_Phot_PM(void) {
    import_array(); 
    // you need add this line to prevent segmentation faults
    // see: http://jamesmontgomery.us/blog/Extensions_in_C.html
    return PyModule_Create(&PSF_Phot_PM);
}


