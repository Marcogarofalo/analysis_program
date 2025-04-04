#define QC3_interface_C

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "QC3_interface.hpp"
#include "abstract.h"
#include "floatobject.h"
#include "import.h"
#include "longobject.h"
#include "object.h"
#include "pyerrors.h"
#include "pythonrun.h"
#include "tupleobject.h"
#include "unicodeobject.h"

double test_python(double a ){
    static PyObject *pName, *pModule;
    static PyObject *pFunc;
    static PyObject *pArgs, *pValue;
    if(!pName){
         printf("importing\n");
           PyRun_SimpleString("import sys");
// //          PyRun_SimpleString("sys.path.insert(0, \"/home/marco/analysis/analysis_program/external/qc3/QC3_swave/\")");
//           PyRun_SimpleString("sys.path.append(\"/home/marco/analysis/analysis_program/build/python_wrapper\")");
            PyRun_SimpleString("sys.path.insert(0, \".\")");
        pName = PyUnicode_DecodeFSDefault("python_function");
        if(!pModule){
            pModule = PyImport_Import(pName);
            //         Py_XDECREF(pName);
        }
    }
    
    
    if(pModule)
    {
        if (! pFunc) pFunc = PyObject_GetAttrString(pModule, "prova");
        /* pFunc is a new reference */
            
        if(pFunc && PyCallable_Check(pFunc))
        {
            
            if(!pArgs) pArgs = PyTuple_New(1);
            PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble(a));
            pValue = PyObject_CallObject(pFunc,pArgs );
            //             Py_XDECREF(pArgs);
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \n");
        }
//         Py_XDECREF(pFunc);
//         Py_XDECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load Module\n");
        return -1;
    }
//     if (Py_FinalizeEx() < 0) {
//         return 120;
//     }
    double r=PyFloat_AsDouble(pValue);
    
    return r;
}


static PyObject *pName, *pModule;
static PyObject *pFunc;
static PyObject *pprint;
static PyObject *pkcot, *pkiso, *pkiso1;
static PyObject *np_array;
static PyObject *pArgs, *pArgs1,  *pkcot_params, *pkiso_params ;
static PyObject *pValue;
static PyObject *pwrite_db_Fmat00;

double python_detQC(  double Estart,double Eend,int steps, double L, double *nnP, int Nkcot , double *kcot_params, int Nkiso ,double *kiso_params ){
//     printf(" python is init: %d\n",Py_IsInitialized());
//     CPyInstance hInstance;
//     printf(" python is init: %d\n",Py_IsInitialized());
   
//     static CPyObject pArgs, pArgs1, pArgs2, pkcot_params, pkiso_params ;
    if(!pName){
        PyRun_SimpleString("import sys");
        // PyRun_SimpleString("sys.path.insert(0, \"/home/marco/analysis/analysis_program/external/qc3/QC3_swave/\")");
        PyRun_SimpleString("sys.path.insert(0, \"../../external/qc3/QC3_swave/\")");
        PyRun_SimpleString("sys.path.append(\"../../external/qc3/QC3_swave/F3\")");
        pName = PyUnicode_FromString("detQC");
        pName = PyUnicode_DecodeFSDefault("detQC");
        printf("importing file\n");
        if(!pModule){
            pModule = PyImport_Import(pName);
            printf("importing module\n");
            //         Py_XDECREF(pName);
        }
        
    }
    
    if(pModule)
    {
        //CPyObject pFunc = PyObject_GetAttrString(pModule, "find_sol");
        if(!pFunc){ pFunc = PyObject_GetAttrString(pModule, "find_sol");  printf("importing find_sol\n"); }
        if(!pprint){ pprint = PyObject_GetAttrString(pModule, "print_find_sol");  printf("importing print_find_sol\n"); }
        
        if(!pkcot){ pkcot = PyObject_GetAttrString(pModule, "kcot_2par"); printf("importing kcot_2par\n"); }
        if(! pkcot || !PyCallable_Check(pkcot))  {printf("cannot load kcot function from python\n"); exit(1);}
            
        if(!pkiso1){ pkiso1 = PyObject_GetAttrString(pModule, "kiso_1par"); printf("importing kiso_1par\n"); }  
        if(! pkiso1 ||  !PyCallable_Check(pkiso1))  {printf("cannot load kiso_constant function from python\n"); exit(1);}
        if(!pkiso){ pkiso = PyObject_GetAttrString(pModule, "kiso_2par"); printf("importing kiso_2par\n"); }
        if(! pkiso ||  !PyCallable_Check(pkiso))  {printf("cannot load kiso function from python\n"); exit(1);}
        
        if(!np_array){ np_array = PyObject_GetAttrString(pModule, "get_np_array");  printf("importing get_np_array\n"); }
        if(pFunc && PyCallable_Check(pFunc))
        {
            
            pArgs = PyTuple_New(9); 
//             printf("alloc pArgs\n");
            
            PyTuple_SetItem(pArgs, 0,  PyFloat_FromDouble(Estart)  );
            PyTuple_SetItem(pArgs, 1,  PyFloat_FromDouble(Eend)    );
            PyTuple_SetItem(pArgs, 2,  PyLong_FromLong(steps)      );
            PyTuple_SetItem(pArgs, 3,  PyFloat_FromDouble(L)       );
            
//             if(!pArgs1 )
            pArgs1 = PyTuple_New(3);   
            PyTuple_SetItem(pArgs1, 0, PyLong_FromLong(nnP[0]));
            PyTuple_SetItem(pArgs1, 1, PyLong_FromLong(nnP[1]));
            PyTuple_SetItem(pArgs1, 2, PyLong_FromLong(nnP[2]));
           
            PyTuple_SetItem(pArgs, 4, PyObject_CallObject(np_array,pArgs1 ) );
            
            
            PyTuple_SetItem(pArgs, 5, pkcot);
//             if(!pkcot_params )
            pkcot_params = PyTuple_New(Nkcot);
            for(int i=0;i<Nkcot;i++){
                PyTuple_SetItem(pkcot_params, i, PyFloat_FromDouble( kcot_params[i] )  );
            }
        
            PyTuple_SetItem(pArgs, 6, pkcot_params);
            
            if(Nkiso==1) PyTuple_SetItem(pArgs, 7, pkiso1);
            else if (Nkiso==2) PyTuple_SetItem(pArgs, 7, pkiso);
            else {printf("no kiso function with %d parameters\n",Nkiso) ; exit(1);}
//             if(!pkiso_params )
            pkiso_params = PyTuple_New(Nkiso);
            for(int i=0;i<Nkiso;i++){
                PyTuple_SetItem(pkiso_params, i, PyFloat_FromDouble( kiso_params[i]  ));
            }
            PyTuple_SetItem(pArgs, 8, pkiso_params);
            
            
            
//             PyObject_CallObject(pprint,pArgs );
            pValue = PyObject_CallObject(pFunc,pArgs );
            
    //         printf("C: find_sol() = %g\n", PyFloat_AsDouble(pValue));
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            printf("Cannot find function \n");
        }
        
    }
    else {
        PyErr_Print();
        printf("ERROR python_detQC: Module not imported\n");
        exit(1);
    }
    

    
    return PyFloat_AsDouble(pValue);
}



void init_python_detQC( ){

    if(!pName){
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.insert(0, \"../../external/qc3/QC3_swave/\")");
        PyRun_SimpleString("sys.path.append(\"../../external/qc3/QC3_swave/F3\")");
        pName = PyUnicode_FromString("detQC");
        pName = PyUnicode_DecodeFSDefault("detQC");
        printf("importing file\n");
        if(!pModule){
            pModule = PyImport_Import(pName);
            printf("importing module\n");
            //         Py_XDECREF(pName);
        }
        
    }
    
    if(pModule)
    {
        //CPyObject pFunc = PyObject_GetAttrString(pModule, "find_sol");
        if(!pFunc){ pFunc = PyObject_GetAttrString(pModule, "find_sol");  printf("importing find_sol\n"); }
        if(!pprint){ pprint = PyObject_GetAttrString(pModule, "print_find_sol");  printf("importing print_find_sol\n"); }
        
//         if(!pkcot){ pkcot = PyObject_GetAttrString(pModule, "kcot_2par"); printf("importing kcot_2par\n"); }
//         if(! pkcot || !PyCallable_Check(pkcot))  {printf("cannot load kcot function from python\n"); exit(1);}
            
//         if(!pkiso1){ pkiso1 = PyObject_GetAttrString(pModule, "kiso_1par"); printf("importing kiso_1par\n"); }  
//         if(! pkiso1 ||  !PyCallable_Check(pkiso1))  {printf("cannot load kiso_constant function from python\n"); exit(1);}
//         if(!pkiso){ pkiso = PyObject_GetAttrString(pModule, "kiso_2par"); printf("importing kiso_2par\n"); }
//         if(! pkiso ||  !PyCallable_Check(pkiso))  {printf("cannot load kiso function from python\n"); exit(1);}
        
        if(!np_array){ np_array = PyObject_GetAttrString(pModule, "get_np_array");  printf("importing get_np_array\n"); }
        if(pFunc && PyCallable_Check(pFunc))
        {
            printf("python function callable\n");
            
    //         printf("C: find_sol() = %g\n", PyFloat_AsDouble(pValue));
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            printf("Cannot find function \n");
        }
        
    }
    else {
        PyErr_Print();
        printf("ERROR python_detQC: Module not imported\n");
        exit(1);
    }
}


void init_python_detQC_kcot_kiso(std::string kcot , std::string kiso, std::string solve_func){

    
    if(pModule)
    {
        pFunc = PyObject_GetAttrString(pModule, solve_func.c_str() );  printf("importing %s\n",kcot.c_str()); 
        if(! pFunc || !PyCallable_Check(pFunc))  {printf("cannot load %s function from python\n",solve_func.c_str()); exit(1);}
          
        pkcot = PyObject_GetAttrString(pModule, kcot.c_str() ); printf("importing %s\n",kcot.c_str()); 
        if(! pkcot || !PyCallable_Check(pkcot))  {printf("cannot load %s function from python\n",kcot.c_str()); exit(1);}

        pkiso = PyObject_GetAttrString(pModule, kiso.c_str() ); printf("importing %s\n",kiso.c_str() ); 
        if(! pkiso || !PyCallable_Check(pkiso))  {printf("cannot load %s function from python\n",kiso.c_str()); exit(1);}
    
    }
    else {
        PyErr_Print();
        printf("ERROR python_detQC: Module not imported\n");
        exit(1);
    }
}



double python_detQC_call(  double Estart,double Eend,int steps, double L, double *nnP, int Nkcot , double *kcot_params, int Nkiso ,double *kiso_params ){
    if(pModule)
    {
        if(pFunc && PyCallable_Check(pFunc))
        {
            
            pArgs = PyTuple_New(9); 
//             printf("alloc pArgs\n");
            
            PyTuple_SetItem(pArgs, 0,  PyFloat_FromDouble(Estart)  );
            PyTuple_SetItem(pArgs, 1,  PyFloat_FromDouble(Eend)    );
            PyTuple_SetItem(pArgs, 2,  PyLong_FromLong(steps)      );
            PyTuple_SetItem(pArgs, 3,  PyFloat_FromDouble(L)       );
            
//             if(!pArgs1 )
            pArgs1 = PyTuple_New(3);   
            PyTuple_SetItem(pArgs1, 0, PyLong_FromLong(nnP[0]));
            PyTuple_SetItem(pArgs1, 1, PyLong_FromLong(nnP[1]));
            PyTuple_SetItem(pArgs1, 2, PyLong_FromLong(nnP[2]));
           
            PyTuple_SetItem(pArgs, 4, PyObject_CallObject(np_array,pArgs1 ) );
            
            
            PyTuple_SetItem(pArgs, 5, pkcot);
            pkcot_params = PyTuple_New(Nkcot);
            for(int i=0;i<Nkcot;i++){
                PyTuple_SetItem(pkcot_params, i, PyFloat_FromDouble( kcot_params[i] )  );
            }
            PyTuple_SetItem(pArgs, 6, pkcot_params);
            
            PyTuple_SetItem(pArgs, 7, pkiso);
            pkiso_params = PyTuple_New(Nkiso);
            for(int i=0;i<Nkiso;i++){
                PyTuple_SetItem(pkiso_params, i, PyFloat_FromDouble( kiso_params[i]  ));
            }
            PyTuple_SetItem(pArgs, 8, pkiso_params);
            
            
            
//             PyObject_CallObject(pprint,pArgs );
            pValue = PyObject_CallObject(pFunc,pArgs );
            
    //         printf("C: find_sol() = %g\n", PyFloat_AsDouble(pValue));
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            printf("Cannot find function \n");
        }
        
    }
    else {
        PyErr_Print();
        printf("ERROR python_detQC: Module not imported\n");
        exit(1);
    }
    

    
    return PyFloat_AsDouble(pValue);
}




void python_detQC_write_database(){
    if(pModule)
    {
        if(!pwrite_db_Fmat00)
            pwrite_db_Fmat00 = PyObject_GetAttrString(pModule, "write_db_Fmat00");   
        if(pwrite_db_Fmat00 && PyCallable_Check(pwrite_db_Fmat00))
            PyObject_CallObject(pwrite_db_Fmat00, PyTuple_New(0) );
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            printf("Cannot find function write_db_Fmat00\n");
        }
    
    }
    else {
        PyErr_Print();
        printf("ERROR python_detQC: Module not imported\n");
        exit(1);
    }
    
}
    
    


void python_detQC_free(){
    Py_XDECREF(pName);
    Py_XDECREF(pModule);
    Py_XDECREF(pFunc);
     
    Py_XDECREF(pprint);
    Py_XDECREF(pkcot);
    Py_XDECREF(pkiso);
    Py_XDECREF(pkiso1);
    
    Py_XDECREF(np_array);
    Py_XDECREF(pArgs);
    Py_XDECREF(pArgs1);
    Py_XDECREF(pkcot_params);
    
    Py_XDECREF(pkiso_params);
    Py_XDECREF(pValue);
    Py_XDECREF(pwrite_db_Fmat00);
    
    
}
