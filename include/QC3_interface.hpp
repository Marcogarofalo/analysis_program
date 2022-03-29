#ifndef PYHELPER_HPP
#define PYHELPER_HPP
// #pragma once

#include <Python.h>
#include <string>
class CPyInstance {
public:
    CPyInstance() {
        Py_Initialize();
    }

    ~CPyInstance() {
        //         Py_Finalize();
        Py_FinalizeEx();
    }
};


class CPyObject {
private:
    PyObject* p;
public:
    CPyObject() : p(NULL) {
    }

    CPyObject(PyObject* _p) : p(_p) {
    }


    ~CPyObject() {
        Release();
    }

    PyObject* getObject() {
        return p;
    }

    PyObject* setObject(PyObject* _p) {
        return (p = _p);
    }

    PyObject* AddRef() {
        if (p) {
            Py_XINCREF(p);
        }
        return p;
    }

    void Release() {
        if (p) {
            Py_XDECREF(p);
        }

        p = NULL;
    }

    PyObject* operator ->() {
        return p;
    }

    bool is() {
        return p ? true : false;
    }

    operator PyObject* () {
        return p;
    }

    PyObject* operator = (PyObject* pp) {
        p = pp;
        return p;
    }

    operator bool() {
        return p ? true : false;
    }
};


double test_python(double a);
double python_detQC(double Estart, double Eend, int steps, double L,
    double* nnP, int Nkcot, double* kcot_params, int Nkiso, double* kiso_params);
void init_python_detQC();
void init_python_detQC_kcot_kiso(std::string kcot, std::string kiso, std::string solve_func = "find_sol");
double python_detQC_call(double Estart, double Eend, int steps,
     double L, double* nnP, int Nkcot, double* kcot_params, int Nkiso, double* kiso_params);
void python_detQC_write_database();
void python_detQC_free();

#endif
