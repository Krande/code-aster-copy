#include <Python.h>
#include "Windows.h"
#include "astercxx.h"

#include "aster_init.h"
#include "aster_numpy.h"
#include "aster_pybind.h"


// wrapper.c


typedef PyObject* (*PyInit_aster_t)(void);

// Pointer to the original PyInit_aster function
PyInit_aster_t original_PyInit_aster = NULL;

// Function to load the original bibc.dll and get the PyInit_aster function
void load_bibc_dll() {
    HMODULE hBibcDll = LoadLibrary("bibc.dll");
    if (!hBibcDll) {
        PyErr_SetString(PyExc_ImportError, "Could not load bibc.dll");
        return;
    }

    original_PyInit_aster = (PyInit_aster_t)GetProcAddress(hBibcDll, "PyInit_aster");
    if (!original_PyInit_aster) {
        PyErr_SetString(PyExc_ImportError, "Could not find PyInit_aster in bibc.dll");
        return;
    }
}

// The wrapper function that calls the original PyInit_aster
PyObject* PyInit_aster() {
    if (!original_PyInit_aster) {
        load_bibc_dll();
    }

    if (original_PyInit_aster) {
        return original_PyInit_aster();
    }

    return NULL;
}

