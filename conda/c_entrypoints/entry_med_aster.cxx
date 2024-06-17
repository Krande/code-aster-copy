#include "Windows.h"
#include "astercxx.h"

#include "aster_init.h"
#include "aster_numpy.h"
#include "aster_pybind.h"


namespace py = pybind11;

// Function pointer type for the original initialization function
typedef PyObject* (*PyInitFunc)();

// Load the original DLL and get the initialization function pointer
extern "C" __declspec(dllexport) PyObject* PyInit_aster_core() {
    static HMODULE originalModule = LoadLibrary("bibc.dll");
    if (!originalModule) {
        PyErr_SetString(PyExc_ImportError, "Could not load original module");
        return nullptr;
    }

    static PyInitFunc pyInitFunc = (PyInitFunc)GetProcAddress(originalModule, "PyInit_med_aster");
    if (!pyInitFunc) {
        PyErr_SetString(PyExc_ImportError, "Could not find PyInit_med_aster in original module");
        return nullptr;
    }

    return pyInitFunc();
}
