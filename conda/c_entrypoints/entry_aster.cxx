#include "entry_aster.h"


namespace py = pybind11;

// Function pointer type for the original initialization function
typedef PyObject* (*PyInitFunc)();

// Global variable to store the original module's handle and function pointer
static HMODULE originalModule = nullptr;
static PyInitFunc pyInitFunc = nullptr;

// Proxy function to initialize the module
extern "C" __declspec(dllexport) PyObject* PyInit_aster() {
    if (!originalModule) {
        originalModule = LoadLibrary("bibc.dll");
        if (!originalModule) {
            PyErr_SetString(PyExc_ImportError, "Could not load original module");
            return nullptr;
        }

        pyInitFunc = (PyInitFunc)GetProcAddress(originalModule, "PyInit_aster");
        if (!pyInitFunc) {
            PyErr_SetString(PyExc_ImportError, "Could not find PyInit_aster in original module");
            return nullptr;
        }
    }

    return pyInitFunc();
}
