#include "Windows.h"
#include "astercxx.h"

#include "aster_init.h"
#include "aster_numpy.h"
#include "aster_pybind.h"


namespace py = pybind11;

// Function pointer type for the original initialization function
typedef PyObject* (*PyInitFunc)();

// Proxy function to initialize the module
extern "C" __declspec(dllexport) PyObject* PyInit_libaster() {
    static HMODULE originalModule = LoadLibrary("bibcxx.dll");

    if (!originalModule) {
        PyErr_SetString(PyExc_ImportError, "Could not load original module");
        return nullptr;
    }

    static PyInitFunc pyInitFunc = (PyInitFunc)GetProcAddress(originalModule, "PyInit_libaster");
    if (!pyInitFunc) {
        PyErr_SetString(PyExc_ImportError, "Could not find PyInit_libaster in original module");
        return nullptr;
    }

    return pyInitFunc();
}
