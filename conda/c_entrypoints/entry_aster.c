//
// Created by Kristoffer on 14.06.2024.
//
// This file is part of Code_Aster.

#include "aster_module.h"

/* Initialization function for the module (*must* be called initaster) */
static char aster_module_documentation[] = "C implementation of the Python aster module\n"
                                           "\n";

static struct PyModuleDef aster_def = { PyModuleDef_HEAD_INIT,
                                        "aster",
                                        aster_module_documentation,
                                        -1,
                                        aster_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL };

PyObject *PyInit_aster( void ) {
    PyObject *aster = (PyObject *)0;

    /* Create the module and add the functions */
    aster = PyModule_Create( &aster_def );

    init_etape_stack();
    return aster;
}