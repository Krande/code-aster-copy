/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#define CODEASTER_IMPORT_ARRAY 1
#include "aster_numpy.h"

extern void calc_SPEC_OSCI( int, double *, double *, int, double *, int, double *, double * );

static PyObject *SPEC_OSCI( PyObject *self, PyObject *args ) {
    PyObject *OX, *OY, *OF, *OA;
    PyArrayObject *Vx, *Vy, *Vf, *Va, *Sp;
    npy_intp dims[3];
    int nbpts, len_f, len_a;

    if ( !PyArg_ParseTuple( args, "OOOO:SPEC_OSCI", &OX, &OY, &OF, &OA ) )
        return NULL;

    Vx = (PyArrayObject *)PyArray_FROM_OTF( OX, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
    Vy = (PyArrayObject *)PyArray_FROM_OTF( OY, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
    Vf = (PyArrayObject *)PyArray_FROM_OTF( OF, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
    Va = (PyArrayObject *)PyArray_FROM_OTF( OA, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );

    /* verification des arguments */
    if ( !( Vx && Vy && Va && Vf ) )
        return NULL;
    if ( !( PyArray_NDIM( Vx ) == 1 && PyArray_NDIM( Vy ) == 1 && PyArray_NDIM( Vf ) == 1 &&
            PyArray_NDIM( Va ) == 1 ) ) {
        PyErr_SetString( PyExc_TypeError, "On attend des objets de dimension 1" );
        return NULL;
    }
    if ( PyArray_DIM( Vx, 0 ) != PyArray_DIM( Vy, 0 ) ) {
        PyErr_SetString( PyExc_TypeError, "Vx et Vy n'ont pas le meme cardinal" );
        return NULL;
    }
    if ( PyArray_DIM( Vx, 0 ) <= 1 ) {
        PyErr_SetString( PyExc_TypeError, "Vx et Vy n'ont pas assez de valeurs" );
        return NULL;
    }

    nbpts = PyArray_DIM( Vx, 0 );
    len_f = PyArray_DIM( Vf, 0 );
    len_a = PyArray_DIM( Va, 0 );

    dims[0] = (npy_intp)len_a;
    dims[1] = 3;
    dims[2] = (npy_intp)len_f;

    Sp = (PyArrayObject *)PyArray_SimpleNew( 3, dims, NPY_DOUBLE );
    calc_SPEC_OSCI( nbpts, (double *)PyArray_DATA( Vx ), (double *)PyArray_DATA( Vy ), len_f,
                    (double *)PyArray_DATA( Vf ), len_a, (double *)PyArray_DATA( Va ),
                    (double *)PyArray_DATA( Sp ) );

    Py_DECREF( Vx );
    Py_DECREF( Vy );
    Py_DECREF( Vf );
    Py_DECREF( Va );

    return PyArray_Return( Sp );
}

#ifdef ASTER_DEBUG_FONCTIONS
/* utile pour le remplissage des contiguous array */
static PyObject *_INFO( PyObject *self, PyObject *args ) {
    PyObject *OX;
    PyArrayObject *Vx;
    int ndim, i, j, k, sum, ind;
    double *val;

    if ( !PyArg_ParseTuple( args, "O:_INFO", &OX ) )
        return NULL;

    Vx = (PyArrayObject *)PyArray_ContiguousFromAny( OX, NPY_DOUBLE, 3, 3 );

    if ( !Vx )
        return NULL;

    ndim = PyArray_NDIM( Vx );
    printf( "Tableau de dimension : %d\n", ndim );
    sum = 1;
    for ( i = 0; i < ndim; i++ ) {
        printf( "dimensions[%d]=%d        strides[%d]=%d\n", i, PyArray_DIM( Vx, i ), i,
                PyArray_STRIDE( Vx, i ) );
        sum = sum * PyArray_DIM( Vx, i );
    }

    if ( ndim == 3 ) {
        printf( "DUMP du tableau :\n-----------------\n" );
        for ( i = 0; i < PyArray_DIM( Vx, 0 ); i++ ) {
            for ( j = 0; j < PyArray_DIM( Vx, 1 ); j++ ) {
                for ( k = 0; k < PyArray_DIM( Vx, 2 ); k++ ) {
                    val = (double *)( PyArray_DATA( Vx ) + i * PyArray_STRIDE( Vx, 0 ) +
                                      j * PyArray_STRIDE( Vx, 1 ) + k * PyArray_STRIDE( Vx, 2 ) );
                    printf( "Vx[%d,%d,%d]=%lf\n", i, j, k, *val );
                }
            }
        }
        printf( "\nDUMP du data :\n--------------\n" );
        val = PyArray_DATA( Vx );
        for ( i = 0; i < PyArray_DIM( Vx, 0 ); i++ ) {
            for ( j = 0; j < PyArray_DIM( Vx, 1 ); j++ ) {
                for ( k = 0; k < PyArray_DIM( Vx, 2 ); k++ ) {
                    ind = i * PyArray_DIM( Vx, 1 ) * PyArray_DIM( Vx, 2 ) +
                          j * PyArray_DIM( Vx, 2 ) + k;
                    printf( "Vx[%d,%d,%d / %d]=%lf\n", i, j, k, ind, val[ind] );
                }
            }
        }
    }

    Py_DECREF( Vx );
    Py_INCREF( Py_None );
    return Py_None;
}
#endif

#ifndef ASTER_WITHOUT_PYMOD
static PyMethodDef methods[] = {
    {"SPEC_OSCI", SPEC_OSCI, METH_VARARGS, "Operation SPEC_OSCI de CALC_FONCTION"},
#ifdef ASTER_DEBUG_FONCTIONS
    {"_INFO", _INFO, METH_VARARGS, "Just for test !"},
#endif
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef aster_fonctions_def = {
    PyModuleDef_HEAD_INIT, "aster_fonctions", NULL, -1, methods, NULL, NULL, NULL, NULL};

PyObject *PyInit_aster_fonctions( void ) {
    PyObject *aster_fonctions = PyModule_Create( &aster_fonctions_def );
    import_array();
    return aster_fonctions;
}
#endif
