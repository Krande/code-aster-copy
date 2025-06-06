/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
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

/* person_in_charge: mathieu.courtois at edf.fr */

#include "aster_utils.h"

#include "aster_fort_utils.h"
#include "aster_module.h"

/*
 * Fonctions de manipulation des chaines de caractères pour échange avec le Fortran.
 */
STRING_SIZE FStrlen( _IN const char *fstr, _IN const STRING_SIZE flen ) {
    /* Retourne la longueur (dernier caractère non blanc) de "fstr".
     */
    STRING_SIZE n;
    _check_string_length( flen );
    n = flen;
    while ( n >= 1 && fstr[n - 1] == ' ' ) {
        n--;
    }
    return n;
}

char *MakeCStrFromFStr( _IN const char *fstr, _IN const STRING_SIZE flen ) {
    /* Alloue et retourne une chaine C (terminant par \0) étant
     * la copie de la chaine Fortran sans les blancs finaux.
     * La chaine devra etre libérée par l'appelant.
     */
    char *cstr = NULL;
    STRING_SIZE n;

    n = FStrlen( fstr, flen );
    cstr = (char *)malloc( ( n + 1 ) * sizeof( char ) );
    strncpy( cstr, fstr, n );
    cstr[n] = '\0';

    return cstr;
}

void CopyCStrToFStr( _INOUT char *fstr, _IN const char *cstr, _IN const STRING_SIZE flen ) {
    /* Copie une chaine C dans une chaine Fortran déjà allouée (de taille
     * flen) et sans ajout du '\0' à la fin.
     */
    STRING_SIZE i, n;
    n = strlen( cstr );
    if ( n > flen ) {
        n = flen;
    }
    for ( i = 0; i < n; i++ ) {
        fstr[i] = cstr[i];
    }
    while ( i < flen ) {
        fstr[i] = ' ';
        i++;
    }
}

char *MakeFStrFromCStr( _IN const char *cstr, _IN const STRING_SIZE flen ) {
    /* Alloue et retourne une chaine C (complétée par des blancs
     * destinée à être transmise au Fortran, d'où FStr) étant
     * la copie de la chaine C.
     * La chaine devra etre libérée par l'appelant.
     */
    char *fstr = NULL;
    fstr = (char *)malloc( ( flen + 1 ) * sizeof( char ) );
    CopyCStrToFStr( fstr, cstr, flen );
    fstr[flen] = '\0';
    return fstr;
}

void BlankStr( _IN char *fstr, _IN const STRING_SIZE flen ) {
    /* Initialise un blanc une chaine de caractères (sans '\0' à la fin).
     * S'applique à une chaine allouée par le Fortran.
     */
    memset( fstr, ' ', flen );
}

char *MakeBlankFStr( _IN const STRING_SIZE flen ) {
    /* Initialise un blanc une chaine de caractères avec '\0' à la fin
     * (qui peut ainsi être passé au Fortran).
     * Alloue une chaine qui sera passée au Fortran.
     */
    char *fstr;
    fstr = (char *)malloc( ( flen + 1 ) * sizeof( char ) );
    BlankStr( fstr, flen );
    fstr[flen] = '\0';
    return fstr;
}

char *MakeTabFStr( _IN const int size, _IN const STRING_SIZE flen ) {
    /* Alloue un tableau de chaine de caractères Fortran. Chaque chaine
     * est de longueur "flen". Le même "flen" sera utilisé
     * dans SetTabFStr.
     * Alloue un tableau de chaine qui sera passé au Fortran.
     */
    return MakeBlankFStr( size * flen );
}

void SetTabFStr( _IN char *tab, _IN const int index, _IN const char *cstr,
                 _IN const STRING_SIZE flen ) {
    /* Remplit l'indice "index" (de 0 à size-1) du tableau de chaine
     * de caractères "tab" avec la chaine "cstr".
     */
    char *strk = NULL;
    strk = &tab[index * flen];
    CopyCStrToFStr( strk, cstr, flen );
}

/* pour que ce soit clair */
void FreeStr( char *cstr ) { free( cstr ); }

void _check_string_length( const STRING_SIZE flen ) {
    if ( flen > 2147483647 ) {
        printf( "WARNING: The string length seems corrupted. "
                "The value of STRING_SIZE is probably bad : %s (%d bytes)\n"
                "Please contact your support.\n",
                xstr( STRING_SIZE ), (int)sizeof( STRING_SIZE ) );
    }
}

/*
 * Fonctions de conversion
 */
void convc8( _IN const int nval, _IN PyObject *tup, _OUT ASTERDOUBLE *val ) {
    /*
     * tup est un tuple de tuples internes, chaque tuple
     * interne contenant le type et les deux parties du complexe.
     */
    int i = 0;
    int k = 0;
    int conv_un_c8( _IN PyObject * tup, _OUT ASTERDOUBLE * val );
    if ( nval != 0 ) {
        PyObject *v = (PyObject *)0;
        for ( i = 0; i < nval; i++ ) {
            v = PyTuple_GetItem( tup, i );
            k += conv_un_c8( v, val + k );
        }
    }
    return;
}

int conv_un_c8( _IN PyObject *tup, _OUT ASTERDOUBLE *val ) {
    /*
     * Enrichissement des complexes stockes dans val a partir du tuple tup
     */
    char *repres = (char *)0; /* representation "RI" (reelle/imaginaire) ou "MP" (module phase) */
    double x = 0.0;
    double y = 0.0;
    double *rho = &x;
    double *theta = &y;
    if ( PyComplex_Check( tup ) || PyFloat_Check( tup ) || PyLong_Check( tup ) ) {
        /* On est dans le cas d'un objet Python complexe */
        /* representation : partie reelle/partie imaginaire */
        *val = (ASTERDOUBLE)PyComplex_RealAsDouble( tup );
        if ( PyErr_Occurred() ) {
            PyErr_Print();
            MYABORT( "TypeError: Expecting a PyComplex" );
        }
        *( val + 1 ) = (ASTERDOUBLE)PyComplex_ImagAsDouble( tup );
        if ( PyErr_Occurred() ) {
            PyErr_Print();
            MYABORT( "TypeError: Expecting a PyComplex" );
        }
    } else if ( PyTuple_Check( tup ) ) {
        /* On est dans le cas d'un complexe représenté par un triplet : "RI" ou "MP",x,y */
        if ( !PyArg_ParseTuple( tup, "sdd", &repres, &x, &y ) )
            MYABORT( "erreur dans la partie Python" );
        if ( strcmp( repres, "RI" ) == 0 ) {
            /* representation : partie reelle/partie imaginaire */
            *val = (ASTERDOUBLE)x;
            *( val + 1 ) = (ASTERDOUBLE)y;
        } else {
            /* representation RHO,THETA (les angles sont fournis en degres) */
            *val = (ASTERDOUBLE)( *rho * cos( *theta / 180. * CALL_R8PI() ) );
            *( val + 1 ) = (ASTERDOUBLE)( *rho * sin( *theta / 180. * CALL_R8PI() ) );
        }
    } else {
        MYABORT( "erreur dans la partie Python" );
    }
    return 2;
}

void convr8( _IN const int nval, _IN PyObject *tup, _OUT ASTERDOUBLE *val ) {
    /*
     * Convertit un Tuple en tableau de double
     */
    int i;
    PyObject *v = (PyObject *)0;
    if ( nval == 0 )
        return;
    if ( !PyTuple_Check( tup ) ) {
        printf( "tup : " );
        PyObject_Print( tup, stdout, 0 );
        printf( "\n " );
        MYABORT( "erreur sur le type : devrait etre un tuple" );
    }
    for ( i = 0; i < nval; i++ ) {
        v = PyTuple_GetItem( tup, i );
        val[i] = (ASTERDOUBLE)PyFloat_AsDouble( v );
        if ( PyErr_Occurred() ) {
            PyErr_Print();
            MYABORT( "TypeError: Expecting a PyFloat" );
        }
    }
    return;
}

void convert( _IN const int nval, _IN PyObject *tup, _OUT ASTERINTEGER *val ) {
    /*
     * Convertit un Tuple en tableau d entier
     */
    int i;
    PyObject *v = (PyObject *)0;
    if ( nval == 0 )
        return;
    if ( !PyTuple_Check( tup ) ) {
        printf( "tup : " );
        PyObject_Print( tup, stdout, 0 );
        printf( "\n " );
        MYABORT( "erreur sur le type : devrait etre un tuple" );
    }
    for ( i = 0; i < nval; i++ ) {
        v = PyTuple_GetItem( tup, i );
        val[i] = (ASTERINTEGER)PyLong_AsLong( v );
        if ( PyErr_Occurred() ) {
            PyErr_Print();
            MYABORT( "TypeError: Expecting a PyLong" );
        }
    }
    return;
}

void convertxt( _IN const int nval, _IN PyObject *tup, _OUT char *val, _IN STRING_SIZE taille ) {
    /* Convertit un Tuple en tableau de chaines
     * Pour retour au Fortran : le tableau existe deja (val)
     *  nval   : indique le nombre d'elements du tuple a convertir
     *  tup    : est le tuple Python a convertir
     *  val    : est le tableau de chaines Fortran a remplir
     *  taille : indique la taille des chaines
     */
    PyObject *v = (PyObject *)0;
    int i;
    const char *s;
    if ( nval != 0 ) {
        if ( !PyTuple_Check( tup ) ) {
            printf( "tup : " );
            PyObject_Print( tup, stdout, 0 );
            printf( "\n " );
            MYABORT( "erreur sur le type : devrait etre un tuple" );
        }
        for ( i = 0; i < nval; i++ ) {
            v = PyTuple_GetItem( tup, i );
            /*                               v=PySequence_GetItem(tup,i); */
            s = PyUnicode_AsUTF8( v );

            if ( s == NULL ) {
                printf( "s : " );
                PyObject_Print( v, stdout, 0 );
                printf( "\n " );
                MYABORT( "erreur sur le type : devrait etre une string" );
            }

            /* le fortran attend des chaines de caracteres completees par des blancs */
            SetTabFStr( val, i, s, taille );
        }
    }
}

void converltx( _IN const int nval, _IN PyObject *tup, _OUT char *val, _IN STRING_SIZE taille ) {
    /* Convertit une Liste  en tableau de chaines
     * Pour retour au Fortran : le tableau existe deja (val)
     */
    PyObject *v = (PyObject *)0;
    int i;
    const char *s = (char *)0;

    if ( nval != 0 ) {
        if ( !PyList_Check( tup ) ) {
            printf( "tup : " );
            PyObject_Print( tup, stdout, 0 );
            printf( "\n " );
            MYABORT( "erreur sur le type : devrait etre une liste" );
        }
        for ( i = 0; i < nval; i++ ) {
            v = PyList_GetItem( tup, i );
            /* v=PySequence_GetItem(tup,i); */
            s = PyUnicode_AsUTF8( v );
            if ( s == NULL ) {
                printf( "s : " );
                PyObject_Print( v, stdout, 0 );
                printf( "\n " );
                MYABORT( "erreur sur le type : devrait etre une string" );
            }

            /* le fortran attend des chaines de caracteres completees par des blancs */
            SetTabFStr( val, i, s, taille );
        }
    }
    return;
}

/*
 * Fonctions pour créer des listes et tuples
 */
PyObject *MakeTupleString( long nbval, char *kval, STRING_SIZE lkval, ASTERINTEGER *lval ) {
    /*
     *   Entrees:
     *      nbval nombre de chaines dans kval
     *      kval  tableau de nbval chaines FORTRAN
     *      lkval longueur des chaines FORTRAN (compilateur)
     *      lval  longueur des nbval chaines FORTRAN (utilisateur)
     *   Sorties:
     *      RETOUR fonction : tuple de string Python de longueur nbval
     *   Fonction:
     *      Convertir un tableau de chaines FORTRAN en un tuple de string Python de meme
     * longueur
     */
    int i;
    int len;
    char *deb = kval;
    PyObject *tupl;
    tupl = PyTuple_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( lval ) {
            len = (int)lval[i];
        } else {
            len = lkval;
        }
        if ( PyTuple_SetItem( tupl, i, PyUnicode_FromStringAndSize( deb, FStrlen( deb, len ) ) ) ) {
            Py_DECREF( tupl );
            return NULL;
        }
        deb = deb + lkval;
    }
    return tupl;
}

PyObject *MakeListString( long nbval, char *kval, STRING_SIZE lkval ) {
    /*
     *   Entrees:
     *      nbval nombre de chaines dans kval
     *      kval  tableau de nbval chaines FORTRAN
     *      lkval longueur des chaines FORTRAN (compilateur)
     *   Sorties:
     *      RETOUR fonction : tuple de string Python de longueur nbval les espaces terminant la
     *      chaine sont supprimes
     *   Fonction:
     *      Convertir un tableau de chaines FORTRAN en un tuple de string Python de meme
     * longueur
     */
    int i;
    char *deb = kval;
    PyObject *l = PyList_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( PyList_SetItem( l, i, PyUnicode_FromStringAndSize( deb, FStrlen( deb, lkval ) ) ) ) {
            Py_DECREF( l );
            return NULL;
        }
        deb = deb + lkval;
    }
    return l;
}

PyObject *MakeTupleInt( long nbval, ASTERINTEGER *kval ) {
    /*
     *   Entrees:
     *      nbval nombre d'entiers dans kval
     *      kval  tableau de nbval ASTERINTEGER FORTRAN
     *   Sorties:
     *      RETOUR fonction : tuple de int Python de longueur nbval
     *   Fonction:
     *      Convertir un tableau de ASTERINTEGER FORTRAN en un tuple de int Python de meme
     * longueur
     */
    int i;
    PyObject *tupl;
    tupl = PyTuple_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( PyTuple_SetItem( tupl, i, PyLong_FromLong( (long)kval[i] ) ) ) {
            Py_DECREF( tupl );
            return NULL;
        }
    }
    return tupl;
}

PyObject *MakeListInt( long nbval, ASTERINTEGER *kval ) {
    /*
     *   Entrees:
     *       nbval nombre d'entiers dans kval
     *       kval  tableau de nbval ASTERINTEGER FORTRAN
     *   Sorties:
     *       RETOUR fonction : liste de int Python de longueur nbval
     *   Fonction:
     *       Convertir un tableau de ASTERINTEGER FORTRAN en une liste de int Python de meme
     * longueur
     */
    int i;
    PyObject *l = PyList_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( PyList_SetItem( l, i, PyLong_FromLong( (long)kval[i] ) ) ) {
            Py_DECREF( l );
            return NULL;
        }
    }
    return l;
}

PyObject *MakeTupleFloat( long nbval, ASTERDOUBLE *kval ) {
    /*
     *  Entrees:
     *     nbval nombre de reels dans kval
     *     kval  tableau de nbval double FORTRAN
     *  Sorties:
     *     RETOUR fonction : tuple de float Python de longueur nbval
     *  Fonction:
     *     Convertir un tableau de double FORTRAN en un tuple de float Python de meme longueur
     */
    int i;
    PyObject *tupl;
    tupl = PyTuple_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( PyTuple_SetItem( tupl, i, PyFloat_FromDouble( (double)kval[i] ) ) ) {
            Py_DECREF( tupl );
            return NULL;
        }
    }
    return tupl;
}

PyObject *MakeListFloat( long nbval, ASTERDOUBLE *kval ) {
    /*
     * Entrees:
     *      nbval nombre de reels dans kval
     *      kval  tableau de nbval double FORTRAN
     * Sorties:
     *  RETOUR fonction : list de float Python de longueur nbval
     * Fonction:
     *      Convertir un tableau de double FORTRAN en une liste de float Python de meme longueur
     */
    int i;
    PyObject *l = PyTuple_New( (Py_ssize_t)nbval );
    for ( i = 0; i < nbval; i++ ) {
        if ( PyList_SetItem( l, i, PyFloat_FromDouble( (double)kval[i] ) ) ) {
            Py_DECREF( l );
            return NULL;
        }
    }
    return l;
}
