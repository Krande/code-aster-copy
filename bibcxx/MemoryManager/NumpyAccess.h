#ifndef NUMPYACCESS_H_
#define NUMPYACCESS_H_

/**
 * @file NumpyAccess.h
 * @brief Interface to access Jeveux objects with numpy arrays
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

/* person_in_charge: mathieu.courtois@edf.fr */

#include "astercxx.h"

#include "aster_numpy.h"

#include "MemoryManager/JeveuxVector.h"

#include <string>

/* Definition of numpy type for each JeveuxVector */
template < typename T >
struct npy_type {
  public:
    int value;
};
template <>
struct npy_type< JeveuxVectorLong > {
    static const int value = NPY_LONG;
};
template <>
struct npy_type< JeveuxVectorShort > {
    static const int value = NPY_INT;
};
template <>
struct npy_type< JeveuxVectorReal > {
    static const int value = NPY_DOUBLE;
};
template <>
struct npy_type< JeveuxVectorComplex > {
    static const int value = NPY_CDOUBLE;
};
template <>
struct npy_type< JeveuxVectorLogical > {
    static const int value = NPY_BOOL;
};
template <>
struct npy_type< JeveuxVectorChar8 > {
    static const int value = NPY_STRING;
};
template <>
struct npy_type< JeveuxVectorChar16 > {
    static const int value = NPY_STRING;
};
template <>
struct npy_type< JeveuxVectorChar24 > {
    static const int value = NPY_STRING;
};
template <>
struct npy_type< JeveuxVectorChar32 > {
    static const int value = NPY_STRING;
};
template <>
struct npy_type< JeveuxVectorChar80 > {
    static const int value = NPY_STRING;
};
template <>
struct npy_type< ASTERINTEGER > {
    static const int value = NPY_LONG;
};
template <>
struct npy_type< ASTERDOUBLE > {
    static const int value = NPY_DOUBLE;
};
template <>
struct npy_type< ASTERCOMPLEX > {
    static const int value = NPY_CDOUBLE;
};
template <>
struct npy_type< ASTERINTEGER4 > {
    static const int value = NPY_INT;
};

#endif /* NUMPYACCESS_H_ */
