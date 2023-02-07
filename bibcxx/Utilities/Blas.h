#ifndef BLAS_H_
#define BLAS_H_

/**
 * @file Blas.h
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

#include "astercxx.h"

class AsterBLAS {
  public:
    // generic scal
    template < typename T >
    static void scal( const ASTERINTEGER &n, const T &alpha, T *x, const ASTERINTEGER &incx );

    // dscal
    static void scal( const ASTERINTEGER &n, const ASTERDOUBLE &alpha, ASTERDOUBLE *x,
                      const ASTERINTEGER &incx );

    // zscal
    static void scal( const ASTERINTEGER &n, const ASTERCOMPLEX &alpha, ASTERCOMPLEX *x,
                      const ASTERINTEGER &incx );

    // zdscal
    static void scal( const ASTERINTEGER &n, const ASTERDOUBLE &alpha, ASTERCOMPLEX *x,
                      const ASTERINTEGER &incx );
};

template < typename T >
void AsterBLAS::scal( const ASTERINTEGER &n, const T &alpha, T *x, const ASTERINTEGER &incx ) {
    for ( ASTERINTEGER i = 0; i < n; i += incx )
        x[i] *= alpha;
};

#endif /* BLAS_H_ */
