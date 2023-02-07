/**
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

#include "Utilities/Blas.h"

#include "aster_fort_blas.h"

void AsterBLAS::scal( const ASTERINTEGER &n, const ASTERDOUBLE &alpha, ASTERDOUBLE *x,
                      const ASTERINTEGER &incx ) {
    CALL_DSCAL( &n, &alpha, x, &incx );
};

void AsterBLAS::scal( const ASTERINTEGER &n, const ASTERCOMPLEX &alpha, ASTERCOMPLEX *x,
                      const ASTERINTEGER &incx ) {
    CALL_ZSCAL( &n, &alpha, x, &incx );
};

void AsterBLAS::scal( const ASTERINTEGER &n, const ASTERDOUBLE &alpha, ASTERCOMPLEX *x,
                      const ASTERINTEGER &incx ) {
    CALL_ZDSCAL( &n, &alpha, x, &incx );
};
