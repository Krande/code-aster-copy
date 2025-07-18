! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
! ===============================================================
! THIS LAPACK 2.0 ROUTINE IS DEPRECATED
! DO NOT USE IT : YOU SHOULD PREFER UP-TO-DATE LAPACK ROUTINE
!
! BUT DO NOT REMOVE IT :
! THE PRESENT ROUTINE IS MANDATORY FOR ARPACK LIBRARY
! WHICH STICKS TO LAPACK 2.0 VERSION
! ==============================================================
subroutine ar_dlarfg(n, alpha, x, incx, tau)
!
!     SUBROUTINE LAPACK CALCULANT UN REFLECTEUR H TEL QUE DECRIT
!     CI DESSOUS.
!-----------------------------------------------------------------------
!  -- LAPACK AUXILIARY ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     SEPTEMBER 30, 1994
!
!  PURPOSE
!  =======
!
!  DLARFG GENERATES A REAL ELEMENTARY REFLECTOR H OF ORDER N, SUCH
!  THAT
!
!        H * ( ALPHA ) = ( BETA ),   H' * H = I.
!            (   X   )   (   0  )
!
!  WHERE ALPHA AND BETA ARE SCALARS, AND X IS AN (N-1)-ELEMENT REAL
!  VECTOR. H IS REPRESENTED IN THE FORM
!
!        H = I - TAU * ( 1 ) * ( 1 V' ) ,
!                      ( V )
!
!  WHERE TAU IS A REAL SCALAR AND V IS A REAL (N-1)-ELEMENT
!  VECTOR.
!
!  IF THE ELEMENTS OF X ARE ALL ZERO, THEN TAU = 0 AND H IS TAKEN TO BE
!  THE UNIT MATRIX.
!
!  OTHERWISE  1 <= TAU <= 2.
!
!  ARGUMENTS
!  =========
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE ELEMENTARY REFLECTOR.
!
!  ALPHA   (INPUT/OUTPUT) REAL*8
!          ON ENTRY, THE VALUE ALPHA.
!          ON EXIT, IT IS OVERWRITTEN WITH THE VALUE BETA.
!
!  X       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION
!                         (1+(N-2)*ABS(INCX))
!          ON ENTRY, THE VECTOR X.
!          ON EXIT, IT IS OVERWRITTEN WITH THE VECTOR V.
!
!  INCX    (INPUT) INTEGER
!          THE INCREMENT BETWEEN ELEMENTS OF X. INCX > 0.
!
!  TAU     (OUTPUT) REAL*8
!          THE VALUE TAU.
!
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE DLAMCH PAR R8PREM ET R8MIEM,
!            REMPLACEMENT DE RETURN PAR GOTO 1000,
!            MODIFICATION DES APPELS BLAS (ROUTINE ASTER BL...),
!            IMPLICIT NONE.
! INTRINSIC FUNCTION
!    ABS, SIGN
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "blas/dlapy2.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
    integer(kind=8) :: incx, n
    real(kind=8) :: alpha, tau
!     ..
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: x(*)
!
!     .. PARAMETERS ..
    real(kind=8) :: one, zero
    parameter(one=1.0d+0, zero=0.0d+0)
!     ..
!     .. LOCAL SCALARS ..
    integer(kind=8) :: j, knt
    real(kind=8) :: beta, rsafmn, safmin, xnorm
    blas_int :: b_incx, b_n
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
    if (n .le. 1) then
        tau = zero
        goto 1000
    end if
!
    b_n = to_blas_int(n-1)
    b_incx = to_blas_int(incx)
    xnorm = dnrm2(b_n, x, b_incx)
!
    if (xnorm .eq. zero) then
!
!        H  =  I
!
        tau = zero
    else
!
!        GENERAL CASE
!
        beta = -sign(dlapy2(alpha, xnorm), alpha)
        safmin = r8miem()/(r8prem()*0.5d0)
        if (abs(beta) .lt. safmin) then
!
!           XNORM, BETA MAY BE INACCURATE, SCALE X AND RECOMPUTE THEM
!
            rsafmn = one/safmin
            knt = 0
10          continue
            knt = knt+1
            b_n = to_blas_int(n-1)
            b_incx = to_blas_int(incx)
            call dscal(b_n, rsafmn, x, b_incx)
            beta = beta*rsafmn
            alpha = alpha*rsafmn
            if (abs(beta) .lt. safmin) goto 10
!
!           NEW BETA IS AT MOST 1, AT LEAST SAFMIN
!
            b_n = to_blas_int(n-1)
            b_incx = to_blas_int(incx)
            xnorm = dnrm2(b_n, x, b_incx)
            beta = -sign(dlapy2(alpha, xnorm), alpha)
            tau = (beta-alpha)/beta
            b_n = to_blas_int(n-1)
            b_incx = to_blas_int(incx)
            call dscal(b_n, one/(alpha-beta), x, b_incx)
!
!           IF ALPHA IS SUBNORMAL, IT MAY LOSE RELATIVE ACCURACY
!
            alpha = beta
            do j = 1, knt
                alpha = alpha*safmin
            end do
        else
            tau = (beta-alpha)/beta
            b_n = to_blas_int(n-1)
            b_incx = to_blas_int(incx)
            call dscal(b_n, one/(alpha-beta), x, b_incx)
            alpha = beta
        end if
    end if
!
1000 continue
    call matfpe(1)
!
!     END OF DLARFG
!
end subroutine
