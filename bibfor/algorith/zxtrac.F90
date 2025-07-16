! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
subroutine zxtrac(interp, prec, crit, nbinst, ti, &
                  temps, y, neq, xtract, ier)
    implicit none
#include "asterfort/utmess.h"
#include "blas/zcopy.h"
    integer(kind=8) :: nbinst, neq, ier
    real(kind=8) :: prec, ti(*)
    complex(kind=8) :: xtract(*), y(*)
    character(len=*) :: interp, crit
!     EXTRACTION DANS UN TABLEAU CONTENANT DES VECTEURS A DES INSTANTS
!     SUCESSIFS DU VECTEUR EVENTUELLEMENT INTERPOLLE A L INSTANT SOUHAIT
!-----------------------------------------------------------------------
! IN  : INTERP  : TYPE D'INTERPOLATION
! IN  : PREC    : PRECISION DU TEST
! IN  : CRIT    : CRITERE 'ABSOLU' OU 'RELATIF'
! IN  : NBINST  : DIMENSION DE LA LISTE DES INSTANTS
! IN  : TI      : LISTE DES INSTANTS
! IN  : TEMPS   : TEMPS A INTERPOLER
! IN  : Y       : TABLEAU DE VECTEURS A DES INSTANTS DONNES
! IN  : NEQ     : DIMENSION DES VECTEURS
! OUT : XTRACT  : VECTEUR INTERPOLE AU TEMPS TEMPS
! OUT : IER     : CODE RETOUR, = 0 : IL Y A EU EXTRACTION
!-----------------------------------------------------------------------
    real(kind=8) :: prec2
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: temps
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    ier = 0
!
!     --- RECUPERATION DU CHAMP ---
!
    prec2 = prec
    if (crit(1:7) .eq. 'RELATIF') prec2 = prec*ti(1)
    if (abs(temps-ti(1)) .le. prec2) then
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, y(1), b_incx, xtract, b_incy)
        goto 999
    end if
    if (crit(1:7) .eq. 'RELATIF') prec2 = prec*ti(nbinst)
    if (abs(temps-ti(nbinst)) .le. prec2) then
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, y((nbinst-1)*neq+1), b_incx, xtract, b_incy)
        goto 999
    end if
!
    if (temps .lt. ti(1)) then
        ier = ier+1
        goto 999
    end if
    if (temps .gt. ti(nbinst)) then
        ier = ier+1
        goto 999
    end if
    if (interp(1:3) .eq. 'NON') then
!
!        --- PAS D'INTERPOLATION ---
        do i = 2, nbinst-1
            if (crit(1:7) .eq. 'RELATIF') prec2 = prec*ti(i)
            if (abs(temps-ti(i)) .le. prec2) then
                b_n = to_blas_int(neq)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, y((i-1)*neq+1), b_incx, xtract, b_incy)
                goto 999
            end if
        end do
        ier = ier+1
    else
!        ON INTERDIT L'INTERPOLATION POUR LES COMPLEXES
        call utmess('E', 'ALGORITH11_79')
    end if
!
999 continue
end subroutine
