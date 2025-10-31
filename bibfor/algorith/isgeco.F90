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
subroutine isgeco(icod1, icod2, ndim, iopt, icod)
!    P. RICHARD     DATE 19/02/91
!-----------------------------------------------------------------------
!  BUT:  GERER L'ADDITION OU LA SOUSTRACTION DES DEUX ENTIER CODES SUR
!   LES 7 PREMIERES PUISSANCE ENTIERES DE DEUX
    implicit none
!
!   SI IOPT=1     ADDITION ICOD1+ICOD2
!   SI IOPT=-1     ICOD1-ICOD2
!-----------------------------------------------------------------------
!
! ICOD1    /I/: PREMIER ENTIER CODE
! ICOD2    /I/: DEUXIEME  ENTIER CODE
! NDIM     /I/: NOMBRE DE PUISSANCE A DECODER
! IOPT     /I/: OPTION DE CALCUL
! ICOD     /O/: RESULTAT
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
#include "asterfort/iscode.h"
#include "asterfort/isdeco.h"
    integer(kind=8) :: i, ik, iopt, nbcpmx, ndim
!-----------------------------------------------------------------------
    parameter(nbcpmx=320)
    integer(kind=8) :: icod1(1), icod2(1), icod(1)
    integer(kind=8) :: idec1(nbcpmx), idec2(nbcpmx), idec(nbcpmx)
!
!-----------------------------------------------------------------------
!
    call isdeco(icod1, idec1, ndim)
    call isdeco(icod2, idec2, ndim)
!
    if (iopt .eq. 1) then
        do i = 1, ndim
            ik = idec1(i)+idec2(i)
            if (ik .gt. 0) then
                idec(i) = 1
            else
                idec(i) = 0
            end if
        end do
    end if
!
    if (iopt .eq. -1) then
        do i = 1, ndim
            ik = idec1(i)-idec2(i)
            if (ik .gt. 0) then
                idec(i) = 1
            else
                idec(i) = 0
            end if
        end do
    end if
!
    call iscode(idec, icod, ndim)
!
end subroutine
