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
function ulnume()
    implicit none
    integer(kind=8) :: ulnume
!     ------------------------------------------------------------------
!     RETOURNE UN NUMERO D'UNITE LOGIQUE NON UTILISE
!              -1 SI AUCUN DE DISPONIBLE
!     LA RECHERCHE EST LIMITEE A L'INTERVALLE 70-99
! person_in_charge: j-pierre.lefebvre at edf.fr
!
#include "asterfort/ulinit.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: mxf
    parameter(mxf=100)
    character(len=1) :: typefi(mxf), accefi(mxf), etatfi(mxf), modifi(mxf)
    character(len=16) :: ddname(mxf)
    character(len=255) :: namefi(mxf)
    integer(kind=8) :: first, unitfi(mxf), nbfile
    common/asgfi1/first, unitfi, nbfile
    common/asgfi2/namefi, ddname, typefi, accefi, etatfi, modifi
!
    integer(kind=8) :: i, ival, k
!
    if (first .ne. 17111990) call ulinit()
!
    ival = -1
    do i = 99, 70, -1
        do k = 1, nbfile
            if (unitfi(k) .eq. i) then
                goto 1
            end if
        end do
        ival = i
        goto 2
1       continue
    end do
    call utmess('A', 'UTILITAI5_10')
2   continue
    ulnume = ival
end function
