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
function meiden(scal, ncmp, i1, i3, nec, &
                i2, i4)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
    aster_logical :: meiden
    character(len=4) :: scal
    integer(kind=8) :: ncmp, i1, i3, nec, i2, i4
! ----------------------------------------------------------------------
!     ENTREES:
!        SCAL : R, I , C, K8, K16, K24
!        NCMP : NOMBRE DE COMPOSANTES DES GRANDEURS
!          I1 : ADRESSE DANS ZR OU ZI ... DU DEBUT DE LA 1ERE GRANDEUR
!          I3 : ADRESSE DANS ZR OU ZI ... DU DEBUT DE LA 2EME GRANDEUR
!        NEC  : NOMBRE D'ENTIERS CODES
!          I2 : ADRESSE DANS ZI DU DEBUT DU DG DE LA 1ERE GRANDEUR
!          I4 : ADRESSE DANS ZI DU DEBUT DU DG DE LA 2EME GRANDEUR
!
!     SORTIES:
!     MEIDEN : VRAI SI LES 2 GRANDEURS SONT IDENTIQUES.
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
!
! DEB-------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iec
!-----------------------------------------------------------------------
    meiden = .false.
!
!     -- ON TESTE D'ABORD L'EGALITE DES DESCIPTEUR GRANDEUR:
    do iec = 1, nec
        if (zi(i2+iec) .ne. zi(i4+iec)) goto 999
    end do
!
!     -- ON TESTE ENSUITE LES VALEURS:
    if (scal(1:1) .eq. 'I') then
        do i = 1, ncmp
            if (zi(i1+i) .ne. zi(i3+i)) goto 999
        end do
    else if (scal(1:1) .eq. 'R') then
        do i = 1, ncmp
            if (zr(i1+i) .ne. zr(i3+i)) goto 999
        end do
    else if (scal(1:1) .eq. 'C') then
        do i = 1, ncmp
            if (zc(i1+i) .ne. zc(i3+i)) goto 999
        end do
    else if (scal(1:3) .eq. 'K8 ') then
        do i = 1, ncmp
            if (zk8(i1+i) .ne. zk8(i3+i)) goto 999
        end do
    else if (scal(1:3) .eq. 'K16') then
        do i = 1, ncmp
            if (zk16(i1+i) .ne. zk16(i3+i)) goto 999
        end do
    else if (scal(1:3) .eq. 'K24') then
        do i = 1, ncmp
            if (zk24(i1+i) .ne. zk24(i3+i)) goto 999
        end do
    else
        ASSERT(.false.)
    end if
    meiden = .true.
999 continue
end function
