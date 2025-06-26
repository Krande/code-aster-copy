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
subroutine ulimpr(impr)
    implicit none
    integer(kind=8) :: impr
! person_in_charge: j-pierre.lefebvre at edf.fr
!
!     IMPRESSION DES TABLES DECRIVANT LES UNITES LOGIQUE OUVERTES
!
    integer(kind=8) :: mxf
    parameter(mxf=100)
    character(len=1) :: typefi(mxf), accefi(mxf), etatfi(mxf), modifi(mxf)
    character(len=16) :: ddname(mxf)
    character(len=255) :: namefi(mxf)
    integer(kind=8) :: first, unitfi(mxf), nbfile
    common/asgfi1/first, unitfi, nbfile
    common/asgfi2/namefi, ddname, typefi, accefi, etatfi, modifi
!
    integer(kind=8) :: i
    character(len=8) :: ktyp, kacc, keta
!
    write (impr, 999) 'LA TABLE A CONTENU JUSQU''A ', nbfile,&
     &                ' ASSOCIATION(S)'
    do i = 1, nbfile
        write (impr, 1000) namefi(i)
        ktyp = '?'
        if (typefi(i) .eq. 'A') then
            ktyp = 'ASCII'
        else if (typefi(i) .eq. 'B') then
            ktyp = 'BINARY'
        else if (typefi(i) .eq. 'L') then
            ktyp = 'LIBRE'
        end if
        kacc = '?'
        if (accefi(i) .eq. 'N') then
            kacc = 'NEW'
        else if (accefi(i) .eq. 'O') then
            kacc = 'OLD'
        else if (accefi(i) .eq. 'A') then
            kacc = 'APPEND'
        end if
        keta = '?'
        if (etatfi(i) .eq. 'O') then
            keta = 'OPEN'
        else if (etatfi(i) .eq. 'F') then
            keta = 'CLOSE'
        else if (etatfi(i) .eq. 'R') then
            keta = 'RESERVE '
        end if
        write (impr, 1001) ddname(i), unitfi(i), ktyp, kacc, keta, modifi(i)
    end do
!
999 format(a, i4, a)
1000 format(1x, a)
1001 format(6x, a16, i3, 3(1x, a8), 1x, a1)
end subroutine
