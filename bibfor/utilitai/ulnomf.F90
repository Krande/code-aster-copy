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
function ulnomf(nomfic, kacc, typef)
    implicit none
    integer(kind=8) :: ulnomf
    character(len=*) :: nomfic, kacc, typef
!     ------------------------------------------------------------------
!     RETOURNE LE NUMERO D'UNITE LOGIQUE ASSOCIE AU NOM SYSTEME
!              -1 SI AUCUN DE DISPONIBLE
!     RENVOIE LE TYPE D'ACCES AU FICHIER ASSOCIE DANS L'ARGUMENT KACC
!     RENVOIE LE TYPE DE FICHIER ASSOCIE DANS L'ARGUMENT TYPEF
! person_in_charge: j-pierre.lefebvre at edf.fr
    integer(kind=8) :: mxf
    parameter(mxf=100)
    character(len=1) :: typefi(mxf), accefi(mxf), etatfi(mxf), modifi(mxf)
    character(len=16) :: ddname(mxf)
    character(len=255) :: namefi(mxf)
    integer(kind=8) :: first, unitfi(mxf), nbfile
    common/asgfi1/first, unitfi, nbfile
    common/asgfi2/namefi, ddname, typefi, accefi, etatfi, modifi
!
    integer(kind=8) :: ival, k
!
    ival = -1
    kacc = '?'
    typef = '?'
    do k = 1, mxf-1
        if (typefi(k) .ne. '?') then
            if (namefi(k) .eq. nomfic) then
                ival = unitfi(k)
                typef = typefi(k)
                kacc = accefi(k)
                goto 2
            end if
        end if
    end do
2   continue
    ulnomf = ival
end function
