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
subroutine dimvoi(nvtot, nvoima, nscoma, touvoi, dimvlo)
!
!    MAILLE M0 DE CONNECTIVITE ADCOM0
!    CALCULE LA TAILLE DES DONNEES DU VOISINAGE LOCAL D UN ELEMENT
!    A PARTIR DE TOUVOI
!    IN : DIM,M0,ADCOM0,IATYMA,
!         NVTOT,NVOIMA,NSCOMA,TOUVOI
!    OUT : DIMVLO
!
    implicit none
    integer(kind=8) :: nvtot, nvoima, nscoma
    integer(kind=8) :: touvoi(1:nvoima, 1:nscoma+2)
    integer(kind=8) :: dimvlo
    integer(kind=8) :: iv, nsco
!
    dimvlo = 1+5*nvtot
    if (nvtot .ge. 1) then
        do iv = 1, nvtot
            nsco = touvoi(iv, 2)
            dimvlo = dimvlo+2*nsco
        end do
    end if
!
end subroutine
