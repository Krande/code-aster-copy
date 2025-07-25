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
subroutine isdeco(icod, idec, ndim)
!    P. RICHARD     DATE 06/11/90
!-----------------------------------------------------------------------
!  BUT: DECODER UN ENTIER CODE SUR LES 30 PREMIERES PUISSANCES
!          DE DEUX ( PAS DE PUISSANCE 0)
    implicit none
!-----------------------------------------------------------------------
!
!  ICOD(*)  /I/: ENTIER CODE :
!                ICOD(1) : 30 1ERES CMPS CODE SUR LES PUISS DE 2:1 A 30
!                ICOD(2) : 30 CMPS SUIV CODE SUR LES PUISS DE 2:1 A 30
!                ...
!  IDEC     /O/: VECTEUR DES NDIM PREMIERES CMPS
!  NDIM     /I/: NOMBRE DE CMPS A DECODER
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, necmax
    integer(kind=8) :: idec(ndim), icod(*)
    integer(kind=8) :: nec, iec, i, ipui, k
    parameter(necmax=10)
    integer(kind=8) :: ifin(necmax)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    nec = (ndim-1)/30+1
!
! --- IFIN DONNE POUR CHAQUE ENTIER CODE LE NOMBRE MAX DE CMPS
! --- QUE L'ON PEUT TROUVER SUR CET ENTIER :
!     ------------------------------------
    do iec = 1, nec-1
        ifin(iec) = 30
    end do
    ifin(nec) = ndim-30*(nec-1)
!
    k = 0
    do iec = 1, nec
        ipui = 1
        do i = 1, ifin(iec)
            k = k+1
            ipui = ipui*2
            if (iand(icod(iec), ipui) .eq. ipui) then
                idec(k) = 1
            else
                idec(k) = 0
            end if
        end do
    end do
!
end subroutine
