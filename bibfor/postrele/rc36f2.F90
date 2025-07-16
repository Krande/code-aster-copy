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
subroutine rc36f2(nbsigr, nocc, saltij, i1a4, isk, &
                  isl, nk, nl, n0)
    implicit none
#include "asterf_types.h"
    integer(kind=8) :: nbsigr, nocc(*), i1a4, isk, isl, nk, nl, n0
    real(kind=8) :: saltij(*)
!
!     MISE A ZERO DES LIGNES ET COLONNES DANS SALT POUR LES
!     SITUATION K ET L SI NOCC = 0
!
!     ------------------------------------------------------------------
    integer(kind=8) :: k, l, i1
    aster_logical :: colona, colonb, lignea, ligneb
!     ------------------------------------------------------------------
!
    colona = .false.
    colonb = .false.
    lignea = .false.
    ligneb = .false.
!
    if (i1a4 .eq. 1 .or. i1a4 .eq. 2) then
        nocc(2*(isl-1)+1) = nl-n0
        if (nocc(2*(isl-1)+1) .eq. 0) colona = .true.
    else if (i1a4 .eq. 3 .or. i1a4 .eq. 4) then
        nocc(2*(isl-1)+2) = nl-n0
        if (nocc(2*(isl-1)+2) .eq. 0) colonb = .true.
    end if
!
    if (i1a4 .eq. 1 .or. i1a4 .eq. 3) then
        nocc(2*(isk-1)+1) = nk-n0
        if (nocc(2*(isk-1)+1) .eq. 0) lignea = .true.
    else if (i1a4 .eq. 2 .or. i1a4 .eq. 4) then
        nocc(2*(isk-1)+2) = nk-n0
        if (nocc(2*(isk-1)+2) .eq. 0) ligneb = .true.
    end if
!
    if (colona) then
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            saltij(i1+4*(isl-1)+1) = 0.d0
            saltij(i1+4*(isl-1)+2) = 0.d0
        end do
        i1 = 4*nbsigr*(isl-1)
        do l = 1, nbsigr
            saltij(i1+4*(l-1)+1) = 0.d0
            saltij(i1+4*(l-1)+3) = 0.d0
        end do
    end if
!
    if (colonb) then
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            saltij(i1+4*(isl-1)+3) = 0.d0
            saltij(i1+4*(isl-1)+4) = 0.d0
        end do
        i1 = 4*nbsigr*(isl-1)
        do l = 1, nbsigr
            saltij(i1+4*(l-1)+2) = 0.d0
            saltij(i1+4*(l-1)+4) = 0.d0
        end do
    end if
!
    if (lignea) then
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            saltij(i1+4*(isk-1)+1) = 0.d0
            saltij(i1+4*(isk-1)+2) = 0.d0
        end do
        i1 = 4*nbsigr*(isk-1)
        do l = 1, nbsigr
            saltij(i1+4*(l-1)+1) = 0.d0
            saltij(i1+4*(l-1)+3) = 0.d0
        end do
    end if
!
    if (ligneb) then
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
            saltij(i1+4*(isk-1)+3) = 0.d0
            saltij(i1+4*(isk-1)+4) = 0.d0
        end do
        i1 = 4*nbsigr*(isk-1)
        do l = 1, nbsigr
            saltij(i1+4*(l-1)+2) = 0.d0
            saltij(i1+4*(l-1)+4) = 0.d0
        end do
    end if
!
end subroutine
