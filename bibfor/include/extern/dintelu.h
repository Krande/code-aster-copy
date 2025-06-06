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

interface
    subroutine dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                       gammas, gammac, clacier, eys, typdiag, uc, &
                       ntot, dnsinf, dnssup, nrd, mrd, ndemi) bind(C)
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int64_t), intent(in) :: typco
        real(c_double), intent(in) :: alphacc
        real(c_double), intent(in) :: ht
        real(c_double), intent(in) :: bw
        real(c_double), intent(in) :: enrobi
        real(c_double), intent(in) :: enrobs
        real(c_double), intent(in) :: facier
        real(c_double), intent(in) :: fbeton
        real(c_double), intent(in) :: gammas
        real(c_double), intent(in) :: gammac
        integer(c_int64_t), intent(in) :: clacier
        real(c_double), intent(in) :: eys
        integer(c_int64_t), intent(in) :: typdiag
        integer(c_int64_t), intent(in) :: uc
        integer(c_int64_t), intent(inout) :: ntot
        real(c_double), intent(in), optional :: dnsinf
        real(c_double), intent(in), optional :: dnssup
        real(c_double), intent(out), optional :: nrd(1:ntot)
        real(c_double), intent(out), optional :: mrd(1:ntot)
        integer(c_int64_t), intent(out), optional :: ndemi
    end subroutine dintelu
end interface
