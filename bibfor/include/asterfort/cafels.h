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
!
interface
    subroutine cafels(cequi, effm, effn, ht, bw, &
                      enrobi, enrobs, scmaxi, scmaxs, ssmax, &
                      ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                      dnsinf, dnssup, sigmsi, sigmss, &
                      sigmci, sigmcs, &
                      alpha, pivot, etat, ierr)
        real(kind=8) :: cequi
        real(kind=8) :: effm
        real(kind=8) :: effn
        real(kind=8) :: ht
        real(kind=8) :: bw
        real(kind=8) :: enrobi
        real(kind=8) :: enrobs
        real(kind=8) :: scmaxi
        real(kind=8) :: scmaxs
        real(kind=8) :: ssmax
        integer(kind=8) :: ferrcomp
        integer(kind=8) :: precs
        integer(kind=8) :: ferrsyme
        real(kind=8) :: slsyme
        integer(kind=8) :: uc
        integer(kind=8) :: um
        real(kind=8) :: dnsinf
        real(kind=8) :: dnssup
        real(kind=8) :: sigmsi
        real(kind=8) :: sigmss
        real(kind=8) :: sigmci
        real(kind=8) :: sigmcs
        real(kind=8) :: alpha
        integer(kind=8) :: pivot
        integer(kind=8) :: etat
        integer(kind=8) :: ierr
    end subroutine cafels
end interface
