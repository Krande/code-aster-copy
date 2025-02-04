! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
       subroutine cafeluiter(typco, alphacc, effm, effn, ht, bw,&
                      enrobi, enrobs, facier, fbeton, gammas, gammac,&
                      clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                      condns, astend, ascomp, sstend, sscomp, ectend, eccomp,&
                      alpha, pivot, etat, ierr)
    integer :: typco
    real(kind=8) :: alphacc
    real(kind=8) :: effm
    real(kind=8) :: effn
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer :: clacier
    real(kind=8) :: eys
    integer :: typdiag
    integer :: ferrcomp
    integer :: precs
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: uc
    integer :: um
    logical :: condns
    real(kind=8) :: astend
    real(kind=8) :: ascomp
    real(kind=8) :: sstend
    real(kind=8) :: sscomp
    real(kind=8) :: ectend
    real(kind=8) :: eccomp
    real(kind=8) :: alpha
    integer :: pivot
    integer :: etat
    integer :: ierr
    end subroutine cafeluiter
end interface
