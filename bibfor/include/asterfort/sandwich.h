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
     subroutine sandwich(enrobi, enrobs, facier, fbeton, gammas, gammac,&
                    thiter, epiter, aphiter, cond109,&
                    flongi, ftrnsv, ferrcomp, ferrsyme, slsyme,&
                    epucisa, ferrmin, rholmin, rhotmin, compress, uc, um,&
                    ht, effrts, dnsits, ierr)
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    real(kind=8) :: thiter
    real(kind=8) :: epiter
    real(kind=8) :: aphiter
    integer :: cond109
    integer :: flongi
    integer :: ftrnsv
    integer :: ferrcomp
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: epucisa
    integer :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer :: compress
    integer :: uc
    integer :: um
    real(kind=8) :: ht
    real(kind=8) :: effrts(8)
    real(kind=8) :: dnsits(6)
    integer :: ierr
    end subroutine sandwich
end interface
