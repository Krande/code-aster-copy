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
    subroutine glbels(typco, cequi, effrts, ht, bw,&
                  enrobyi, enrobys, enrobzi, enrobzs,&
                  facier, fbeton, sigcyi, sigcys, sigczi, sigczs, sigs,&
                  precs, flongi, ftrnsv, ferrsyme, slsyme, ferrcomp,&
                  epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                  dnsits, ierr)
    integer :: typco
    real(kind=8) :: cequi
    real(kind=8) :: effrts(6)
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobyi
    real(kind=8) :: enrobys
    real(kind=8) :: enrobzi
    real(kind=8) :: enrobzs    
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: sigcyi
    real(kind=8) :: sigcys
    real(kind=8) :: sigczi
    real(kind=8) :: sigczs
    real(kind=8) :: sigs
    integer :: precs
    integer :: flongi
    integer :: ftrnsv
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: ferrcomp
    integer :: epucisa
    integer :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer :: compress
    integer :: uc
    integer :: um
    real(kind=8) :: dnsits(6)
    integer :: ierr
    end subroutine glbels
end interface
