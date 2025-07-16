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
    subroutine glbels(typco, cequi, effrts, ht, bw,&
                  enrobyi, enrobys, enrobzi, enrobzs,&
                  facier, fbeton, sigcyi, sigcys, sigczi, sigczs, sigs,&
                  precs, ferrsyme, slsyme, ferrcomp,&
                  epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                  dnsits, ierrl, ierrt)
    integer(kind=8) :: typco
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
    integer(kind=8) :: precs
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: epucisa
    integer(kind=8) :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer(kind=8) :: compress
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: dnsits(6)
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt
    end subroutine glbels
end interface
