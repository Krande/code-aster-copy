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
    subroutine clcplq(typcmb, typco, nb, precs, &
                      flongi, ftrnsv, ferrsyme, slsyme, ferrcomp, epucisa,&
                      ferrmin, rholmin, rhotmin, compress, cequi,&
                      enrobi, enrobs, sigs, sigci, sigcs,&
                      alphacc, gammas, gammac, facier, eys, typdiag,&
                      fbeton, clacier, uc, um,&
                      wmaxi, wmaxs, sigelsqp, kt, phixi, phixs, phiyi, phiys,&
                      ht, effrts, dnsits, ierr)
        integer :: typcmb
        integer :: typco
        integer :: nb
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
        real(kind=8) :: cequi
        real(kind=8) :: enrobi
        real(kind=8) :: enrobs
        real(kind=8) :: sigs
        real(kind=8) :: sigci
        real(kind=8) :: sigcs    
        real(kind=8) :: alphacc
        real(kind=8) :: gammas
        real(kind=8) :: gammac
        real(kind=8) :: facier
        real(kind=8) :: eys
        integer :: typdiag
        real(kind=8) :: fbeton
        integer :: clacier
        integer :: uc
        integer :: um
        real(kind=8) :: wmaxi
        real(kind=8) :: wmaxs
        real(kind=8) :: sigelsqp
        real(kind=8) :: kt
        real(kind=8) :: phixi
        real(kind=8) :: phixs
        real(kind=8) :: phiyi
        real(kind=8) :: phiys
        real(kind=8) :: ht
        real(kind=8) :: effrts(8)
        real(kind=8) :: dnsits(6)
        integer :: ierr
    end subroutine clcplq
end interface
