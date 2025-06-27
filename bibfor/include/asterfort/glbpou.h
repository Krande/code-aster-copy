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
    subroutine glbpou(typcmb, typco, cequi, effrts, ht, bw,&
                  enrobyi, enrobys, enrobzi, enrobzs,&
                  facier, fbeton, sigelsqp, kt, eys,&
                  alphacc, clacier, gammas, gammac, typdiag,&
                  sigcyi, sigcys, sigczi, sigczs, sigs,&
                  wmaxyi, wmaxys, wmaxzi, wmaxzs,&
                  phiyi, phiys, phizi, phizs,&
                  precs, ferrsyme, slsyme, ferrcomp,&
                  epucisa, ferrmin, rholmin, rhotmin, compress, uc, um, &
                  rhoacier, areinf, ashear, astirr, rhocrit, datcrit, lcrit,&
                  dnsits, dnsvol, construc, ierrl, ierrt)
    integer(kind=8) :: typcmb
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
    real(kind=8) :: sigelsqp
    real(kind=8) :: kt
    real(kind=8) :: eys
    real(kind=8) :: alphacc
    integer(kind=8) :: clacier
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer(kind=8) :: typdiag
    real(kind=8) :: sigcyi
    real(kind=8) :: sigcys
    real(kind=8) :: sigczi
    real(kind=8) :: sigczs
    real(kind=8) :: sigs
    real(kind=8) :: wmaxyi
    real(kind=8) :: wmaxys
    real(kind=8) :: wmaxzi
    real(kind=8) :: wmaxzs
    real(kind=8) :: phiyi
    real(kind=8) :: phiys
    real(kind=8) :: phizi
    real(kind=8) :: phizs
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
    real(kind=8) :: rhoacier
    real(kind=8) :: areinf
    real(kind=8) :: ashear
    real(kind=8) :: astirr
    real(kind=8) :: rhocrit
    real(kind=8) :: datcrit
    real(kind=8) :: lcrit
    real(kind=8) :: dnsits(6)
    real(kind=8) :: dnsvol
    real(kind=8) :: construc
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt
    end subroutine glbpou
end interface
