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
    subroutine sandcas4(effrts, ht, enrobi, enrobs, facier, fbeton, gammas, gammac,&
                    thiter, epiter, cond109, ferrcomp, ferrsyme, slsyme, uc, um,&
                    dnsxi, dnsxs, dnsyi, dnsys, etsxi, etsxs, etsyi, etsys,&
                    snsxi, snsxs, snsyi, snsys, ncmaxi, ncmini, ncmaxs, ncmins,&
                    t_inf, t_sup, theta_inf, theta_sup, ierr)
    real(kind=8) :: effrts(6)
    real(kind=8) :: ht
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    real(kind=8) :: thiter
    real(kind=8) :: epiter
    integer(kind=8) :: cond109
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: dnsxi
    real(kind=8) :: dnsxs
    real(kind=8) :: dnsyi
    real(kind=8) :: dnsys
    integer(kind=8) :: etsxi
    integer(kind=8) :: etsxs
    integer(kind=8) :: etsyi
    integer(kind=8) :: etsys
    real(kind=8) :: snsxi
    real(kind=8) :: snsxs
    real(kind=8) :: snsyi
    real(kind=8) :: snsys
    real(kind=8) :: ncmaxi
    real(kind=8) :: ncmini
    real(kind=8) :: ncmaxs
    real(kind=8) :: ncmins
    real(kind=8) :: t_inf
    real(kind=8) :: t_sup
    real(kind=8) :: theta_inf
    real(kind=8) :: theta_sup
    integer(kind=8) :: ierr
    end subroutine sandcas4
end interface
