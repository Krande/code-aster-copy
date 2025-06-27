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
    subroutine xmoimp(nh8, nh20, np6, np15, np5,&
                      np13, nt4, nt10, ncpq4, ncpq8,&
                      ncpt3, ncpt6, ndpq4, ndpq8, ndpt3,&
                      ndpt6, nf4, nf8, nf3, nf6,&
                      npf2, npf3, naxt3, naxq4, naxq8,&
                      naxt6, nax2, nax3, nth8, ntp6,&
                      ntp5, ntt4, ntpq4, ntpt3, ntaq4,&
                      ntat3, ntf4, ntf3, ntpf2, ntax2,&
                      nhyq8, nhyt6, nhymq8, nhymt6, nhysq8,&
                      nhyst6, nhydq8, nhydt6, nphm, nhe20,&
                      npe15, npy13, nte10, nhem20, npem15,&
                      npym13, ntem10, nhes20, npes15, npys13,&
                      ntes10, nhed20,nped15, npyd13,&
                      nted10, nbhm, nchm)
        integer(kind=8) :: nh8(15)
        integer(kind=8) :: nh20(7)
        integer(kind=8) :: np6(15)
        integer(kind=8) :: np15(7)
        integer(kind=8) :: np5(15)
        integer(kind=8) :: np13(7)
        integer(kind=8) :: nt4(15)
        integer(kind=8) :: nt10(7)
        integer(kind=8) :: ncpq4(15)
        integer(kind=8) :: ncpq8(7)
        integer(kind=8) :: ncpt3(15)
        integer(kind=8) :: ncpt6(7)
        integer(kind=8) :: ndpq4(15)
        integer(kind=8) :: ndpq8(7)
        integer(kind=8) :: ndpt3(15)
        integer(kind=8) :: ndpt6(7)
        integer(kind=8) :: nf4(11)
        integer(kind=8) :: nf8(7)
        integer(kind=8) :: nf3(11)
        integer(kind=8) :: nf6(7)
        integer(kind=8) :: npf2(11)
        integer(kind=8) :: npf3(7)
        integer(kind=8) :: naxt3(7)
        integer(kind=8) :: naxq4(7)
        integer(kind=8) :: naxq8(7)
        integer(kind=8) :: naxt6(7)
        integer(kind=8) :: nax2(7)
        integer(kind=8) :: nax3(7)
        integer(kind=8) :: nth8(7)
        integer(kind=8) :: ntp6(7)
        integer(kind=8) :: ntp5(7)
        integer(kind=8) :: ntt4(7)
        integer(kind=8) :: ntpq4(7)
        integer(kind=8) :: ntpt3(7)
        integer(kind=8) :: ntaq4(7)
        integer(kind=8) :: ntat3(7)
        integer(kind=8) :: ntf4(7)
        integer(kind=8) :: ntf3(7)
        integer(kind=8) :: ntpf2(7)
        integer(kind=8) :: ntax2(7)
        integer(kind=8) :: nhyq8(17)
        integer(kind=8) :: nhyt6(17)
        integer(kind=8) :: nhymq8(7)
        integer(kind=8) :: nhymt6(7)
        integer(kind=8) :: nhysq8(7)
        integer(kind=8) :: nhyst6(7)
        integer(kind=8) :: nhydq8(7)
        integer(kind=8) :: nhydt6(7)
        integer(kind=8) :: nphm(17)
        integer(kind=8) :: nhe20(17)
        integer(kind=8) :: nhem20(7)
        integer(kind=8) :: nhed20(7)
        integer(kind=8) :: nhes20(7)
        integer(kind=8) :: npe15(17)
        integer(kind=8) :: npem15(7)
        integer(kind=8) :: npes15(7)
        integer(kind=8) :: nped15(7)
        integer(kind=8) :: npy13(17)
        integer(kind=8) :: npym13(7)
        integer(kind=8) :: npys13(7)
        integer(kind=8) :: npyd13(7)
        integer(kind=8) :: nte10(17)
        integer(kind=8) :: ntem10(7)
        integer(kind=8) :: ntes10(7)
        integer(kind=8) :: nted10(7)
        integer(kind=8) :: nbhm(17)
        integer(kind=8) :: nchm(17)
    end subroutine xmoimp
end interface 
