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

subroutine verifepsa(fami, kpg, ksp, poum, &
                     epsa, iepsa_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/rcvarc.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    character(len=*), intent(in) :: poum
    real(kind=8), intent(out) :: epsa(6)
    integer(kind=8), optional, intent(out) :: iepsa_(6)
!
! --------------------------------------------------------------------------------------------------
!
! Get anelastic deformation (defined as external state variable)
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  poum         : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature => epsa is increment
! In  ndim         : space dimension
! Out epsa         : anelastic strain
! Out iepsa_       : 0 if defined
!                    1 if not
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, iretm(6), iretp(6), iret(6)
    real(kind=8) :: defam(6), defap(6)
    character(len=6), parameter :: name_epsa(6) = (/'EPSAXX', 'EPSAYY', 'EPSAZZ', &
                                                    'EPSAXY', 'EPSAXZ', 'EPSAYZ'/)
!
! --------------------------------------------------------------------------------------------------
!
    defam = 0.d0
    defap = 0.d0
    epsa = 0.d0
!
! - Get anelastic strains
!
    if (poum .eq. 'T' .or. poum .eq. '-') then
        do i = 1, 6
            call rcvarc(' ', name_epsa(i), '-', fami, kpg, &
                        ksp, defam(i), iretm(i))
            if (iretm(i) .ne. 0) then
                defam(i) = 0.d0
            end if
        end do
    end if
!
    if (poum .eq. 'T' .or. poum .eq. '+') then
        do i = 1, 6
            call rcvarc(' ', name_epsa(i), '+', fami, kpg, &
                        ksp, defap(i), iretp(i))
            if (iretp(i) .ne. 0) then
                defap(i) = 0.d0
            end if
        end do
    end if
!
! - Compute poum strains
!
    if (poum .eq. 'T') then
        do i = 1, 6
            iret(i) = iretm(i)+iretp(i)
            if (iret(i) .eq. 0) then
                epsa(i) = defap(i)-defam(i)
            end if
        end do
!
    else if (poum .eq. '-') then
        do i = 1, 6
            iret(i) = iretm(i)
            if (iret(i) .eq. 0) then
                epsa(i) = defam(i)
            end if
        end do
!
    else if (poum .eq. '+') then
        do i = 1, 6
            iret(i) = iretp(i)
            if (iret(i) .eq. 0) then
                epsa(i) = defap(i)
            end if
        end do
!
    else
        ASSERT(((poum .eq. 'T') .or. (poum .eq. '-') .or. (poum .eq. '+')))
    end if
!
! - Output errors
!
    if (present(iepsa_)) then
        do i = 1, 6
            iepsa_(i) = iret(i)
        end do
    end if
!
end subroutine
