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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmiclb(fami, kpg, ksp, option, rela_comp,&
                      imate, xlong0, aire, tmoins, tplus,&
                      dlong0, effnom, vim, effnop, vip,&
                      klv, fono, epsm, carcri, codret)
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=16) :: option
        character(len=16) :: rela_comp
        integer(kind=8) :: imate
        real(kind=8) :: xlong0
        real(kind=8) :: aire
        real(kind=8) :: tmoins
        real(kind=8) :: tplus
        real(kind=8) :: dlong0
        real(kind=8) :: effnom
        real(kind=8) :: vim(*)
        real(kind=8) :: effnop
        real(kind=8) :: vip(*)
        real(kind=8) :: klv(21)
        real(kind=8) :: fono(6)
        real(kind=8) :: epsm
        real(kind=8) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: codret
    end subroutine nmiclb
end interface
