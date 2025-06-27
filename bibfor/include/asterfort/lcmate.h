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
interface
    subroutine lcmate(fami, kpg, ksp, comp, mod, &
                      imat, nmat, tempd, tempf, tref, impexp, &
                      typma, hsr, materd, materf, matcst, &
                      nbcomm, cpmono, angmas, pgl, itmax, &
                      toler, ndt, ndi, nr, crit, &
                      nvi, vind, nfs, nsg, toutms, &
                      nhsr, numhsr, sigd, mult_comp_)
        integer(kind=8) :: nmat
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=16) :: comp(*)
        character(len=8) :: mod
        integer(kind=8) :: imat
        real(kind=8) :: tempd
        real(kind=8) :: tempf
        real(kind=8) :: tref
        integer(kind=8) :: impexp
        character(len=8) :: typma
        real(kind=8) :: hsr(*)
        real(kind=8) :: materd(nmat, 2)
        real(kind=8) :: materf(nmat, 2)
        character(len=3) :: matcst
        integer(kind=8) :: nbcomm(*)
        character(len=24) :: cpmono(*)
        real(kind=8) :: angmas(3)
        real(kind=8) :: pgl(3, 3)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: ndt
        integer(kind=8) :: ndi
        integer(kind=8) :: nr
        real(kind=8) :: crit(*)
        integer(kind=8), intent(in) :: nvi
        real(kind=8) :: vind(*)
        integer(kind=8) :: nfs
        integer(kind=8) :: nsg
        real(kind=8) :: toutms(*)
        integer(kind=8) :: nhsr
        integer(kind=8) :: numhsr(*)
        real(kind=8) :: sigd(6)
        character(len=16), optional, intent(in) :: mult_comp_
    end subroutine lcmate
end interface
