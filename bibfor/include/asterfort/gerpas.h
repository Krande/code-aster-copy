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
    subroutine gerpas(fami, kpg, ksp, rela_comp, mod,&
                      imat, matcst, nbcomm, cpmono, nbphas,&
                      nvi, nmat, y, pas, itmax,&
                      eps, toly, cothe, coeff, dcothe,&
                      dcoeff, coel, pgl, angmas, neps,&
                      epsd, detot, x, nfs, nsg,&
                      nhsr, numhsr, hsr, iret)
        integer(kind=8) :: nhsr
        integer(kind=8) :: nsg
        integer(kind=8) :: nfs
        integer(kind=8) :: neps
        integer(kind=8) :: nmat
        integer(kind=8) :: nvi
        integer(kind=8) :: nbphas
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=16) :: rela_comp
        character(len=8) :: mod
        integer(kind=8) :: imat
        character(len=3) :: matcst
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=24) :: cpmono(5*nmat+1)
        real(kind=8) :: y(nvi)
        real(kind=8) :: pas
        integer(kind=8) :: itmax
        real(kind=8) :: eps
        real(kind=8) :: toly
        real(kind=8) :: cothe(nmat)
        real(kind=8) :: coeff(nmat)
        real(kind=8) :: dcothe(nmat)
        real(kind=8) :: dcoeff(nmat)
        real(kind=8) :: coel(nmat)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: angmas(3)
        real(kind=8) :: epsd(neps)
        real(kind=8) :: detot(neps)
        real(kind=8) :: x
        integer(kind=8) :: numhsr(*)
        real(kind=8) :: hsr(nsg, nsg, nhsr)
        integer(kind=8) :: iret
    end subroutine gerpas
end interface
