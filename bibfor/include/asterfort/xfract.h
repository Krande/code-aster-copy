! --------------------------------------------------------------------
! Copyright (C) 1991 - 2019 - EDF R&D - www.code-aster.org
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
    subroutine xfract(ds_thm, nvec, nnop, nnops, nddls, nddlm,&
                      ndim, pla, deplp, deplm,&
                      ffp, ffc, dffc, saut, gradpf,&
                      q1, q2, dpf, q1m, q2m, sautm,&
                      gradpfm, pf, ffp2, psup, pinf,&
                      job, jmate,&
                      t, dimuel, lamb, jheavn, ncompn,&
                      ifiss, nfiss, nfh, ifa, jheafa,&
                      ncomph, contac, depl0, depl1, lambm, pfm)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8), intent(in) :: nvec
        integer(kind=8), intent(in) :: nddlm
        integer(kind=8), intent(in) :: nnop
        integer(kind=8), intent(in) :: nnops
        integer(kind=8), intent(in) :: nddls
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: pla(27)
        integer(kind=8), intent(in) :: dimuel
        real(kind=8), intent(in) :: deplm(dimuel)
        real(kind=8), intent(in) :: deplp(dimuel)
        real(kind=8), intent(in) :: ffp(27)
        real(kind=8), intent(in) :: ffc(16)
        real(kind=8), intent(in) :: dffc(16,3)
        real(kind=8), intent(in) :: ffp2(27)
        character(len=8), intent(in) :: job
        integer(kind=8), intent(in) :: jmate
        real(kind=8), intent(out) :: lamb(3)
        real(kind=8), intent(out) :: t
        real(kind=8), intent(out) :: psup
        real(kind=8), intent(out) :: pinf
        real(kind=8), intent(out) :: saut(3)
        real(kind=8), intent(out) :: gradpf(3)
        real(kind=8), intent(out) :: q1
        real(kind=8), intent(out) :: q2
        real(kind=8), intent(out) :: dpf
        real(kind=8), intent(out) :: q1m
        real(kind=8), intent(out) :: q2m
        real(kind=8), intent(out) :: sautm(3)
        real(kind=8), intent(out) :: gradpfm(3)
        real(kind=8), intent(out) :: pf
        integer(kind=8), intent(in) :: nfiss
        integer(kind=8), intent(in) :: jheavn
        integer(kind=8) :: ifiss
        integer(kind=8), intent(in) :: ncompn
        integer(kind=8) :: nfh
        integer(kind=8) :: ifa
        integer(kind=8) :: jheafa
        integer(kind=8), intent(in) :: ncomph
        integer(kind=8), intent(in) :: contac
        real(kind=8), optional, intent(in) :: depl0(dimuel)
        real(kind=8), optional, intent(in) :: depl1(dimuel)
        real(kind=8), optional, intent(out) :: lambm(3)
        real(kind=8), optional, intent(out) :: pfm
    end subroutine xfract
end interface
