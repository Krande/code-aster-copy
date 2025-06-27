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
    subroutine xvechm(ds_thm, nnops, ddls, ddlm, ndim, pla,&
                      saut, sautm, nd, ffc, w11, w11m, jac,&
                      q1, dt, ta, q1m, ta1, q2, q2m, dffc,&
                      rho11, gradpf, rho11m, gradpfm, ffp2,&
                      vect, ffp,&
                      nnop, delta, lamb, am, r, p, psup,&
                      pinf, pf, ncompn, jheavn, ifiss, nfiss,&
                      nfh, ifa, jheafa, ncomph)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: nnops
        integer(kind=8) :: nnop
        integer(kind=8) :: ddls
        integer(kind=8) :: ddlm
        integer(kind=8) :: ndim
        integer(kind=8) :: pla(27)
        real(kind=8) :: saut(3)
        real(kind=8) :: sautm(3)
        real(kind=8) :: nd(3)
        real(kind=8) :: ffc(16)
        real(kind=8) :: w11
        real(kind=8) :: w11m
        real(kind=8) :: jac
        real(kind=8) :: q1
        real(kind=8) :: dt
        real(kind=8) :: ta
        real(kind=8) :: q1m
        real(kind=8) :: ta1
        real(kind=8) :: q2
        real(kind=8) :: q2m
        real(kind=8) :: dffc(16,3)
        real(kind=8) :: rho11
        real(kind=8) :: gradpf(3)
        real(kind=8) :: rho11m
        real(kind=8) :: gradpfm(3)
        real(kind=8) :: ffp2(27)
        real(kind=8) :: vect(560)
        real(kind=8) :: ffp(27)
        real(kind=8) :: delta(6)
        real(kind=8) :: lamb(3)
        real(kind=8) :: am(3)
        real(kind=8) :: r
        real(kind=8) :: p(3,3)
        real(kind=8) :: psup
        real(kind=8) :: pinf
        real(kind=8) :: pf
        integer(kind=8) :: ncompn
        integer(kind=8) :: jheavn
        integer(kind=8) :: nfiss
        integer(kind=8) :: ifiss
        integer(kind=8) :: nfh
        integer(kind=8) :: ifa
        integer(kind=8) :: jheafa
        integer(kind=8) :: ncomph
    end subroutine xvechm
end interface
