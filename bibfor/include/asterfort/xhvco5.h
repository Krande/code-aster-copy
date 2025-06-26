! --------------------------------------------------------------------
! Copyright (C) 1991 - 2017 - EDF R&D - www.code-aster.org
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
    subroutine xhvco5(ndim, nnop, nnops, pla, nd,&
                      tau1, tau2, mu, nddls, jac,&
                      ffc, ffp, nddlm, wsaut,&
                      saut, vect, ifiss, nfiss, nfh,&
                      ifa, jheafa, ncomph, jheavn, ncompn, pf)
        integer(kind=8) :: ndim
        integer(kind=8) :: nnop
        integer(kind=8) :: nnops
        integer(kind=8) :: pla(27)
        real(kind=8) :: nd(3)
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        real(kind=8) :: mu(3)
        integer(kind=8) :: nddls
        real(kind=8) :: jac
        real(kind=8) :: ffc(16)
        real(kind=8) :: ffp(27)
        integer(kind=8) :: nddlm
        real(kind=8) :: wsaut(3)
        real(kind=8) :: saut(3)
        real(kind=8) :: vect(560)
        integer(kind=8) :: ifiss
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
        integer(kind=8) :: ifa
        integer(kind=8) :: jheafa
        integer(kind=8) :: ncomph
        integer(kind=8) :: jheavn
        integer(kind=8) :: ncompn
        real(kind=8) :: pf
    end subroutine xhvco5
end interface
