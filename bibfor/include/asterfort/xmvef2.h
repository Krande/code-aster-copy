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
    subroutine xmvef2(ndim, nno, nnos, ffp, jac,&
                      seuil, reac12, singu, fk, nfh,&
                      coeffp, coeffr, mu, algofr, nd,&
                      ddls, ddlm, idepl, pb, vtmp)
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        integer(kind=8) :: nnos
        real(kind=8) :: ffp(27)
        real(kind=8) :: jac
        real(kind=8) :: seuil
        real(kind=8) :: reac12(3)
        integer(kind=8) :: singu
        integer(kind=8) :: nfh
        real(kind=8) :: coeffp
        real(kind=8) :: coeffr
        real(kind=8) :: mu
        integer(kind=8) :: algofr
        real(kind=8) :: nd(3)
        integer(kind=8) :: ddls
        integer(kind=8) :: ddlm
        integer(kind=8) :: idepl
        real(kind=8) :: pb(3)
        real(kind=8) :: vtmp(400)
        real(kind=8) :: fk(27,3,3)
    end subroutine xmvef2
end interface
