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
#include "asterf_types.h"
!
interface
    subroutine xmmsa1(algofr, ndim, nno, nnos, nnol,&
                      pla, ffc, ffp, idepd, idepm,&
                      nfh, nd, tau1, tau2, singu,&
                      fk, lact, ddls, ddlm, coeffr,&
                      coeffp, p, adher, knp, ptknp,&
                      ik)
        integer(kind=8) :: algofr
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        integer(kind=8) :: nnos
        integer(kind=8) :: nnol
        integer(kind=8) :: pla(27)
        real(kind=8) :: ffc(8)
        real(kind=8) :: ffp(27)
        integer(kind=8) :: idepd
        integer(kind=8) :: idepm
        integer(kind=8) :: nfh
        real(kind=8) :: nd(3)
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        integer(kind=8) :: singu
        real(kind=8) :: rr
        integer(kind=8) :: lact(8)
        integer(kind=8) :: ddls
        integer(kind=8) :: ddlm
        real(kind=8) :: coeffr
        real(kind=8) :: coeffp
        real(kind=8) :: p(3, 3)
        aster_logical :: adher
        real(kind=8) :: knp(3, 3)
        real(kind=8) :: ptknp(3, 3)
        real(kind=8) :: ik(3, 3)
        real(kind=8) :: fk(27,3,3)
    end subroutine xmmsa1
end interface
