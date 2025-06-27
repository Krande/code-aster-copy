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
    subroutine xmfrot(algofr, coeffr, coeffp, ddlm, ddls,&
                      ffc, ffp, idepd, idepm, indco,&
                      jac, lact, mmat, mu, nd,&
                      ndim, nfh, nfiss, nno, nnol,&
                      nnos, nvit, pla, seuil,&
                      singu, fk, tau1, tau2)
        integer(kind=8) :: algofr
        real(kind=8) :: coeffr
        real(kind=8) :: coeffp
        integer(kind=8) :: ddlm
        integer(kind=8) :: ddls
        real(kind=8) :: ffc(8)
        real(kind=8) :: ffp(27)
        integer(kind=8) :: idepd
        integer(kind=8) :: idepm
        integer(kind=8) :: indco
        real(kind=8) :: jac
        integer(kind=8) :: lact(8)
        real(kind=8) :: mmat(216, 216)
        real(kind=8) :: mu
        real(kind=8) :: nd(3)
        integer(kind=8) :: ndim
        integer(kind=8) :: nfh
        integer(kind=8) :: nfiss
        integer(kind=8) :: nno
        integer(kind=8) :: nnol
        integer(kind=8) :: nnos
        integer(kind=8) :: nvit
        integer(kind=8) :: pla(27)
        real(kind=8) :: seuil
        integer(kind=8) :: singu
        real(kind=8) :: tau1(3)
        real(kind=8) :: tau2(3)
        real(kind=8) :: fk(27,3,3)
    end subroutine xmfrot
end interface
