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
    subroutine sigvmc(fami, nno, ndim, nbsig, npg, &
                      jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                      nodeCoor, nodeDisp, &
                      time, anglNaut, jvMaterCode, nharm, &
                      sigm)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: nno, ndim, nbsig, npg
        integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
        real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
        real(kind=8), intent(in) :: time, anglNaut(3)
        integer(kind=8), intent(in) :: jvMaterCode
        real(kind=8), intent(in) :: nharm
        real(kind=8), intent(out) :: sigm(nbsig*npg)
    end subroutine sigvmc
end interface
