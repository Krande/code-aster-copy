! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
    subroutine ngfint(option, typmod, ndim, nddl, neps, &
                      npg, w, b, compor, fami, &
                      mat, angmas, lgpg, carcri, instam, &
                      instap, ddlm, ddld, ni2ldc, sigmam, &
                      vim, sigmap, vip, fint, matr, &
                      lMatr, lVect, lSigm, &
                      codret)
        character(len=8)  :: typmod(2)
        character(len=*)  :: fami
        character(len=16) :: option, compor(COMPOR_SIZE)
        integer :: ndim, nddl, neps, npg, mat, lgpg
        real(kind=8) :: w(neps, npg), ni2ldc(neps, npg), b(neps, npg, nddl)
        real(kind=8) :: angmas(3), carcri(CARCRI_SIZE), instam, instap
        real(kind=8) :: ddlm(nddl), ddld(nddl)
        real(kind=8) :: sigmam(neps, npg), sigmap(neps, npg)
        real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), matr(nddl, nddl), fint(nddl)
        aster_logical, intent(in)    :: lMatr, lVect, lSigm
        integer, intent(out)         :: codret
    end subroutine ngfint
end interface
