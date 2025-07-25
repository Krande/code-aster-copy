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
    subroutine nglgic(fami, option, typmod, ndim, nno, &
                      nnob, npg, nddl, iw, vff, &
                      vffb, idff, idffb, geomi, compor, &
                      mate, lgpg, carcri, angmas, instm, &
                      instp, matsym, ddlm, ddld, siefm, &
                      vim, siefp, vip, fint, matr, &
                      lMatr, lVect, lSigm, lVari, &
                      codret)
        aster_logical, intent(in)    :: matsym
        character(len=8), intent(in) :: typmod(2)
        character(len=*), intent(in) :: fami
        character(len=16), intent(in):: option, compor(COMPOR_SIZE)
        integer(kind=8), intent(in)          :: ndim, nno, nnob, npg, nddl, lgpg
        integer(kind=8), intent(in)          :: mate, iw, idff, idffb
        real(kind=8), intent(in)     :: geomi(ndim, nno), carcri(CARCRI_SIZE), instm, instp
        real(kind=8), intent(in)     :: vff(nno, npg), vffb(nnob, npg)
        real(kind=8), intent(in)     :: angmas(3), ddlm(nddl), ddld(nddl), siefm(3*ndim+4, npg)
        real(kind=8), intent(in)     :: vim(lgpg, npg)
  real(kind=8), intent(out)    :: fint(nddl), matr(nddl, nddl), siefp(3*ndim+4, npg), vip(lgpg, npg)
        aster_logical, intent(in)    :: lMatr, lVect, lSigm, lVari
        integer(kind=8), intent(out)         :: codret
    end subroutine nglgic
end interface
