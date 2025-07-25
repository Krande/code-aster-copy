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
#include "asterf_types.h"
!
interface
    subroutine nifilg(ndim, nnod, nnog, nnop, npg,&
                      iw, vffd, vffg, vffp, idffd,&
                      vu, vg, vp, geomi, typmod,&
                      option, mate, compor, lgpg, carcri,&
                      instm, instp, ddlm, ddld, angmas,&
                      sigm, vim, sigp, vip,&
                      lMatr, lVect, lSigm, lVari,&
                      vect, matr, matsym, codret)
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nnop
        integer(kind=8) :: nnog
        integer(kind=8) :: nnod
        integer(kind=8) :: ndim
        integer(kind=8) :: iw
        real(kind=8) :: vffd(nnod, npg)
        real(kind=8) :: vffg(nnog, npg)
        real(kind=8) :: vffp(nnop, npg)
        integer(kind=8) :: idffd
        integer(kind=8) :: vu(3, 27)
        integer(kind=8) :: vg(27)
        integer(kind=8) :: vp(27)
        real(kind=8) :: geomi(ndim, nnod)
        character(len=8) :: typmod(*)
        character(len=16) :: option
        integer(kind=8) :: mate
        character(len=16) :: compor(*)
        real(kind=8) :: carcri(*)
        real(kind=8) :: instm
        real(kind=8) :: instp
        real(kind=8) :: ddlm(*)
        real(kind=8) :: ddld(*)
        real(kind=8) :: angmas(*)
        real(kind=8) :: sigm(2*ndim+1, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigp(2*ndim+1, npg)
        real(kind=8) :: vip(lgpg, npg)
        aster_logical :: resi
        aster_logical :: rigi
        real(kind=8) :: vect(*)
        real(kind=8) :: matr(*)
        aster_logical :: matsym
        aster_logical, intent(in) :: lMatr, lVect, lSigm, lVari
        integer(kind=8) :: codret
    end subroutine nifilg
end interface
