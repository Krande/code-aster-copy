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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmgvno(fami, ndim, nno1, nno2, npg, &
                      iw, vff1, vff2, idfde1, idfde2, &
                      geom, typmod, option, mat, compor, &
                      lgpg, carcri, instam, instap, ddlm, &
                      ddld, angmas, sigm, vim, sigp, &
                      vip, matr, vect, codret)
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        character(len=*) :: fami
        integer(kind=8) :: iw
        real(kind=8) :: vff1(nno1, npg)
        real(kind=8) :: vff2(nno2, npg)
        integer(kind=8) :: idfde1
        integer(kind=8) :: idfde2
        real(kind=8) :: geom(ndim, nno1)
        character(len=8) :: typmod(2)
        character(len=16) :: option
        integer(kind=8) :: mat
        character(len=16) :: compor(COMPOR_SIZE)
        real(kind=8) :: carcri(CARCRI_SIZE)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: ddlm(*)
        real(kind=8) :: ddld(*)
        real(kind=8) :: angmas(3)
        real(kind=8) :: sigm(2*ndim+1, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigp(2*ndim+1, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: matr(*)
        real(kind=8) :: vect(*)
        integer(kind=8) :: codret
    end subroutine nmgvno
end interface
