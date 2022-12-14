! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine nmgrla(option  , typmod  ,&
                      fami    , imate   ,&
                      ndim    , nno     , npg     , lgpg     ,&
                      ipoids  , ivf     , vff      , idfde,&
                      compor  , carcri  , mult_comp,&
                      instam  , instap  ,&
                      geomInit, dispPrev, dispIncr,&
                      angmas  , sigmPrev, sigmCurr,&
                      vim     , vip     ,&
                      matsym  , matuu   , vectu   ,&
                      codret)
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(*)
        character(len=*), intent(in) :: fami
        integer, intent(in) :: imate
        integer, intent(in) :: ndim, nno, npg, lgpg
        integer, intent(in) :: ipoids, ivf, idfde
        real(kind=8), intent(in) :: vff(*)
        character(len=16), intent(in) :: compor(*)
        real(kind=8), intent(in) :: carcri(*)
        character(len=16), intent(in) :: mult_comp
        real(kind=8), intent(in) :: instam, instap, angmas(*)
        real(kind=8), intent(in) :: geomInit(ndim, nno)
        real(kind=8), intent(inout) :: dispPrev(ndim*nno),  dispIncr(ndim*nno)
        real(kind=8), intent(inout) :: sigmPrev(2*ndim, npg), sigmCurr(2*ndim, npg)
        real(kind=8), intent(inout) :: vim(lgpg, npg), vip(lgpg, npg)
        aster_logical, intent(in) :: matsym
        real(kind=8), intent(inout) :: matuu(*)
        real(kind=8), intent(inout) :: vectu(ndim*nno)
        integer, intent(inout) :: codret
    end subroutine nmgrla
end interface
