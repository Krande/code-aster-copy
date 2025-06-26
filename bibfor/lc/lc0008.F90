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

subroutine lc0008(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  nvi, dsidep, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/lcmaza.h"
#include "asterfort/lcmazarsmu.h"
#include "asterfort/utmess.h"
!
! aslint: disable=W1504,W0104
!
    character(len=*), intent(in)    :: fami
    integer(kind=8), intent(in)             :: kpg
    integer(kind=8), intent(in)             :: ksp
    integer(kind=8), intent(in)             :: ndim
    integer(kind=8), intent(in)             :: imate
    character(len=16), intent(in)   :: compor(*)
    real(kind=8), intent(in)        :: carcri(*)
    real(kind=8), intent(in)        :: instam
    real(kind=8), intent(in)        :: instap
    real(kind=8), intent(in)        :: epsm(*)
    real(kind=8), intent(in)        :: deps(*)
    real(kind=8), intent(in)        :: sigm(*)
    real(kind=8), intent(in)        :: vim(*)
    character(len=16), intent(in)   :: option
    real(kind=8), intent(in)        :: angmas(*)
    character(len=8), intent(in)    :: typmod(*)
    integer(kind=8), intent(in)             :: icomp
    integer(kind=8), intent(in)             :: nvi
!
    real(kind=8), intent(out)       :: sigp(*)
    real(kind=8), intent(out)       :: vip(*)
    real(kind=8), intent(out)       :: dsidep(*)
    integer(kind=8), intent(out)            :: codret
!
! --------------------------------------------------------------------------------------------------
!
!                   Behaviour : MAZARS
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical    :: ismazarsmu, iscplane, iscouplage
    character(len=6) :: model
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    ismazarsmu = compor(1) (1:11) .eq. 'MAZARS_UNIL'
    iscplane = typmod(1) (1:6) .eq. 'C_PLAN'
!
    if (ismazarsmu) then
        iscouplage = (option(6:9) .eq. 'COUP')
        if (iscouplage) then
            call utmess('F', 'ALGORITH4_10', sk=compor(1))
        end if
        ! Loi de mazars unilatérale : Mu_Model
        if (iscplane) then
            model = 'C_PLAN'
            call lcmazarsmu(fami, kpg, ksp, ndim, imate, model, epsm, &
                            deps, vim, option, sigp, vip, dsidep)
        else
            model = '3D'
            call lcmazarsmu(fami, kpg, ksp, ndim, imate, model, epsm, &
                            deps, vim, option, sigp, vip, dsidep)
        end if
    else
        ! Loi de mazars "classique"
        call lcmaza(fami, kpg, ksp, ndim, typmod, &
                    imate, compor, epsm, deps, vim, &
                    option, sigp, vip, dsidep)
    end if
end subroutine
