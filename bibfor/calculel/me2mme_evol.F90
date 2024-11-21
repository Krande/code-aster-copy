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
subroutine me2mme_evol(modelZ, caraElemZ, mateZ, matecoZ, nharm, jvBase, &
                       iLoad, loadName, ligrel_calcZ, inst_prev, inst_curr, &
                       inst_theta, resuElem, vectElem)
!
    use loadCompute_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/load_neum_prep.h"
#include "LoadTypes_type.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=*), intent(in) :: caraElemZ
    character(len=*), intent(in) :: mateZ
    character(len=*), intent(in) :: matecoZ
    integer, intent(in) :: nharm
    character(len=1), intent(in) :: jvBase
    integer, intent(in) :: iLoad
    character(len=8), intent(in) :: loadName
    character(len=*), intent(in) :: ligrel_calcZ
    real(kind=8), intent(in) :: inst_prev
    real(kind=8), intent(in) :: inst_curr
    real(kind=8), intent(in) :: inst_theta
    character(len=19), intent(inout) :: resuElem
    character(len=19), intent(in) :: vectElem
!
! --------------------------------------------------------------------------------------------------
!
! CALC_VECT_ELEM
!
! EVOL_CHAR loads
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbout = 1
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN), lpaout(nbout)
    character(len=19) :: lchin(LOAD_NEUM_NBMAXIN), lchout(nbout)
    character(len=4), parameter :: loadApply = "Dead"
    integer :: nb_in_prep
    character(len=24) :: model, caraElem, mate, mateco, ligrel_calc
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    caraElem = caraElemZ
    mate = mateZ
    mateco = matecoZ
    ligrel_calc = ligrel_calcZ
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Preparing input fields
    call load_neum_prep(model, caraElem, mate, mateco, 'Dead', inst_prev, &
                        inst_curr, inst_theta, LOAD_NEUM_NBMAXIN, nb_in_prep, lchin, &
                        lpain, nharm=nharm)

! - Composite dead Neumann loads (EVOL_CHAR)
    call compEvolChar(model, caraElem, inst_prev, jvBase, &
                      iLoad, loadName, loadApply, ligrel_calc, &
                      nb_in_prep, lpain, lchin, &
                      resuElem, vectElem)
!
end subroutine
