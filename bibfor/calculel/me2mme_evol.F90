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
subroutine me2mme_evol(modelZ, caraElemZ, materFieldZ, matecoZ, nharm, jvBase, &
                       iLoad, loadName, ligrelCalcZ, timePrev, timeCurr, &
                       timeTheta, resuElem, vectElem)
!
    use loadMecaCompute_module
    use loadMecaCompute_type
    !
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, matecoZ
    integer(kind=8), intent(in) :: nharm
    character(len=1), intent(in) :: jvBase
    integer(kind=8), intent(in) :: iLoad
    character(len=8), intent(in) :: loadName
    character(len=*), intent(in) :: ligrelCalcZ
    real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
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
    integer(kind=8), parameter :: nbout = 1
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN), lpaout(nbout)
    character(len=19) :: lchin(LOAD_NEUM_NBMAXIN), lchout(nbout)
    aster_logical, parameter :: applySuiv = ASTER_FALSE
    integer(kind=8) :: nbFieldIn
    character(len=13) :: loadPreObject
    character(len=24) :: loadLigrel
    character(len=24) :: varcCurr
    character(len=24) :: dispPrev, dispCumuInst, strxPrev, viteCurr, acceCurr
    character(len=24) :: compor
!
! --------------------------------------------------------------------------------------------------
!
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "

! - Input fields (useless)
    dispPrev = " "
    dispCumuInst = " "
    strxPrev = " "
    viteCurr = " "
    acceCurr = " "
    compor = " "

! - Preparing input fields
    call prepGeneralFields(modelZ, caraElemZ, matecoZ, &
                           nharm, &
                           varcCurr, dispPrev, dispCumuInst, &
                           nbFieldIn, lpain, lchin)

! - Composite dead Neumann loads (EVOL_CHAR)
    loadPreObject = loadName(1:8)//'.CHME'
    loadLigrel = loadPreObject(1:13)//'.LIGRE'
    call compLoadEvolVect(modelZ, caraElemZ, materFieldZ, compor, &
                          timePrev, jvBase, &
                          applySuiv, iLoad, &
                          timePrev, timeCurr, timeTheta, &
                          loadPreObject, loadLigrel, ligrelCalcZ, &
                          nbFieldIn, lpain, lchin, &
                          resuElem, vectElem, &
                          dispPrev, dispCumuInst, strxPrev, viteCurr, acceCurr)
!
end subroutine
