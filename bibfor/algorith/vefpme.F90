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
subroutine vefpme(stop, &
                  modelZ, caraElemZ, materFieldZ, matecoZ, &
                  loadNameJvZ, loadInfoJvZ, &
                  timePara, &
                  dispPrevZ, dispCumuInstZ, &
                  varcCurrZ, &
                  vectElemZ, ligrelCalcZ_, jvBase_)
!
    use loadMecaCompute_module
    use loadMecaCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/load_list_info.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
#include "LoadTypes_type.h"
!
    character(len=1), intent(in) :: stop
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, matecoZ
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    real(kind=8), intent(in) :: timePara(3)
    character(len=*), intent(in) :: dispPrevZ, dispCumuInstZ
    character(len=*), intent(in) :: varcCurrZ
    character(len=*), intent(inout) :: vectElemZ
    character(len=*), optional, intent(in) :: ligrelCalcZ_
    character(len=1), optional, intent(in) :: jvBase_
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Dead and "driven" loads
!
! --------------------------------------------------------------------------------------------------
!
! In  stop              : continue or stop computation if no loads on elements
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  mateco            : name of coded material
! In  caraElem          : name of elementary characteristics (field)
! In  loadNameJv        : name of object for list of loads name
! In  loadInfoJv        : name of object for list of loads info
! In  timePara          : times informations
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! In  varcCurr          : external state variables for current time
! IO  vectElem          : name of elementary vectors
! In  ligrelCalc        : LIGREL to compute loads
! In  jvBase            : JEVEUX base to create vector
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN)
    character(len=24) :: lchin(LOAD_NEUM_NBMAXIN)
    aster_logical, parameter :: applyPilo = ASTER_TRUE
    aster_logical :: applySuiv
    integer(kind=8) :: nbLoad, iLoad, loadNume, nbFieldInGene
    integer(kind=8) :: nharm
    real(kind=8) :: timePrev, timeCurr, timeTheta
    character(len=8) :: loadName, model
    character(len=13) :: loadPreObject
    character(len=24) :: ligrelCalc, loadLigrel
    character(len=24) :: vectElem, resuElem
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    aster_logical :: noLoadInList
    character(len=1) :: jvBase
    character(len=24) :: varcCurr
    character(len=24) :: dispPrev, dispCumuInst, strxPrev, viteCurr, acceCurr
    character(len=24) :: compor
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    resuElem = '&&VEFPME.0000000'
    model = modelZ(1:8)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrelCalc)
    if (present(ligrelCalcZ_)) then
        ligrelCalc = ligrelCalcZ_(1:24)
    end if
    jvBase = 'V'
    if (present(jvBase_)) then
        jvBase = jvBase_
    end if
    lpain = " "
    lchin = " "

! - Input fields (required)
    dispPrev = dispPrevZ
    dispCumuInst = dispCumuInstZ
    varcCurr = varcCurrZ

! - Input fields (useless)
    nharm = -1
    strxPrev = " "
    viteCurr = " "
    acceCurr = " "
    compor = " "

! - Time stepping
    timePrev = timePara(1)
    timeCurr = timePara(1)+timePara(2)
    timeTheta = timePara(3)

! - Result name for vectElem
    vectElem = vectElemZ
    if (vectElem .eq. ' ') then
        vectElem = '&&VEFPME'
    end if

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJvZ, loadInfoJvZ)

! - Allocate result
    call detrsd('VECT_ELEM', vectElem)
    call vemare(jvBase, vectElem, model)
    call reajre(vectElem, ' ', jvBase)
    if (noLoadInList) then
        goto 99
    end if

! - Preparing input fields
    call prepGeneralFields(model, caraElemZ, matecoZ, &
                           nharm, &
                           varcCurr, dispPrev, dispCumuInst, &
                           nbFieldInGene, lpain, lchin)

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)
        loadPreObject = loadName(1:8)//'.CHME'
        loadLigrel = loadPreObject(1:13)//'.LIGRE'

! ----- Standard dead Neumann loads
        if ((loadNume .eq. 5) .or. (loadNume .eq. 8)) then
            applySuiv = ASTER_FALSE
            call compLoadVect(applyPilo, applySuiv, &
                              iLoad, loadNume, &
                              model, materFieldZ, compor, &
                              strxPrev, viteCurr, acceCurr, &
                              timePrev, timeCurr, timeTheta, &
                              loadPreObject, loadLigrel, &
                              stop, nbFieldInGene, lpain, lchin, &
                              ligrelCalc, jvBase, resuElem, vectElem)
        end if

! ----- Standard undead Neumann loads
        if ((loadNume .eq. 9) .or. (loadNume .eq. 11)) then
            applySuiv = ASTER_TRUE
            call compLoadVect(applyPilo, applySuiv, &
                              iLoad, loadNume, &
                              model, materFieldZ, compor, &
                              strxPrev, viteCurr, acceCurr, &
                              timePrev, timeCurr, timeTheta, &
                              loadPreObject, loadLigrel, &
                              stop, nbFieldInGene, lpain, lchin, &
                              ligrelCalc, jvBase, resuElem, vectElem)
        end if
    end do
!
99  continue
!
    vectElemZ = vectElem(1:19)//'.RELR'
!
    call jedema()
end subroutine
