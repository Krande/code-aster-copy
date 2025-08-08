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
subroutine vecgme(stop, &
                  modelZ, caraElemZ, materFieldZ, matecoZ, comporZ, &
                  loadNameJvZ, loadInfoJvZ, &
                  timePrev, timeCurr, &
                  dispPrevZ, dispCumuInstZ, &
                  viteCurrZ, acceCurrZ, strxPrevZ, &
                  vectElemZ, &
                  varcCurrZ_, nharm_, &
                  ligrelCalcZ_, jvBase_)
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
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, matecoZ, comporZ
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    real(kind=8), intent(in) :: timePrev, timeCurr
    character(len=*), intent(in) :: dispPrevZ, dispCumuInstZ
    character(len=*), intent(in) :: viteCurrZ, acceCurrZ, strxPrevZ
    character(len=*), intent(inout) :: vectElemZ
    character(len=*), optional, intent(in) :: varcCurrZ_
    integer(kind=8), optional, intent(in) :: nharm_
    character(len=*), optional, intent(in) :: ligrelCalcZ_
    character(len=1), optional, intent(in) :: jvBase_
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Undead loads - Depending on geometry or speed - Vector
!
! --------------------------------------------------------------------------------------------------
!
! In  stop              : continue or stop computation if no loads on elements
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  materField        : name of material characteristics (field)
! In  mateco            : name of coded material
! In  compor            : name of non-linear behaviour definition (field)
! In  loadNameJv        : name of object for list of loads name
! In  loadInfoJv        : name of object for list of loads info
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  strxPrev          : fibers information at beginning of current time
! IO  vectElem          : name of elementary vectors
! In  varcCurr          : external state variables for current time
! In  nharm_            : Fourier mode
! In  ligrelCalc        : LIGREL to compute loads
! In  jvBase            : JEVEUX base to create vector
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN)
    character(len=24) :: lchin(LOAD_NEUM_NBMAXIN)
    aster_logical, parameter :: applyPilo = ASTER_FALSE, applySuiv = ASTER_TRUE
    integer(kind=8) :: nbLoad, iLoad, loadNume, nbFieldInGene
    integer(kind=8) :: nharm
    real(kind=8), parameter :: timeTheta = 0.d0
    character(len=8) :: loadName
    character(len=13) :: loadPreObject
    character(len=24) :: ligrelCalc, loadLigrel
    character(len=19) :: vectElem, resuElem
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
    resuElem = '&&VECGME.0000000'
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=ligrelCalc)
    if (present(ligrelCalcZ_)) then
        ligrelCalc = ligrelCalcZ_
    end if
    jvBase = 'V'
    if (present(jvBase_)) then
        jvBase = jvBase_
    end if
    lpain = " "
    lchin = " "

! - Input fields (optional)
    varcCurr = ' '
    if (present(varcCurrZ_)) then
        varcCurr = varcCurrZ_
    end if
    nharm = -1
    if (present(nharm_)) then
        nharm = nharm_
    end if

! - Input fields (required)
    dispPrev = dispPrevZ
    dispCumuInst = dispCumuInstZ
    strxPrev = strxPrevZ
    viteCurr = viteCurrZ
    acceCurr = acceCurrZ
    compor = comporZ

! - Result name for vectElem
    vectElem = vectElemZ
    if (vectElem .eq. ' ') then
        vectElem = '&&VECGME'
    end if

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJvZ, loadInfoJvZ)

! - Allocate result
    call detrsd('VECT_ELEM', vectElem)
    call vemare(jvBase, vectElem, modelZ)
    call reajre(vectElem, ' ', jvBase)
    if (noLoadInList) then
        goto 99
    end if

! - Preparing input fields
    call prepGeneralFields(modelZ, caraElemZ, matecoZ, &
                           nharm, &
                           varcCurr, dispPrev, dispCumuInst, &
                           nbFieldInGene, lpain, lchin)

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)
        loadPreObject = loadName(1:8)//'.CHME'
        loadLigrel = loadPreObject(1:13)//'.LIGRE'

        if (loadNume .eq. 4) then
! --------- Standard dead Neumann loads
            call compLoadVect(applyPilo, applySuiv, &
                              iLoad, loadNume, &
                              modelZ, materFieldZ, compor, &
                              strxPrev, viteCurr, acceCurr, &
                              timePrev, timeCurr, timeTheta, &
                              loadPreObject, loadLigrel, &
                              stop, nbFieldInGene, lpain, lchin, &
                              ligrelCalc, jvBase, resuElem, vectElem)

! --------- Composite dead Neumann loads (EVOL_CHAR)
            call compLoadEvolVect(modelZ, caraElemZ, materFieldZ, compor, &
                                  timeCurr, jvBase, &
                                  applySuiv, iLoad, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadPreObject, loadLigrel, ligrelCalc, &
                                  nbFieldInGene, lpain, lchin, &
                                  resuElem, vectElem, &
                                  dispPrev, dispCumuInst, &
                                  strxPrev, viteCurr, acceCurr)
        end if
    end do
!
99  continue
!
    vectElemZ = vectElem//'.RELR'
!
    call jedema()
end subroutine
