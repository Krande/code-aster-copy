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
subroutine vechth(typeTher, &
                  modelZ, matecoZ, &
                  loadNameJvZ, loadInfoJvZ, &
                  timeCurr, &
                  vectElemZ, &
                  varcCurrZ_, timeMapZ_, tempPrevZ_, timeMoveZ_, &
                  jvBase_)
!
    use loadTherCompute_module
    use loadTherCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/load_list_info.h"
#include "asterfort/vemare.h"
#include "asterfort/reajre.h"
#include "LoadTypes_type.h"
!
    character(len=4), intent(in) :: typeTher
    character(len=*), intent(in) :: modelZ, matecoZ
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    real(kind=8), intent(in) :: timeCurr
    character(len=*), intent(inout) :: vectElemZ
    character(len=*), optional, intent(in) :: varcCurrZ_, timeMapZ_, tempPrevZ_, timeMoveZ_
    character(len=1), optional, intent(in) :: jvBase_
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads (thermic)
!
! Neumann loads elementary vectors (second member)
!
! --------------------------------------------------------------------------------------------------
!
! In  typeTher          : type of thermics
!                         'MOVE' for moving sources
!                         'STAT' if not
! In  model             : name of the model
! In  mateco            : name of coded material
! In  loadNameJv        : name of object for list of loads name
! In  loadInfoJv        : name of object for list of loads info
! In  timeCurr          : current time
! In  timeMap           : time (<CARTE>)
! In  varcCurr          : command variable for current time
! In  tempPrev          : previous temperature
! In  timeMove          : modified time (<CARTE>) for THER_NON_LINE_MO
! IO  vectElem          : name of vectElem result
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: lpain(LOAD_NEUT_NBMAXIN)
    character(len=24) :: lchin(LOAD_NEUT_NBMAXIN)
    integer(kind=8) :: nbLoad, iLoad, loadNume, nbFieldInGene
    character(len=8) :: loadName
    aster_logical :: noLoadInList
    character(len=1) :: jvBase
    character(len=24) :: vectElem, resuElem
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24) :: timeMap, tempPrev, timeMove, tempIter, varcCurr
    character(len=13) :: loadPreObject
    character(len=24) :: loadLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    resuElem = '&&VECHTH.0000000'
    jvBase = 'V'
    if (present(jvBase_)) then
        jvBase = jvBase_
    end if
    lpain = " "
    lchin = " "

! - Get fields
    timeMap = " "
    if (present(timeMapZ_)) then
        timeMap = timeMapZ_
    end if
    tempPrev = " "
    if (present(tempPrevZ_)) then
        tempPrev = tempPrevZ_
    end if
    timeMove = ' '
    if (present(timeMoveZ_)) then
        ASSERT(typeTher .eq. 'MOVE')
        timeMove = timeMoveZ_
    end if
    varcCurr = ' '
    if (present(varcCurrZ_)) then
        varcCurr = varcCurrZ_
    end if
    tempIter = " "

! - Name of elementary vectors
    vectElem = vectElemZ
    if (vectElem .eq. ' ') then
        vectElem = '&&VECHTH'
    end if

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJvZ, loadInfoJvZ)

! - Allocate result
    call vemare(jvBase, vectElem, modelZ)
    call reajre(vectElem, ' ', jvBase)
    if (noLoadInList) then
        goto 99
    end if

! - Preparing input fields
    call prepGeneralFields(modelZ, matecoZ, &
                           varcCurr, tempPrev, tempIter, &
                           nbFieldInGene, lpain, lchin)

! - Computation
    do iLoad = 1, nbLoad
        loadName = listLoadName(iLoad) (1:8)
        loadNume = listLoadInfo(nbLoad+iLoad+1)
        loadPreObject = loadName(1:8)//'.CHTH'
        loadLigrel = loadPreObject(1:13)//'.LIGRE'

        if (loadNume .gt. 0) then
! --------- Standard Neumann loads
            call compLoadVect(typeTher, &
                              modelZ, timeMap, timeMove, &
                              iLoad, loadNume, &
                              loadPreObject, loadLigrel, &
                              nbFieldInGene, lpain, lchin, &
                              jvBase, resuElem, vectElem)

! --------- Composite Neumann loads (EVOL_CHAR)
            call compLoadEvolVect(typeTher, &
                                  timeCurr, modelZ, timeMap, timeMove, &
                                  iLoad, loadPreObject, loadLigrel, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElem, vectElem)
        end if
    end do
!
99  continue
!
    vectElemZ = vectElem(1:19)//'.RELR'
!
    call jedema()
end subroutine
