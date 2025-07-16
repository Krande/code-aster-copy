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
subroutine merxth(l_stat, &
                  modelZ, caraElemZ, matecoZ, &
                  loadNameJvZ, loadInfoJvZ, &
                  tpsthe, timeMapZ, &
                  tempIterZ, comporTherZ, varcCurrZ, &
                  matrElemZ, jvBase)
!
    use loadTherCompute_module
    use loadTherCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/gcnco2.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/load_list_info.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/ther_mtan.h"
#include "LoadTypes_type.h"
!
    aster_logical, intent(in) :: l_stat
    character(len=*), intent(in) :: modelZ, caraElemZ, matecoZ
    character(len=*), intent(in) :: loadNameJvZ, loadInfoJvZ
    real(kind=8), intent(in) :: tpsthe(6)
    character(len=*), intent(in) :: comporTherZ, timeMapZ
    character(len=*), intent(in) :: tempIterZ, varcCurrZ
    character(len=*), intent(in) :: matrElemZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Tangent matrix (non-linear) - Material and loads
!
! --------------------------------------------------------------------------------------------------
!
! In  l_stat            : .true. if stationnary
! In  model             : name of the model
! In  caraElem          : name of elementary characteristics (field)
! In  mateco            : name of coded material
! In  loadNameJv        : name of object for list of loads name
! In  loadInfoJv        : name of object for list of loads info
! In  tpsthe            : parameters for time
! In  timeMap           : time (<CARTE>)
! In  tempIter          : temperature field at current Newton iteration
! In  comporTher        : name of comportment definition (field)
! In  varcCurr          : command variable for current time
! In  matrElem          : name of matrElem result
! In  jvBase            : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: lpain(LOAD_NEUT_NBMAXIN)
    character(len=24) :: lchin(LOAD_NEUT_NBMAXIN)
    integer(kind=8) :: nbLoad, iLoad, loadNume, nbFieldInGene
    character(len=8) :: loadName, newnom
    aster_logical :: noLoadInList
    character(len=24), pointer :: listLoadName(:) => null()
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24) :: matrElem, resuElem
    real(kind=8) :: timeCurr, timePara(2), theta
    character(len=24) :: timeMap, tempPrev, tempIter, varcCurr
    character(len=13) :: loadPreObject
    character(len=24) :: loadLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lchin = " "

! - Get fields
    timeMap = timeMapZ
    tempIter = tempIterZ
    varcCurr = varcCurrZ
    tempPrev = " "

! - Get time parameters
    timeCurr = tpsthe(1)
    timePara(1) = tpsthe(3)
    timePara(2) = tpsthe(2)
    theta = timePara(1)

! - Name of elementary matrices
    matrElem = matrElemZ
    if (matrElem .eq. ' ') then
        matrElem = '&&MERXTH'
    end if

! - Allocate result
    call detrsd('MATR_ELEM', matrElem)
    call memare(jvBase, matrElem, modelZ, 'MTAN_THER')
    call reajre(matrElem, ' ', jvBase)

! - Generate new RESU_ELEM name
    resuElem = matrElem(1:8)//'.0000000'
    newnom = resuElem(9:16)
    call gcnco2(newnom)
    resuElem(10:16) = newnom(2:8)

! - Tangent matrix - Compute material part (rigidity and mass)
    call ther_mtan(l_stat, &
                   modelZ, caraElemZ, matecoZ, &
                   timePara, varcCurr, &
                   comporTherZ, tempIter, resuElem, &
                   matrElem, jvBase)

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        loadNameJvZ, loadInfoJvZ)

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
            call compLoadMatr(l_stat, theta, &
                              modelZ, timeMap, &
                              loadNume, &
                              loadPreObject, loadLigrel, &
                              nbFieldInGene, lpain, lchin, &
                              jvBase, resuElem, matrElem)

! --------- Composite Neumann loads (EVOL_CHAR)
            call compLoadEvolMatr(l_stat, theta, timeCurr, &
                                  modelZ, timeMap, &
                                  loadPreObject, loadLigrel, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElem, matrElem)
        end if
    end do

    call jedema()
!
end subroutine
