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
subroutine mecgme(stop, &
                  modelZ, caraElemZ, materFieldZ, matecoZ, comporZ, &
                  listLoadZ, &
                  timePrev, timeCurr, &
                  dispPrevZ, dispCumuInstZ, &
                  matrElemZ, &
                  varcCurrZ_, nharm_, &
                  ligrelCalcZ_, jvBase_)
!
    use loadMecaCompute_module
    use loadMecaCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/load_list_info.h"
#include "asterfort/lxliis.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/wkvect.h"
#include "LoadTypes_type.h"
!
    character(len=1), intent(in) :: stop
    character(len=*), intent(in) :: modelZ, caraElemZ
    character(len=*), intent(in) :: materFieldZ, matecoZ, comporZ
    character(len=*), intent(in) :: listLoadZ
    real(kind=8), intent(in) :: timePrev, timeCurr
    character(len=*), intent(in) :: dispPrevZ, dispCumuInstZ
    character(len=*), intent(in) :: matrElemZ
    character(len=*), optional, intent(in) :: varcCurrZ_
    integer(kind=8), optional, intent(in) :: nharm_
    character(len=*), optional, intent(in) :: ligrelCalcZ_
    character(len=1), optional, intent(in) :: jvBase_
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Undead loads - Depending on geometry or speed - Matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  stop              : continue or stop computation if no loads on elements
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  materField        : name of material characteristics (field)
! In  mateco            : name of coded material
! In  compor            : name of non-linear behaviour definition (field)
! In  listLoad          : name of object for list of loads! In  timePrev          : previous time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! In  matrElem          : name of elementary matrix
! In  varcCurr          : external state variables for current time
! In  nharm             : Fourier mode
! In  ligrelCalc        : LIGREL to compute loads
! In  jvBase            : JEVEUX base to create objects
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: lpain(LOAD_NEUM_NBMAXIN)
    character(len=24) :: lchin(LOAD_NEUM_NBMAXIN)
    aster_logical, parameter :: applyPilo = ASTER_FALSE, applySuiv = ASTER_TRUE
    integer(kind=8) :: nbLoad, iLoad, loadNume, nbFieldIn
    integer(kind=8) :: nharm
    real(kind=8), parameter :: timeTheta = 0.d0
    aster_logical :: noLoadInList
    character(len=24), pointer :: listLoadName(:) => null(), listLoadFunc(:) => null()
    character(len=24) :: loadFuncJv
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=8) :: loadName, loadFunc
    character(len=13) :: loadPreObject
    character(len=24) :: matrCoefJv
    real(kind=8), pointer :: matrCoef(:) => null()
    real(kind=8) :: coefVale
    character(len=24) :: ligrelCalc, loadLigrel
    integer(kind=8) :: iret, ier, iResuElem, iMatr
    aster_logical :: l_first_matr, hasLoadFunc
    integer(kind=8) ::  nbResuElem
    character(len=24), pointer :: listResuElem(:) => null()
    character(len=1) :: jvBase
    character(len=24) :: varcCurr
    character(len=24) :: dispPrev, dispCumuInst, strxPrev, viteCurr, acceCurr
    character(len=24) :: compor
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
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
    compor = comporZ

! - Input fields (useless)
    strxPrev = " "
    viteCurr = " "
    acceCurr = " "

! - Get loads
    call load_list_info(noLoadInList, nbLoad, listLoadName, listLoadInfo, &
                        list_load_=listLoadZ)
    if (noLoadInList) then
        goto 99
    end if
    loadFuncJv = listLoadZ(1:19)//'.FCHA'

! - Allocate result
    call jeexin(matrElemZ(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        l_first_matr = ASTER_TRUE
        call memare(jvBase, matrElemZ, modelZ, 'CHAR_MECA')
        call reajre(matrElemZ, ' ', jvBase)
    else
        l_first_matr = ASTER_FALSE
        call jelira(matrElemZ(1:19)//'.RELR', 'LONUTI', nbResuElem)
        if (nbResuElem .gt. 0) then
            call jeveuo(matrElemZ(1:19)//'.RELR', 'L', vk24=listResuElem)
        end if
    end if

! - Preparing input fields
    call prepGeneralFields(modelZ, caraElemZ, matecoZ, &
                           nharm, &
                           varcCurr, dispPrev, dispCumuInst, &
                           nbFieldIn, lpain, lchin)
! - Computation
    if (l_first_matr) then
        do iLoad = 1, nbLoad
            iMatr = 0
            loadName = listLoadName(iLoad) (1:8)
            loadNume = listLoadInfo(nbLoad+iLoad+1)
            loadPreObject = loadName(1:8)//'.CHME'
            loadLigrel = loadPreObject(1:13)//'.LIGRE'

            if (loadNume .eq. 4) then
                call compLoadMatr(applyPilo, applySuiv, &
                                  iLoad, iMatr, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadLigrel, loadPreObject, &
                                  stop, nbFieldIn, lpain, lchin, &
                                  ligrelCalc, jvBase, matrElemZ)

                call compLoadEvolMatr(timeCurr, jvBase, &
                                      applySuiv, iLoad, iMatr, &
                                      modelZ, materFieldZ, compor, &
                                      strxPrev, viteCurr, acceCurr, &
                                      timePrev, timeCurr, timeTheta, &
                                      loadPreObject, loadLigrel, ligrelCalc, &
                                      nbFieldIn, lpain, lchin, &
                                      matrElemZ)

            end if
        end do
    else
        do iResuElem = 1, nbResuElem
            if (listResuElem(iResuElem) (10:10) .eq. 'G') then
                call lxliis(listResuElem(iResuElem) (7:8), iLoad, ier)
                iMatr = -iResuElem
                loadName = listLoadName(iLoad) (1:8)
                loadNume = listLoadInfo(nbLoad+iLoad+1)
                loadPreObject = loadName(1:8)//'.CHME'
                loadLigrel = loadPreObject(1:13)//'.LIGRE'

                if (loadNume .eq. 4) then
                    call compLoadMatr(applyPilo, applySuiv, &
                                      iLoad, iMatr, loadNume, &
                                      modelZ, materFieldZ, compor, &
                                      strxPrev, viteCurr, acceCurr, &
                                      timePrev, timeCurr, timeTheta, &
                                      loadLigrel, loadPreObject, &
                                      stop, nbFieldIn, lpain, lchin, &
                                      ligrelCalc, jvBase, matrElemZ)

                    call compLoadEvolMatr(timeCurr, jvBase, &
                                          applySuiv, iLoad, iMatr, &
                                          modelZ, materFieldZ, compor, &
                                          strxPrev, viteCurr, acceCurr, &
                                          timePrev, timeCurr, timeTheta, &
                                          loadPreObject, loadLigrel, ligrelCalc, &
                                          nbFieldIn, lpain, lchin, &
                                          matrElemZ)
                end if
            end if
        end do
    end if

! - Get number of resu_elem for undead loads
    call jeexin(matrElemZ(1:19)//'.RELR', iret)
    ASSERT(iret .ne. 0)
    call jelira(matrElemZ(1:19)//'.RELR', 'LONUTI', nbResuElem)
    if (nbResuElem .eq. 0) then
        goto 99
    else
        call jeveuo(matrElemZ(1:19)//'.RELR', 'L', vk24=listResuElem)
        if (listResuElem(1) (7:8) .eq. '00') then
            goto 99
        end if
    end if

! - Access to function
    matrCoefJv = matrElemZ(1:15)//'.COEF'
    call jeexin(loadFuncJv, iret)
    hasLoadFunc = ASTER_FALSE
    if (iret .gt. 0) then
        hasLoadFunc = .true.
        call jeveuo(loadFuncJv, 'L', vk24=listLoadFunc)
    end if

! - Create list of coefficients
    call jedetr(matrCoefJv)
    call wkvect(matrCoefJv, 'V V R', nbResuElem, vr=matrCoef)
    do iResuElem = 1, nbResuElem
        if (hasLoadFunc) then
            call lxliis(listResuElem(iResuElem) (7:8), iLoad, iret)
            loadFunc = listLoadFunc(iLoad) (1:8)
            ASSERT(iLoad .gt. 0)
            call fointe('F ', loadFunc, 1, ['INST'], [timeCurr], coefVale, iret)
        else
            coefVale = 1.d0
        end if
        matrCoef(iResuElem) = coefVale
    end do
!
99  continue
!
    call jedema()
end subroutine
