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
subroutine me2mme(modelZ, listLoadZ, &
                  materFieldZ, matecoZ, caraElemZ, &
                  timeCurr, vectElemZ, numeHarm, jvBaseZ)
!
    use ExternalStateVariableComp_module
    use listLoad_module
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/meharm.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/vechme.h"
#include "asterfort/vedime.h"
#include "LoadTypes_type.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: modelZ, caraElemZ, vectElemZ, listLoadZ, materFieldZ, matecoZ, jvBaseZ
    real(kind=8) :: timeCurr
    integer(kind=8) :: numeHarm
!
! --------------------------------------------------------------------------------------------------
!
! Compute RHS for macro-element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: jvBase
    character(len=2) :: codret
    character(len=8) :: model, caraElem
    character(len=24) :: vectElem
    character(len=24) :: modelLigrel, chtime, listLoad
    character(len=24) :: loadInfoJv, loadNameJv
    character(len=24) :: varcCurr, varcRefe
    real(kind=8) :: timePara(3)
    aster_logical, parameter :: varcLine = ASTER_TRUE
    integer(kind=8) :: nbFaceVite, nbPlaneWave, nbWave
    character(len=8), pointer :: planeWave(:) => null()
    character(len=8), pointer :: faceVite(:) => null()
    character(len=8), pointer :: wave(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelZ
    caraElem = caraElemZ
    vectElem = vectElemZ
    jvBase = jvBaseZ
    listLoad = listLoadZ
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)

! - Access to loads
    loadNameJv = listLoad(1:19)//".LCHA"
    loadInfoJv = listLoad(1:19)//".INFC"

! - Special loads
    call detectMecaNeumLoad(listLoad, LOAD_NEUM_WAVE, &
                            nbWave, wave)
    call detectMecaNeumLoad(listLoad, LOAD_NEUM_PWAVE, &
                            nbPlaneWave, planeWave)
    call detectMecaNeumLoad(listLoad, LOAD_NEUM_VITE_FACE, &
                            nbFaceVite, faceVite)

    if (nbWave .ne. 0 .or. nbPlaneWave .ne. 0 .or. nbFaceVite .ne. 0) then
        call utmess("F", "SOUSTRUC_19")
    end if

! - Prepare time map
    timePara = 0.d0
    timePara(1) = timeCurr
    chtime = '&&ME2MME.CH_INST_R'
    call mecact('V', chtime, 'MODELE', modelLigrel, 'INST_R  ', &
                ncmp=1, nomcmp='INST', sr=timeCurr)

! - Prepare vect_elem
    call jedetr(vectElem(1:19)//'.RELR')

! - Compute external state variables
    varcCurr = '&&ME2MME.VARC'
    varcRefe = '&&ME2MME.VARC.REF'
    call vrcins(model, materFieldZ, caraElem, timeCurr, varcCurr, codret)
    call vrcref(model, materFieldZ, caraElem, varcRefe)

! - Compute elementary vectors for Neumann loads
    call vechme("S", &
                modelZ, caraElemZ, materFieldZ, matecoZ, &
                loadNameJv, loadInfoJv, &
                timePara, &
                vectElem, &
                varcCurr)

! - Compute elementary vectors for Dirichlet loads
    call vedime(model, loadNameJv, loadInfoJv, timeCurr, "R", vectElem, &
                lCumul_=ASTER_TRUE)

! - Compute elementary vectors for external state variables
    call varcCompElem(varcLine, &
                      numeHarm, modelZ, caraElemZ, materFieldZ, matecoZ, &
                      chtime, varcRefe, varcCurr, &
                      jvBase, vectElem, &
                      lCumul_=ASTER_TRUE)
!
    call detrsd('CHAMP_GD', varcCurr)
    call detrsd('CHAMP_GD', varcRefe)
    AS_DEALLOCATE(vk8=planeWave)
    AS_DEALLOCATE(vk8=wave)
    AS_DEALLOCATE(vk8=faceVite)
!
    call jedema()
end subroutine
