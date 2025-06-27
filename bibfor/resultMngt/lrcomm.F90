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
subroutine lrcomm(lReuse, resultName, model, caraElem, fieldMate, lLireResu_, lVeriVari_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmdoco.h"
#include "asterfort/nmdorc.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/vrcomp.h"
#include "asterfort/vrcom2.h"
#include "asterfort/comp_info.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
    aster_logical, intent(in) :: lReuse
    character(len=8), intent(in) :: resultName
    character(len=8), intent(in):: model, fieldMate, caraElem
    aster_logical, intent(in), optional :: lLireResu_, lVeriVari_
!
! --------------------------------------------------------------------------------------------------
!
! LIRE_RESU and CREA_RESU
!
! Non-linear behaviour management
!
! --------------------------------------------------------------------------------------------------
!
! In  lReuse           : flag to re-use datastructure
! In  resultName       : name of results datastructure
! In  model            : name of model
! In  caraElem         : name of elementary characteristics field
! In  fieldMate        : name of material field
! In  lLireResu        : flag for LIRE_RESU command
! In  lVeriVari        : flag to check consistency of VARI_ELGA fields with COMPORTEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, nbOcc
    integer(kind=8) :: ifm, niv
    character(len=1) :: typeStop
    character(len=19) :: modelLigrel
    character(len=19) :: fieldToSave, fieldPrevious
    character(len=24) :: compor, carcri, variElga
    aster_logical :: lInitialState, lKeywfactCompor, lLireResu, verbose, lVeriVari, lModiVari
    integer(kind=8) :: storeNb, iStore, storeNume
    integer(kind=8), pointer :: storeList(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    verbose = niv .ge. 2
!
! - Initializations
!
    lInitialState = ASTER_FALSE
    lLireResu = ASTER_FALSE
    if (present(lLireResu_)) then
        lLireResu = lLireResu_
    end if
    lVeriVari = ASTER_TRUE
    if (present(lVeriVari_)) then
        lVeriVari = lVeriVari_
    end if
!
! - Name of constant fields
!
    compor = '&&LRCOMM.COMPOR'
    carcri = '&&LRCOMM.CARCRI'
!
! - detecter si COMPORTEMENT est renseigné
!
    lKeywfactCompor = ASTER_FALSE
    call getfac('COMPORTEMENT', nbOcc)
    if (nbOcc .gt. 0) then
        lKeywfactCompor = ASTER_TRUE
    end if
!
! - Get number of storing slots in results datastructure
!
    call rs_get_liststore(resultName, storeNb)
    if (storeNb .le. 0) then
        call utmess('F', 'RESULT2_97')
    end if
!
! - Get list of storing slots in results datastructure
!
    AS_ALLOCATE(vi=storeList, size=storeNb)
    call rs_get_liststore(resultName, storeNb, storeList)
!
! - Manage behaviour field
!
    if (model .ne. ' ') then
! ----- Read objects for constitutive laws
        call nmdorc(model, fieldMate, lInitialState, compor, carcri)
        if (niv .ge. 2) then
            call comp_info(model, compor)
        end if
        do iStore = 1, storeNb
            storeNume = storeList(iStore)
            call rsexch(' ', resultName, 'COMPORTEMENT', storeNume, fieldToSave, iret)
            if (lReuse .and. (.not. lKeywfactCompor)) then
! ------------- reuse sans COMPORTEMENT
                if (iret .eq. 100) then
                    if (storeNume .eq. 1) then
! --------------------- 1er instant : pas de précédent instant : => ELAS
                        call copisd('CHAMP_GD', 'G', compor, fieldToSave)
                    else
! --------------------- Get at previous storing index
                        call rsexch(' ', resultName, 'COMPORTEMENT', storeNume-1, &
                                    fieldPrevious, iret)
! --------------------- copier la carte COMPORTEMENT du précédent instant
                        call copisd('CHAMP_GD', 'G', fieldPrevious, fieldToSave)
                    end if
                    call rsnoch(resultName, 'COMPORTEMENT', storeNume)
                end if
            else
! ------------- Autres cas (sans reuse, reuse avec compor), ou lire_resu
                if (iret .le. 100) then
! ----------------- copier compor (ELAS ou celui donné par utilisateurs)
                    call copisd('CHAMP_GD', 'G', compor, fieldToSave)
                    call rsnoch(resultName, 'COMPORTEMENT', storeNume)
                end if
            end if
        end do
    end if
!
! - Check comportment
!
    do iStore = 1, storeNb
        storeNume = storeList(iStore)
! ----- Get internal state variables
        call rsexch(' ', resultName, 'VARI_ELGA', storeNume, variElga, iret)
        if (iret .eq. 0) then
            if (model .eq. ' ') then
                call utmess('F', 'RESULT2_17')
            end if
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
            if (lReuse .and. (.not. lKeywfactCompor)) then
                ! pour le cas reuse sans COMPORTEMENT
                ! vérifier entre VARI_ELGA existant et carte de comportement du resu
                call rsexch(' ', resultName, 'COMPORTEMENT', storeNume, compor, iret)
            end if
            call nmdoco(model, caraElem, compor)
            if (lVeriVari) then
! ------------- Check consistency of internal state variables with behaviour
                typeStop = 'A'
                call vrcomp(compor, variElga, modelLigrel, iret, typeStop_=typeStop, &
                            verbose_=verbose, lModiVari_=lModiVari)
                if (iret .eq. 1) then
                    call utmess(typeStop, 'RESULT1_1')
                end if
! ------------- Change internal variable field
                if (lModiVari) then
                    call vrcom2(compor, variElga, modelLigrel, lLireResu)
                end if
            end if
        end if
    end do
!
! - Clean
!
    AS_DEALLOCATE(vi=storeList)
!
    call jedema()
!
end subroutine
