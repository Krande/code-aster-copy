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
subroutine ccvepo(option, model, &
                  postCompResu, &
                  postCompPoux)
!
    use listLoad_module
    use FED_module
    use mesh_module
    use postComp_type
    use postComp_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cochre.h"
#include "asterfort/copich.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medome_once.h"
#include "asterfort/rs_get_listload.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: model
    type(POST_COMP_RESU), intent(in) :: postCompResu
    type(POST_COMP_POUX), intent(out) :: postCompPoux
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  Specific checks for POUX beams
!
! --------------------------------------------------------------------------------------------------
!
!  ROUTINE PERMETTANT DE SAVOIR SI DES POUTRES SONT DANS LE LIGREL
!  REDUIT ET DE VERIFIER LES CHARGES REPARTIES
!
! In  option            : option to compute
! In  model             : model
! In  resultIn          : name of datastructure for input results
! In  resultType        : type of results datastructure
! In  nbStore           : number of storing indexes
! In  listStore         : list of storing indexes
! Out postCompPoux      : datastructure to manage loads for POUX beams
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    character(len=16), parameter :: factorKeyword = " "
    integer(kind=8), parameter :: iOccZero = 0
    integer(kind=8) :: ierd, jvPara, numeStore0
    character(len=8) :: answer, modelOnce, mesh
    character(len=16) :: modalType
    character(len=19) :: massMatr, listLoad
    character(len=24) :: ligrel, noojb, modelLigrel, chdepl, listLoad24
    integer(kind=8) :: nbLoad, iexcit
    character(len=24), pointer :: listLoadName(:) => null()
    character(len=24), parameter :: listCellJv = "&&EXLIMA.LISTCELL"
    integer(kind=8) :: nbCell
    integer(kind=8), pointer :: listCell(:) => null()
    integer(kind=8) :: nbPouxLoad
    aster_logical :: partialMesh
    aster_logical :: lLoadsFromUser, lLoadsHarmo
    character(len=8) :: resultIn
    character(len=16) :: resultType
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initialisations
    nbLoad = 0
    nbPouxLoad = 0
    listLoad = '&&CCVEPO.LISCHA'
    numeStore0 = postCompResu%listStore(1)
    resultIn = postCompResu%resultIn
    resultType = postCompResu%resultType

! - We need a model
    if (model .eq. " ") then
        call utmess('F', 'CALCCHAMP1_10')
    end if
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Get ligrel
    partialMesh = hasCellsDefinedFromCmd(factorKeyword, iOccZero)
    if (partialMesh) then
! ----- Generate name of ligrel
        noojb = '12345678.LIGR000000.LIEL'
        call gnomsd(' ', noojb, 14, 19)
        ligrel = noojb(1:19)

! ----- Get list of cells from user
        call getelem(mesh, factorKeyword, iOccZero, 'F', listCellJv, nbCell, model=model)
        ASSERT(nbCell .ne. 0)
        call jeveuo(listCellJv, 'L', vi=listCell)

! ----- Create FED (ligrel) from list of cells
        call createFEDFromList(model, jvBase, ligrel, &
                               nbCell, listCell)

! ----- Cleaning
        call jedetr(listCellJv)

    else
        ligrel = modelLigrel
    end if

! - Modal type
    modalType = ' '
    if (resultType .eq. 'MODE_MECA') then
        call rsadpa(resultIn, 'L', 1, 'TYPE_MODE', 1, 0, sjv=jvPara)
        modalType = zk16(jvPara)
    end if

! - Detect POUX on FED
    call dismoi('EXI_POUX', ligrel, 'LIGREL', repk=answer)
    postCompPoux%lPoux = answer .eq. "OUI"

    if (postCompPoux%lPoux) then
! ----- Get option used to compute mass matrix
        if ((resultType .eq. 'MODE_MECA' .and. modalType(1:8) .eq. 'MODE_DYN') .or. &
            resultType .eq. 'DYNA_TRANS' .or. &
            resultType .eq. 'MODE_ACOU' .or. &
            resultType .eq. 'DYNA_HARMO') then
            postCompPoux%optionMass = 'MASS_MECA'
            call dismoi('REF_MASS_PREM', resultIn, 'RESU_DYNA', &
                        repk=massMatr, arret='C', ier=ierd)
            if (massMatr .ne. ' ') then
                call dismoi('SUR_OPTION', massMatr, 'MATR_ASSE', &
                            repk=postCompPoux%optionMass, arret='C', ier=ierd)
            end if

! --------- Create M.Gamma vector
            call rsexch('F', resultIn, 'DEPL', 1, chdepl, ierd)
            call copich('V', chdepl, postCompPoux%chdynr)
        end if

! ----- Get list of loads
        if (option .eq. "REAC_NODA") then
            call medome_once(resultIn, postCompResu%listStore, postCompResu%nbStore, &
                             list_load_=listLoad, model_=modelOnce)
        else
            call rs_get_listload(resultIn, numeStore0, listLoad, iexcit)
            if (iexcit .eq. 1) then
                call getListLoadsUser(postCompResu, model, &
                                      lLoadsFromUser, listLoad24, lLoadsHarmo)
                listLoad = listLoad24(1:19)
            end if
        end if
        call getNbLoadsFromList(listLoad, nbLoad)

! ----- Detect distributed load on beam
        if (nbLoad .ne. 0) then
            call jeveuo(listLoad(1:19)//'.LCHA', 'L', vk24=listLoadName)
            call cochre(listLoadName, nbLoad, nbPouxLoad, postCompPoux%loadIndx)
            if (nbPouxLoad .gt. 1) then
                call utmess('F', 'CALCULEL2_92')
            end if
        end if

    end if
!
    call jedema()
!
end subroutine
