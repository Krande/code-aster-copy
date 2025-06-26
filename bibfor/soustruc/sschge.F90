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
subroutine sschge(macrElemZ)
!
    use listLoad_module
    use listLoad_type
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/me2mme.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/ss2mm2.h"
#include "asterfort/ssvau1.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: macrElemZ
!
! --------------------------------------------------------------------------------------------------
!     BUT: TRAITER LE MOT CLEF "CAS_CHARGE"
!             DE L'OPERATEUR MACR_ELEM_STAT
!     LES CHARGEMENTS SERONT CONDENSES LORS DE L'ASSEMBLAGE.
!
!     IN: NOMACR : NOM DU MACR_ELEM_STAT
!
!     OUT: LES OBJETS SUIVANTS DU MACR_ELEM_STAT SONT CALCULES:
!             NOMACR.LICA(NOMCAS)
!             NOMACR.LICH(NOMCAS)
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: phenom = "MECA"
    character(len=1), parameter :: jvBase = 'V'
    character(len=8) :: answer, macrElem
    character(len=14) :: numeDof
    character(len=8) :: model, materField, caraElem
    character(len=24) :: mateco
    character(len=8) :: caseName
    character(len=19) :: vectAsse, vectElem
    character(len=24), parameter :: listLoad = "&&SSCHGE_LISTLOAD"
    character(len=16), parameter :: factKeyword = 'CAS_CHARGE'
    integer(kind=8) :: jvLica
    integer(kind=8) :: n1, n2
    integer(kind=8) :: nbDofExte, nbDofInte, nbDof
    integer(kind=8) :: nbLoad, nbLoadMax
    integer(kind=8) :: nbCase, caseNume, iCase
    real(kind=8) :: time
    integer(kind=8), parameter :: numeHarm = 0
    character(len=8), pointer :: refm(:) => null(), macrElemLich(:) => null()
    character(len=8), pointer :: loadName(:) => null()
    real(kind=8), pointer :: macrElemLica(:) => null()
    integer(kind=8), pointer :: desm(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    type(ListLoad_Prep) :: listLoadPrep
    aster_logical, parameter :: kineExcl = ASTER_FALSE, diriExcl = ASTER_FALSE
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    macrElem = macrElemZ

! - Get parameters from super-element
    call jeveuo(macrElem//'.REFM', 'L', vk8=refm)
    model = refm(1)
    materField = refm(3)
    if (materField .ne. '        ') then
        call rcmfmc(materField, mateco, l_ther_=ASTER_FALSE)
    else
        mateco = ' '
    end if
    caraElem = refm(4)
    numeDof = refm(5)
    ASSERT(numeDof(1:8) .eq. macrElem)

! - Get size of super-element
    call jeveuo(macrElem//'.DESM', 'E', vi=desm)
    nbDofExte = desm(4)
    nbDofInte = desm(5)
    nbDof = nbDofExte+nbDofInte

    call jelira(macrElem//'.LICH', 'LONMAX', nbLoadMax)
    nbLoadMax = nbLoadMax-1
    vectElem = '&&VECEL            '
    vectAsse = macrElem//'.CHARMECA'

    listLoadPrep%model = model

    call getfac(factKeyword, nbCase)
    do iCase = 1, nbCase
! ----- Get current case name
        call getvtx(factKeyword, 'NOM_CAS', iocc=iCase, scal=caseName, nbret=n1)

! ----- Creation objets de la collection
        call jecroc(jexnom(macrElem//'.LICA', caseName))
        call jenonu(jexnom(macrElem//'.LICA', caseName), caseNume)
        call jecroc(jexnom(macrElem//'.LICH', caseName))
        call jeveuo(jexnum(macrElem//'.LICH', caseNume), 'E', vk8=macrElemLich)

! ----- Update list of loads
        call getvtx(factKeyword, 'SUIV', iocc=iCase, scal=answer, nbret=n1)
        if (answer .eq. 'OUI') then
            macrElemLich(1) = 'OUI_SUIV'
        else
            macrElemLich(1) = 'NON_SUIV'
        end if
        call getvid(factKeyword, 'CHARGE', iocc=iCase, nbval=0, nbret=n1)
        nbLoad = abs(n1)
        if (nbLoad .gt. nbLoadMax) then
            call utmess('F', 'SOUSTRUC_40')
        end if
        AS_ALLOCATE(vk8=loadName, size=nbLoad)
        call getvid(factKeyword, 'CHARGE', iocc=iCase, nbval=nbLoad, vect=loadName, nbret=n2)
        macrElemLich(2:nbLoad+1) = loadName(1:nbLoad)

! ----- Create list of loads
        call creaListLoadFromList(phenom, listLoadPrep, &
                                  listLoad, jvBase, &
                                  nbLoad, loadName, &
                                  kineExcl, diriExcl)

! ----- Get time
        call getvr8(factKeyword, 'INST', iocc=iCase, scal=time, nbret=n2)

! ----- Compute elementary vectors for loads
        call me2mme(model, listLoad, &
                    materField, mateco, caraElem, &
                    time, vectElem, numeHarm, jvBase)
        call ss2mm2(model, vectElem, caseName)

! ----- Assemble elementary vectors for loads
        call assvec('V', vectAsse, 1, vectElem, [1.d0], numeDof)

! ----- Copy value of load for case
        call jeveuo(vectAsse//'.VALE', 'L', vr=vale)
        call jeveuo(jexnum(macrElem//'.LICA', caseNume), 'E', vr=macrElemLica)
        call jeveuo(jexnum(macrElem//'.LICA', caseNume), 'E', jvLica)
        macrElemLica(1:nbDof) = vale(1:nbDof)

!       -- CONDENSATION DE .LICA(1:NDDLT) DANS .LICA(NDDLT+1,2*NDDLT) :
        call ssvau1(macrElem, jvLica, jvLica+nbDof)

!       -- ON COMPTE LES CAS DE CHARGE EFFECTIVEMENT CALCULES:
        desm(7) = caseNume
!
        call detrsd('VECT_ELEM', vectElem)
        call detrsd('CHAMP_GD', vectAsse)
        call detrsd('LISTE_CHARGES', listLoad)
        AS_DEALLOCATE(vk8=loadName)
    end do
!
    call jedema()
end subroutine
