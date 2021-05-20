! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine meamme(optionz,&
                  modelz, nbLoad, listLoadK24,&
                  matez, matecoz, caraElemz,&
                  time, basez,&
                  matrRigiz, matrMassz,&
                  matrElemz, listElemCalcz,&
                  variz_, comporz_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnl.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
!
character(len=*), intent(in) :: optionz
character(len=*), intent(in) :: modelz
integer, intent(in) :: nbLoad
character(len=24), pointer :: listLoadK24(:) 
character(len=*), intent(in) :: matez, matecoz, caraElemz
real(kind=8), intent(in) :: time
character(len=*), intent(in) :: basez
character(len=*), intent(in) :: matrRigiz, matrMassz, matrElemz
character(len=*), intent(in) :: listElemCalcz
character(len=*), intent(in), optional :: variz_, comporz_
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix for AMOR_MECA / RIGI_MECA_HYST
!
! NB: careful, compute Dirichlet [B] matrix too when RIGI_MECA_HYST
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  model            : name of model
! In  nbLoad           : number of loads
! In  listLoadK24      : pointer to the name of loads
! In  mate             : name of material characteristics (field)
! In  mateco           : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  time             : current time
! In  base             : JEVEUX base to create matrElem
! In  matrRigi         : elementary rigidity matrix
! In  matrMass         : elementary rigidity mass
! In  listElemCalc     : list of element (LIGREL) where matrElem is computed
! In  matrElem         : elementary matrix
! In  modeFourier      : index of Fourier mode
! In  vari             : internal state variables
! In  compor           : field of behaviour (non-linear cases)
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbFieldInMax = 14, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    character(len=16), parameter :: phenom = 'MECANIQUE'
    integer :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer :: iret
    integer, parameter :: modeFourier = 0
    character(len=16) :: option
    character(len=24), parameter :: chvarc = '&&MEAMME.CHVARC'
    character(len=24) :: compor, listElemCalc
    character(len=8) :: physQuantityName
    character(len=24) :: matrRigi, matrMass
    character(len=24) :: resuElemRigi, resuElemMass
    character(len=24) :: chgeom, chcara(18), chharm
    character(len=1) :: base
    character(len=8) :: model, caraElem
    character(len=24) :: mate, mateco
    character(len=19) :: matrElem, resuElem, ligrel
    integer :: iLoad, indxResuElem
    integer :: nbResuElem, iResuElem, idxResuElemRigi
    integer :: nbSubstruct
    character(len=24), pointer :: rerr(:) => null()
    character(len=24), pointer :: listResuElem(:) => null()
    aster_logical :: hasDirichlet
    character(len=8) :: loadName
    character(len=13) :: loadDescBase
    character(len=19) :: loadMapName, loadLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    option       = optionz
    model        = modelz
    caraElem     = caraElemz
    mate         = matez
    mateco       = matecoz
    matrElem     = matrElemz
    base         = basez
    listElemCalc = listElemCalcz
    matrRigi     = matrRigiz
    matrMass     = matrMassz
    hasDirichlet = option .eq. 'RIGI_MECA_HYST'
    lpain  = ' '
    lchin  = ' '
    lpaout = ' '
    lchout = ' '

! - Prepare flags
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi = nbSubstruct)

! - Behaviour when non-linear case, multi-behaviour (PMF) for linear case
    compor = ' '
    if (present(comporz_)) then
        compor = comporz_
    else
        compor = mate(1:8)//'.COMPOR'
    endif

! - Preparation of input fields
    call mecham(option, model, caraElem, modeFourier, chgeom,&
                chcara, chharm, iret)

! - Field for external state variables
    call vrcins(model, mate, caraElem, time, chvarc, codret)

! - Get RESU_ELEM from rigidity matrix
    resuElemRigi    = ' '
    idxResuElemRigi = 0
    if (matrRigi(1:1) .ne. ' ') then
        call jeexin(matrRigi(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24 = listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemRigi    = listResuElem(iResuElem)
                idxResuElemRigi = iResuElem
                call dismoi('NOM_LIGREL', resuElemRigi, 'RESUELEM', repk=ligrel)
                if (ligrel(1:8) .eq. model(1:8)) then
                    goto 20
                endif
            end do
            ASSERT(ASTER_FALSE)
 20         continue
        endif
    endif

! - Get RESU_ELEM from mass matrix
    resuElemMass = ' '
    if (matrMass(1:1) .ne. ' ') then
        call jeexin(matrMass(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrMass(1:19)//'.RELR', 'L', vk24 = listResuElem)
            call jelira(matrMass(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemMass = listResuElem(iResuElem)
                call dismoi('NOM_LIGREL', resuElemMass, 'RESUELEM', repk=ligrel)
                if (ligrel(1:8) .eq. model(1:8)) then
                    goto 40
                endif
            end do
            ASSERT(ASTER_FALSE)
 40         continue
        endif
    endif

! - Prepare RESU_ELEM objects
    call memare(base, matrElem, model, mate, caraElem, 'AMOR_MECA')
    call jeveuo(matrElem//'.RERR', 'E', vk24 = rerr)
    call jedetr(matrElem//'.RELR')
    if (nbSubstruct .gt. 0) then
        rerr(3) = 'OUI_SOUS_STRUC'
    endif

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = matecoz(1:19)
    lpain(3) = 'PCAORIE'
    lchin(3) = chcara(1)(1:19)
    lpain(4) = 'PCADISA'
    lchin(4) = chcara(4)(1:19)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6)(1:19)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7)(1:19)
    lpain(7) = 'PVARCPR'
    lchin(7) = chvarc(1:19)
    lpain(8) = 'PCADISK'
    lchin(8) = chcara(2)(1:19)
    lpain(9) = 'PCINFDI'
    lchin(9) = chcara(15)(1:19)
    lpain(10) = 'PMASSEL'
    lchin(10) = resuElemMass(1:19)
    lpain(11) = 'PCOMPOR'
    lchin(11) = compor(1:19)
    nbFieldIn = 11

! - Add internal state variables
    if (present(variz_)) then
        nbFieldIn = nbFieldIn + 1
        lpain(nbFieldIn) = 'PVARIPG'
        lchin(nbFieldIn) = variz_(1:19)
    endif

! - Get symmetric or unsymmetric rigidity matrix
    if (resuElemRigi .ne. ' ') then
        nbFieldIn = nbFieldIn + 1
        lchin(nbFieldIn) = resuElemRigi(1:19)
        call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
        if (physQuantityName .eq. 'MDNS_R') then
            lpain(nbFieldIn) = 'PRIGINS'
        else
            lpain(nbFieldIn) = 'PRIGIEL'
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24 = listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            if (idxResuElemRigi .lt. nbResuElem) then
                resuElemRigi = listResuElem(idxResuElemRigi+1)
                call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
                if (physQuantityName .eq. 'MDNS_R') then
                    nbFieldIn = nbFieldIn + 1
                    lpain(nbFieldIn) = 'PRIGINS'
                    lchin(nbFieldIn) = resuElemRigi(1:19)
                endif
            endif
        endif
    endif

! - Output fields
    if (option .eq. 'AMOR_MECA') then
        lpaout(1) = 'PMATUUR'
        lpaout(2) = 'PMATUNS'
    else if (option .eq. 'RIGI_MECA_HYST') then
        lpaout(1) = 'PMATUUC'
    else
        ASSERT(ASTER_FALSE)
    endif
    lchout(1) = matrElem(1:8)//'.ME001'
    lchout(2) = matrElem(1:8)//'.ME002'
    nbFieldOut = 2

! - Compute
    call calcul('S',&
                option, listElemCalc,&
                nbFieldIn, lchin, lpain,&
                nbFieldOut, lchout, lpaout,&
                base, 'OUI')

! - Save RESU_ELEM
    call reajre(matrElem, lchout(1), base)
    call reajre(matrElem, lchout(2), base)

! - Dirichlet
    option       = 'MECA_DDLM_R'
    nbFieldIn    = 1
    nbFieldOut   = 1
    resuElem     = matrElem(1:8)//'.XXXXXXX'
    call jelira(matrElem(1:19)//'.RELR', 'LONUTI', indxResuElem)
    indxResuElem = indxResuElem + 1
    if (hasDirichlet) then
        do iLoad = 1, nbLoad
! --------- Current load
            loadName    = listLoadK24(iLoad)(1:8)
            call lisnnl(phenom, loadName, loadDescBase)
            loadMapName = loadDescBase//'.CMULT'
            loadLigrel  = loadDescBase//'.LIGRE'

! --------- Detect if current load is OK
            call jeexin(loadLigrel(1:19)//'.LIEL', iret)
            if (iret .le. 0) cycle
            call exisd('CHAMP_GD', loadMapName, iret)
            if (iret .le. 0) cycle

! --------- Input field
            lpain(1) = 'PDDLMUR'
            lchin(1) = loadMapName

! --------- Generate new RESU_ELEM name
            call codent(indxResuElem, 'D0', resuElem(10:16))

! --------- Output field
            lpaout(1) = 'PMATUUR'
            lchout(1) = resuElem

! --------- Compute
            call calcul('S',&
                        option, loadLigrel,&
                        nbFieldIn, lchin, lpain,&
                        nbFieldOut, lchout, lpaout,&
                        base, 'OUI')

! --------- Save RESU_ELEM
            call reajre(matrElem, resuElem, base)
            indxResuElem = indxResuElem + 1
            if (indxResuElem .eq. 9999999) then
                call utmess('F', 'CHARGES6_82', sk = 'RIGI_MECA_HYST')
            endif

        end do
    endif

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
