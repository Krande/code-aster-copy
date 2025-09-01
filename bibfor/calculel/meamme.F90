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
subroutine meamme(modelz, &
                  matez, matecoz, caraElemz, &
                  time, basez, &
                  matrRigiz, matrMassz, &
                  matrElemz, &
                  variz, comporz, sddyna)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/ndynkk.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/vrcins.h"
!
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: matez, matecoz, caraElemz
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: basez
    character(len=*), intent(in) :: matrRigiz, matrMassz, matrElemz
    character(len=*), intent(in) :: variz, comporz
    character(len=19), intent(in) :: sddyna
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
    integer(kind=8), parameter :: nbFieldInMax = 16, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    character(len=16), parameter :: option = 'AMOR_MECA'
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer(kind=8) :: iret
    integer(kind=8), parameter :: modeFourier = 0
    character(len=24), parameter :: chvarc = '&&MEAMME.CHVARC'
    character(len=24) :: compor, vari
    character(len=8) :: physQuantityName
    character(len=24) :: matrRigi, matrMass
    character(len=24) :: resuElemRigi, resuElemMass
    character(len=24) :: chgeom, chcara(18), chharm
    character(len=1) :: base
    character(len=8) :: model, caraElem, mesh
    character(len=24) :: mate, mateco, amor_flui
    character(len=19) :: matrElem
    integer(kind=8) :: nbResuElem, iResuElem, idxResuElemRigi
    integer(kind=8) :: nbSubstruct
    character(len=24), pointer :: listResuElem(:) => null()
    character(len=19) :: modelLigrel, resuLigrel
    character(len=24), parameter :: nonLinearMap = "&&MEAMMA.NONLIN"
    integer(kind=8), parameter :: nbCmp = 1
    character(len=8), parameter :: cmpName = ('X1')
    integer(kind=8), parameter :: cmpVale = 1
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelz
    caraElem = caraElemz
    mate = matez
    mateco = matecoz
    matrElem = matrElemz
    base = basez
    matrRigi = matrRigiz
    matrMass = matrMassz
    compor = comporz
    vari = variz
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Get parameters
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSubstruct)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Preparation of input fields
    call mecham(option, model, caraElem, modeFourier, chgeom, &
                chcara, chharm, iret)

! - Special map for non-linear cases
    call jedetr(nonLinearMap)
    call mecact('V', nonLinearMap, 'MAILLA', mesh, 'NEUT_I', &
                ncmp=nbCmp, nomcmp=cmpName, si=cmpVale)

! - Special map for fluid damping
    call ndynkk(sddyna, 'AMOR_FLUI', amor_flui)

! - Field for external state variables
    call vrcins(model, mate, caraElem, time, chvarc, codret)

! - Get RESU_ELEM from rigidity matrix
    resuElemRigi = ' '
    idxResuElemRigi = 0
    if (matrRigi(1:1) .ne. ' ') then
        call jeexin(matrRigi(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemRigi = listResuElem(iResuElem)
                idxResuElemRigi = iResuElem
                call dismoi('NOM_LIGREL', resuElemRigi, 'RESUELEM', repk=resuLigrel)
                if (resuLigrel .eq. modelLigrel) then
                    goto 20
                end if
            end do
            ASSERT(ASTER_FALSE)
20          continue
        end if
    end if

! - Get RESU_ELEM from mass matrix
    resuElemMass = ' '
    if (matrMass(1:1) .ne. ' ') then
        call jeexin(matrMass(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrMass(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrMass(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemMass = listResuElem(iResuElem)
                call dismoi('NOM_LIGREL', resuElemMass, 'RESUELEM', repk=resuLigrel)
                if (resuLigrel .eq. modelLigrel) then
                    goto 40
                end if
            end do
            ASSERT(ASTER_FALSE)
40          continue
        end if
    end if

! - Prepare RESU_ELEM objects
    call memare(base, matrElem, model, 'AMOR_MECA', to_aster_logical(nbSubstruct > 0))
    call jedetr(matrElem//'.RELR')

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = matecoz(1:19)
    lpain(3) = 'PCAORIE'
    lchin(3) = chcara(1) (1:19)
    lpain(4) = 'PCADISA'
    lchin(4) = chcara(4) (1:19)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6) (1:19)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7) (1:19)
    lpain(7) = 'PVARCPR'
    lchin(7) = chvarc(1:19)
    lpain(8) = 'PCADISK'
    lchin(8) = chcara(2) (1:19)
    lpain(9) = 'PCINFDI'
    lchin(9) = chcara(15) (1:19)
    lpain(10) = 'PMASSEL'
    lchin(10) = resuElemMass(1:19)
    lpain(11) = 'PCOMPOR'
    lchin(11) = compor(1:19)
    lpain(12) = 'PNONLIN'
    lchin(12) = nonLinearMap(1:19)
    lpain(13) = 'PVARIPG'
    lchin(13) = vari(1:19)
    lpain(14) = 'PAMORFL'
    lchin(14) = amor_flui(1:19)
    nbFieldIn = 14

! - Get symmetric or unsymmetric rigidity matrix
    if (resuElemRigi .ne. ' ') then
        nbFieldIn = nbFieldIn+1
        lchin(nbFieldIn) = resuElemRigi(1:19)
        call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
        if (physQuantityName .eq. 'MDNS_R') then
            lpain(nbFieldIn) = 'PRIGINS'
        else
            lpain(nbFieldIn) = 'PRIGIEL'
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            if (idxResuElemRigi .lt. nbResuElem) then
                resuElemRigi = listResuElem(idxResuElemRigi+1)
                call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
                if (physQuantityName .eq. 'MDNS_R') then
                    nbFieldIn = nbFieldIn+1
                    lpain(nbFieldIn) = 'PRIGINS'
                    lchin(nbFieldIn) = resuElemRigi(1:19)
                end if
            end if
        end if
    end if

! - Output fields
    lpaout(1) = 'PMATUUR'
    lpaout(2) = 'PMATUNS'
    lchout(1) = matrElem(1:8)//'.ME001'
    lchout(2) = matrElem(1:8)//'.ME002'
    nbFieldOut = 2

! - Compute
    ASSERT(nbFieldIn .le. nbFieldInMax)
    ASSERT(nbFieldOut .le. nbFieldOutMax)
    call calcul('S', &
                option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                base, 'OUI')

! - Save RESU_ELEM
    call reajre(matrElem, lchout(1), base)
    call reajre(matrElem, lchout(2), base)

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
