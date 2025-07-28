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
subroutine merime(modelz, nbLoad, listLoadK24, &
                  matez, matecoz, caraElemz, &
                  time, comporMultz, matrElemz, modeFourier, &
                  basez, listElemCalcz, hasExteStatVari_, onlyDirichlet_)
!
!                          DE MERIME...  PROSPER YOUP-LA-BOUM!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/lisnnl.h"
#include "asterfort/mecact.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: modelz
    integer(kind=8), intent(in) :: nbLoad
    character(len=24), pointer :: listLoadK24(:)
    character(len=*), intent(in) :: matez, matecoz, caraElemz
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: comporMultz, matrElemz
    integer(kind=8), intent(in) :: modeFourier
    character(len=*), intent(in) :: basez, listElemCalcz
    aster_logical, intent(in), optional :: hasExteStatVari_, onlyDirichlet_
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix for RIGI_MECA
!
! NB: careful, compute Dirichlet [B] matrix too
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  nbLoad           : number of loads
! In  listLoadK24      : pointer to the name of loads
! In  mate             : name of material characteristics (field)
! In  mateco           : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  time             : current time
! In  comporMult       : name of comportment definition for PMF (field)
! In  matrElem         : elementary matrix
! In  modeFourier      : index of Fourier mode
! In  base             : JEVEUX base to create matrElem
! In  listElemCalc     : list of elements (LIGREL) where matrElem is computed
! In  onlyDirichlet    : flag to compute only Dirichlet [B] matrix
! In  hasExteStatVari  : flag to use external state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 32, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    character(len=16), parameter :: phenom = 'MECANIQUE'
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer(kind=8) :: iret
    character(len=16) :: option
    character(len=24), parameter :: chvarc = '&&MERIME.CHVARC'
    character(len=24), parameter :: chtime = '&&MERIME.CHTIME'
    character(len=24) :: comporMult, listElemCalc
    character(len=24) :: chgeom, chcara(18), chharm
    character(len=1) :: base
    character(len=8) :: model, caraElem
    character(len=24) :: mate, mateco
    character(len=19) :: matrElem, resuElem
    integer(kind=8) :: iLoad, indxResuElem
    integer(kind=8) :: nbSubstruct
    aster_logical :: lxfem, hasFiniteElement, hasExteStatVari, onlyDirichlet
    character(len=8) :: loadName
    character(len=13) :: loadDescBase
    character(len=19) :: loadMapName, loadLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    option = 'RIGI_MECA'
    model = modelz
    caraElem = caraElemz
    mate = matez
    mateco = matecoz
    matrElem = matrElemz
    comporMult = comporMultz
    base = basez
    listElemCalc = listElemCalcz
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Prepare flags
    call exixfe(model, iret)
    lxfem = iret .ne. 0
    onlyDirichlet = ASTER_FALSE
    if (present(onlyDirichlet_)) then
        onlyDirichlet = onlyDirichlet_
    end if
    hasExteStatVari = ASTER_FALSE
    if (present(hasExteStatVari_)) then
        hasExteStatVari = hasExteStatVari_
    end if
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSubstruct)

! - Preparation of input fields
    call mecham(option, model, caraElem, modeFourier, chgeom, &
                chcara, chharm, iret)
    hasFiniteElement = iret .eq. 0

! - Field for time
    call mecact('V', chtime, 'MODELE', model, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=time)

! - Field for external state variables
    if (hasExteStatVari) then
        call vrcins(model, mate, caraElem, time, chvarc, codret)
    end if

! - Prepare RESU_ELEM objects
    call memare(base, matrElem, model, option, to_aster_logical(nbSubstruct > 0))
    call jedetr(matrElem//'.RELR')

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = mateco(1:19)
    lpain(3) = 'PCAORIE'
    lchin(3) = chcara(1) (1:19)
    lpain(4) = 'PCADISK'
    lchin(4) = chcara(2) (1:19)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6) (1:19)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7) (1:19)
    lpain(7) = 'PCASECT'
    lchin(7) = chcara(8) (1:19)
    lpain(8) = 'PCAARPO'
    lchin(8) = chcara(9) (1:19)
    lpain(9) = 'PNBSP_I'
    lchin(9) = chcara(16) (1:19)
    lpain(10) = 'PCOMPOR'
    lchin(10) = comporMult(1:19)
    lpain(11) = 'PHARMON'
    lchin(11) = chharm(1:19)
    lpain(12) = 'PVARCPR'
    lchin(12) = chvarc(1:19)
    lpain(13) = 'PFIBRES'
    lchin(13) = chcara(17) (1:19)
    lpain(14) = 'PCAGNBA'
    lchin(14) = chcara(11) (1:19)
    lpain(15) = 'PCAMASS'
    lchin(15) = chcara(12) (1:19)
    lpain(16) = 'PCAPOUF'
    lchin(16) = chcara(13) (1:19)
    lpain(17) = 'PCAGEPO'
    lchin(17) = chcara(5) (1:19)
    lpain(18) = 'PINSTR'
    lchin(18) = chtime(1:19)
    lpain(19) = 'PCINFDI'
    lchin(19) = chcara(15) (1:19)
    nbFieldIn = 19

! - Add input XFEM fields if required
    if (lxfem) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Output fields
    lpaout(1) = 'PMATUUR'
    lchout(1) = matrElem(1:8)//'.ME001'
    lpaout(2) = 'PMATUNS'
    lchout(2) = matrElem(1:8)//'.ME002'
    nbFieldOut = 2

! - Rigidity
    if (.not. onlyDirichlet) then
        if ((lxfem) .or. ((.not. lxfem) .and. hasFiniteElement)) then
! --------- Compute
            call calcul('S', &
                        option, listElemCalc, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        base, 'OUI')

! --------- Save RESU_ELEM
            call reajre(matrElem, lchout(1), base)
            call reajre(matrElem, lchout(2), base)

        end if
    end if

! - Dirichlet
    option = 'MECA_DDLM_R'
    nbFieldIn = 1
    nbFieldOut = 1
    resuElem = matrElem(1:8)//'.XXXXXXX'
    indxResuElem = 0
    do iLoad = 1, nbLoad
! ----- Current load
        loadName = listLoadK24(iLoad) (1:8)
        call lisnnl(phenom, loadName, loadDescBase)
        loadMapName = loadDescBase//'.CMULT'
        loadLigrel = loadDescBase//'.LIGRE'

! ----- Detect if current load is OK
        call jeexin(loadLigrel(1:19)//'.LIEL', iret)
        if (iret .le. 0) cycle
        call exisd('CHAMP_GD', loadMapName, iret)
        if (iret .le. 0) cycle

! ----- Input field
        lpain(1) = 'PDDLMUR'
        lchin(1) = loadMapName

! ----- Generate new RESU_ELEM name
        call codent(indxResuElem, 'D0', resuElem(10:16))

! ----- Output field
        lpaout(1) = 'PMATUUR'
        lchout(1) = resuElem

! ----- Compute
        call calcul('S', &
                    option, loadLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    base, 'OUI')

! ----- Save RESU_ELEM
        call reajre(matrElem, resuElem, base)
        indxResuElem = indxResuElem+1
        if (indxResuElem .eq. 9999999) then
            call utmess('F', 'CHARGES6_82', sk='RIGI_MECA')
        end if

    end do

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chtime)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
