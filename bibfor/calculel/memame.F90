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
subroutine memame(optionz, modelz, matez, matecoz, caraElemz, time, &
                  comporMultz, matrElemz, basez, listElemCalcz)
!
    use HHO_precalc_module, only: hhoAddInputField
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/vrcins.h"
#include "asterfort/xajcin.h"
!
    character(len=*), intent(in) :: optionz
    character(len=*), intent(in) :: modelz, matez, matecoz, caraElemz
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: comporMultz, matrElemz
    character(len=*), intent(in) :: basez, listElemCalcz
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix MASS_*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option
! In  model            : name of the model
! In  mate             : name of material characteristics (field)
! In  mateco           : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  time             : current time
! In  comporMult       : name of comportment definition for PMF (field)
! In  base             : JEVEUX base to create matrElem
! In  matrElem         : elementary matrix
! In  listElemCalc     : list of elements (LIGREL) where matrElem is computed
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 33, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer(kind=8) :: iret
    integer(kind=8), parameter :: modeFourier = 0
    character(len=16) :: option
    character(len=24), parameter :: chvarc = '&&MERIME.CHVARC'
    character(len=24) :: comporMult, listElemCalc
    character(len=24) :: chgeom, chcara(18), chharm
    character(len=1) :: base
    character(len=8) :: model, caraElem
    character(len=24) :: mate, mateco
    character(len=19) :: matrElem
    integer(kind=8) :: nbSubstruct
    aster_logical :: lxfem, hasFiniteElement
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    option = optionz
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
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSubstruct)

! - Preparation of input fields
    call mecham(option, model, caraElem, modeFourier, chgeom, &
                chcara, chharm, iret)
    hasFiniteElement = iret .eq. 0

! - Field for external state variables
    call vrcins(model, mate, caraElem, time, chvarc, codret)

! - Prepare RESU_ELEM objects
    call jeexin(matrElem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call memare(base, matrElem, model, option, to_aster_logical(nbSubstruct > 0))
    else
        call jedetr(matrElem(1:19)//'.RELR')
    end if

! - Input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = mateco(1:19)
    lpain(3) = 'PCAORIE'
    lchin(3) = chcara(1) (1:19)
    lpain(4) = 'PCADISM'
    lchin(4) = chcara(3) (1:19)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6) (1:19)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7) (1:19)
    lpain(7) = 'PCASECT'
    lchin(7) = chcara(8) (1:19)
    lpain(8) = 'PVARCPR'
    lchin(8) = chvarc(1:19)
    lpain(9) = 'PCAARPO'
    lchin(9) = chcara(9) (1:19)
    lpain(10) = 'PCACABL'
    lchin(10) = chcara(10) (1:19)
    lpain(11) = 'PCAGEPO'
    lchin(11) = chcara(5) (1:19)
    lpain(12) = 'PABSCUR'
    lchin(12) = chgeom(1:8)//'.ABSC_CURV'
    lpain(13) = 'PCAGNBA'
    lchin(13) = chcara(11) (1:19)
    lpain(14) = 'PCAPOUF'
    lchin(14) = chcara(13) (1:19)
    lpain(15) = 'PCOMPOR'
    lchin(15) = comporMult(1:19)
    lpain(16) = 'PNBSP_I'
    lchin(16) = chcara(16) (1:19)
    lpain(17) = 'PFIBRES'
    lchin(17) = chcara(17) (1:19)
    lpain(18) = 'PCINFDI'
    lchin(18) = chcara(15) (1:19)
    nbFieldIn = 18

! - Add input XFEM fields if required
    if (lxfem) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if
!
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)
!
! - Output fields
    lpaout(1) = 'PMATUUR'
    lchout(1) = matrElem(1:15)//'.M01'
    lpaout(2) = 'PMATUNS'
    lchout(2) = matrElem(1:15)//'.M02'
    if (option .eq. 'MASS_MECA') then
        nbFieldOut = 2
    else
        nbFieldOut = 1
    end if

! - Mass
    if (hasFiniteElement) then
! ----- Compute
        call calcul('S', &
                    option, listElemCalc, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    base, 'OUI')

! ----- Save RESU_ELEM
        call reajre(matrElem, lchout(1), base)
        if (nbFieldOut .eq. 2) then
            call reajre(matrElem, lchout(2), base)
        end if

    end if

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
