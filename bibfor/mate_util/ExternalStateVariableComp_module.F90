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
! ==================================================================================================
!
! Module to compute External State Variables
!
! ==================================================================================================
!
module ExternalStateVariableComp_module
! ==================================================================================================
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: varcCompElem
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/gcnco2.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/reajre.h"
#include "asterfort/detrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vemare.h"
#include "asterfort/varcDetect.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! varcCompElem
!
! Compute elementary vector for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine varcCompElem(line, &
                            numeHarm, modelZ, caraElemZ, materFieldZ, matecoZ, &
                            chtimeZ, varcRefeZ, varcZ, &
                            jvBase, vectElemZ, &
                            lCumul_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: line
        integer(kind=8), intent(in) :: numeHarm
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, matecoZ
        character(len=*), intent(in) :: chtimeZ, varcRefeZ, varcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: vectElemZ
        aster_logical, optional, intent(in) :: lCumul_
! ----- Locals
        integer(kind=8), parameter :: nbFieldIn = 18, nbFieldOut = 1
        character(len=8) :: lpain(nbFieldIn), lpaout(nbFieldOut)
        character(len=24) :: lchin(nbFieldIn), lchout(nbFieldOut)
        character(len=8) :: newnom
        character(len=16) :: option
        character(len=24) :: modelLigrel, resuElem
        character(len=24) :: chgeom, chcara(18), chharm
        aster_logical :: lTemp, lHydr, lPtot, lSech, lEpsa, lMeta, lCumul
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(line)

! ----- Initializations
        modelLigrel = modelZ(1:8)//'.MODELE'
        lpain = " "
        lpaout = " "
        lchin = " "
        lchout = " "
        lCumul = ASTER_FALSE
        if (present(lCumul_)) then
            lCumul = lCumul_
        end if

! ----- Linear: only for temperature for the moment
        call varcDetect(materFieldZ, lTemp, lHydr, lPtot, lSech, lEpsa, lMeta)
        if (line) then
            if (lHydr .or. lPtot .or. lSech .or. lEpsa .or. lMeta) then
                call utmess('F', 'SOUSTRUC_18')
            end if
        end if

! ----- Geometry field
        call megeom(modelZ, chgeom)

! ----- Elementary characteristics
        call mecara(caraElemZ, chcara)

! ----- Fourier
        call meharm(modelZ, numeHarm, chharm)

! ----- Input fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = matecoZ
        lpain(3) = 'PCACOQU'
        lchin(3) = chcara(7) (1:19)
        lpain(4) = 'PCAGNPO'
        lchin(4) = chcara(6) (1:19)
        lpain(5) = 'PCADISM'
        lchin(5) = chcara(3) (1:19)
        lpain(6) = 'PCAORIE'
        lchin(6) = chcara(1) (1:19)
        lpain(7) = 'PCAGNBA'
        lchin(7) = chcara(11) (1:19)
        lpain(8) = 'PCAARPO'
        lchin(8) = chcara(9) (1:19)
        lpain(9) = 'PCAMASS'
        lchin(9) = chcara(12) (1:19)
        lpain(10) = 'PCAGEPO'
        lchin(10) = chcara(5) (1:19)
        lpain(11) = 'PNBSP_I'
        lchin(11) = chcara(1) (1:8)//'.CANBSP'
        lpain(12) = 'PFIBRES'
        lchin(12) = chcara(1) (1:8)//'.CAFIBR'
        lpain(13) = 'PHARMON'
        lchin(13) = chharm
        lpain(14) = 'PCINFDI'
        lchin(14) = chcara(15)
        lpain(15) = 'PCADISK'
        lchin(15) = chcara(2)
        lpain(16) = 'PINSTR'
        lchin(16) = chtimeZ
        lpain(17) = 'PVARCRR'
        lchin(17) = varcRefeZ
        lpain(18) = 'PVARCPR'
        lchin(18) = varcZ

! ----- Output field
        lpaout(1) = 'PVECTUR'

! ----- Allocate result
        if (.not. lCumul) then
            call detrsd('VECT_ELEM', vectElemZ)
            call vemare(jvBase, vectElemZ, modelZ)
            call reajre(vectElemZ, ' ', jvBase)
        end if
        newnom = '.0000000'
        resuElem = vectElemZ(1:8)//'.0000000'

! ----- Compute
        if (lTemp) then
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=-1)
            lchout(1) = resuElem
            option = 'CHAR_MECA_TEMP_R'
            call calcul('C', option, modelLigrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, &
                        'OUI')
            call reajre(vectElemZ, resuElem, jvBase)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariableComp_module
