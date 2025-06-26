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
subroutine comp_meta_read(metaPrepBehaviour)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lcinfo.h"
#include "asterfort/getvtx.h"
!
    type(META_PrepBehaviour), intent(inout) :: metaPrepBehaviour
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Read informations from command file
!
! --------------------------------------------------------------------------------------------------
!
! IO  metaPrepBehaviour: datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: factorKeyword
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: metaType, metaLaw
    integer(kind=8) :: nbCompElem, numeComp, nbVari, idummy, idummy2, iret, nbPhase
    character(len=16) :: compElem(2), compCodePY, metaCodePY
    aster_logical :: hasTemper
!
! --------------------------------------------------------------------------------------------------
!
    factorKeyword = metaPrepBehaviour%factorKeyword
    nbFactorKeyword = metaPrepBehaviour%nbFactorKeyword
    hasTemper = ASTER_FALSE

! - Read informations in CALC_META
    do iFactorKeyword = 1, nbFactorKeyword
        call getvtx(factorKeyword, 'RELATION', iocc=iFactorKeyword, scal=metaType, nbret=iret)
        call getvtx(factorKeyword, 'LOI_META', iocc=iFactorKeyword, scal=metaLaw, nbret=iret)

! ----- Create composite
        nbCompElem = 2
        compElem(1) = metaType
        compElem(2) = metaLaw
        call lccree(nbCompElem, compElem, compCodePY)
        nbCompElem = 1
        compElem(1) = metaType
        call lccree(nbCompElem, compElem, metaCodePY)
        if (metaType .eq. "ACIER_REVENU") then
            hasTemper = ASTER_TRUE
        end if

! ----- Get number of variables and index of behaviour
        call lcinfo(compCodePY, numeComp, nbVari, idummy)
        call lcinfo(metaCodePY, idummy, nbPhase, idummy2)

! ----- Save values
        metaPrepBehaviour%paraBehaviour(iFactorKeyword)%metaType = metaType
        metaPrepBehaviour%paraBehaviour(iFactorKeyword)%metaLaw = metaLaw
        metaPrepBehaviour%paraBehaviour(iFactorKeyword)%nbVari = nbVari
        metaPrepBehaviour%paraBehaviour(iFactorKeyword)%nbPhase = nbPhase
        metaPrepBehaviour%paraBehaviour(iFactorKeyword)%numeComp = numeComp

! ----- Clean
        call lcdiscard(compCodePY)
        call lcdiscard(metaCodePY)
    end do

    metaPrepBehaviour%hasTemper = hasTemper
!
end subroutine
