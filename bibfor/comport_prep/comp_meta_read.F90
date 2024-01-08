! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine comp_meta_read(metaPrepPara)
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
    type(META_PrepPara), intent(inout) :: metaPrepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Read informations from command file
!
! --------------------------------------------------------------------------------------------------
!
! IO  metaPrepPara     : datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer :: i_comp, nb_comp
    character(len=16) :: metaType, metaLaw
    integer :: nbCompElem, numeComp, nbVari, idummy, idummy2, iret, nbPhase
    character(len=16) :: compElem(2), compCodePY, metaCodePY
!
! --------------------------------------------------------------------------------------------------
!
    nb_comp = metaPrepPara%nb_comp

! - Read informations in CALC_META
    do i_comp = 1, nb_comp
        call getvtx(factorKeyword, 'RELATION', iocc=i_comp, scal=metaType, nbret=iret)
        call getvtx(factorKeyword, 'LOI_META', iocc=i_comp, scal=metaLaw, nbret=iret)

! ----- Create composite
        nbCompElem = 2
        compElem(1) = metaType
        compElem(2) = metaLaw
        call lccree(nbCompElem, compElem, compCodePY)
        nbCompElem = 1
        compElem(1) = metaType
        call lccree(nbCompElem, compElem, metaCodePY)

! ----- Get number of variables and index of behaviour
        call lcinfo(compCodePY, numeComp, nbVari, idummy)
        call lcinfo(metaCodePY, idummy, nbPhase, idummy2)

! ----- Save values
        metaPrepPara%para(i_comp)%metaType = metaType
        metaPrepPara%para(i_comp)%metaLaw = metaLaw
        metaPrepPara%para(i_comp)%nbVari = nbVari
        metaPrepPara%para(i_comp)%nbPhase = nbPhase
        metaPrepPara%para(i_comp)%numeComp = numeComp

! ----- Clean
        call lcdiscard(compCodePY)
        call lcdiscard(metaCodePY)
    end do
!
end subroutine
