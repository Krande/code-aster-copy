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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romAlgoNLRead(paraAlgo)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/getvtx.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/romBaseGetInfo.h"
#include "asterfort/romAlgoNLPrintInfo.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_AlgoPara), intent(inout) :: paraAlgo
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Solving non-linear problem
!
! Read parameters for algorithm management
!
! --------------------------------------------------------------------------------------------------
!
! IO  paraAlgo       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: coef_pena
    character(len=8) :: baseName
    character(len=16) :: keywf, answer
    character(len=24) :: grnode_int, grnode_sub
    aster_logical :: l_hrom, l_hrom_corref
    type(ROM_DS_Empi) :: base
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM5_41')
    end if
!
! - Initializations
!
    keywf = 'MODELE_REDUIT'
    l_hrom = ASTER_FALSE
    l_hrom_corref = ASTER_FALSE
    grnode_int = ' '
    grnode_sub = ' '
    coef_pena = 0.d0
!
! - Read parameters
!
    call getvid(keywf, 'BASE_PRIMAL', iocc=1, scal=baseName)
    call getvtx(keywf, 'DOMAINE_REDUIT', iocc=1, scal=answer)
    l_hrom = answer .eq. 'OUI'
    if (l_hrom) then
        call getvtx(keywf, 'GROUP_NO_INTERF', iocc=1, scal=grnode_int)
        call getvtx(keywf, 'CORR_COMPLET', iocc=1, scal=answer)
        l_hrom_corref = answer .eq. 'OUI'
        if (l_hrom_corref) then
            call getvtx(keywf, 'GROUP_NO_ENCASTRE', iocc=1, scal=grnode_sub)
            call getvr8(keywf, 'COEF_PENA', iocc=1, scal=coef_pena)
        end if
    end if
!
! - Get informations about base
!
    call romBaseGetInfo(baseName, base)
!
! - Save parameters in datastructure
!
    paraAlgo%l_rom = ASTER_TRUE
    paraAlgo%ds_empi = base
    paraAlgo%l_hrom = l_hrom
    paraAlgo%grnode_int = grnode_int
    paraAlgo%l_hrom_corref = l_hrom_corref
    paraAlgo%grnode_sub = grnode_sub
    paraAlgo%vale_pena = coef_pena
!
! - Debug
!
    if (niv .ge. 2) then
        call romAlgoNLPrintInfo(paraAlgo)
    end if
!
end subroutine
