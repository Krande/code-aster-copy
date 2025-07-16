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
subroutine resuGetLoads(model, resultType, listLoadResu)
!
    use listLoad_type
    use listLoad_module
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/utmess.h"
#include "asterfort/nmdoch.h"
#include "asterfort/ntdoch.h"
#include "asterfort/copisd.h"
!
    character(len=8), intent(in) :: model
    character(len=16), intent(in) :: resultType
    character(len=24), intent(out) :: listLoadResu
!
! --------------------------------------------------------------------------------------------------
!
! LIRE_RESU and CREA_RESU
!
! Get loads
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  resultType       : type of results datastructure (EVOL_NOLI, EVOL_THER, )
! Out listLoadResu     : name of datastructure for loads
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: jvBase = "V"
    integer(kind=8) :: nbOcc
    character(len=24), parameter :: listLoad = '&&LRCOMM.LISTLOAD'
    type(ListLoad_Prep) :: listLoadPrep
!
! --------------------------------------------------------------------------------------------------
!
    listLoadResu = ' '
    call getfac('EXCIT', nbOcc)
    if (nbOcc .gt. 0) then

! ----- Read from command file
        listLoadPrep%model = model
        if (resultType .eq. 'EVOL_ELAS' .or. resultType .eq. 'EVOL_NOLI') then
            call nmdoch(listLoadPrep, listLoad, jvBase)
        else if (resultType .eq. 'EVOL_THER') then
            call ntdoch(listLoadPrep, listLoad, jvBase)
        else
            call utmess('A', 'RESULT2_16', sk=resultType)
        end if

! ----- Generate name of datastructure to save in result datastructure
        call nameListLoad(listLoadResu)
        call copisd('LISTE_CHARGES', 'G', listLoad, listLoadResu)
    end if
!
end subroutine
