! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
!
subroutine verins(ds_posttimestep,instin)
use NonLin_Datastructure_type
implicit none
#include "asterfort/utmess.h"
    type(NL_DS_PostTimeStep) :: ds_posttimestep
    real(kind=8) :: instin, valr(1)
    type(NL_DS_SelectList) :: selectList
    integer :: iinst
!
! Vérifier les instants demandés par MODE_VIBR sont postérieurs 
! à instant initial de DYNA_NON_LINE 
!
        if (ds_posttimestep%l_mode_vibr) then
            selectList = ds_posttimestep%mode_vibr%selector
            valr(1) = instin
            do iinst = 1, selectList%nb_value
                if (selectList%list_value(iinst) <= instin) then
                    call utmess('A', 'UTILITAI8_75',nr=1,valr=valr)
                endif
            end do
        endif  
end subroutine