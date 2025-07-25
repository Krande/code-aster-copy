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
subroutine romEvalCoef(ds_multipara, l_init, i_mode_coef_, i_coef_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/romCoefInfo.h"
#include "asterfort/romEvalCoefFunc.h"
#include "asterfort/romEvalCoefPrep.h"
!
    type(ROM_DS_MultiPara), intent(inout) :: ds_multipara
    aster_logical, intent(in) :: l_init
    integer(kind=8), optional, intent(in) :: i_mode_coef_
    integer(kind=8), optional, intent(in) :: i_coef_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Evaluation of coefficients
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_multipara     : datastructure for multiparametric problems
! In  l_init           : .true. if first evaluation
! In  i_mode_coef      : index of mode to compute coefficients
! In  i_coef           : index of coefficient
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_matr, i_vect, nb_matr, nb_vect, i_coef_list, nb_vari_coef
    integer(kind=8) ::  i_coef, i_mode_coef
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        if (l_init) then
            call utmess('I', 'ROM5_93')
        else
            call utmess('I', 'ROM5_44', si=i_coef_)
        end if
    end if
!
! - Get parameters
!
    nb_vari_coef = ds_multipara%nb_vari_coef
!
! - Get index in global list of coefficients
!
    if (l_init) then
        i_coef = 1
        i_mode_coef = 0
        i_coef_list = 0
    else
        ASSERT(present(i_coef_))
        ASSERT(present(i_mode_coef_))
        i_coef = i_coef_
        i_mode_coef = i_mode_coef_
        i_coef_list = nb_vari_coef*(i_mode_coef-1)+i_coef
    end if
!
! - Prepare EVALCOEF datastructure
!
    call romEvalCoefPrep(i_coef_list, ds_multipara)
!
! - Evaluate coefficients for matrix
!
    nb_matr = ds_multipara%nb_matr
    do i_matr = 1, nb_matr
        call romEvalCoefFunc(ds_multipara%evalcoef, ds_multipara%matr_coef(i_matr), i_coef)
        if (niv .ge. 2) then
            call romCoefInfo('M', &
                             ds_multipara%matr_name(i_matr), &
                             i_coef, &
                             ds_multipara%matr_coef(i_matr))
        end if
    end do
!
! - Evaluate coefficients for second member
!
    nb_vect = ds_multipara%nb_vect
    do i_vect = 1, nb_vect
        call romEvalCoefFunc(ds_multipara%evalcoef, ds_multipara%vect_coef(i_vect), i_coef_list)
        if (niv .ge. 2) then
            call romCoefInfo('V', &
                             ds_multipara%vect_name(i_vect), &
                             i_coef, &
                             ds_multipara%vect_coef(i_vect))
        end if
    end do
!
end subroutine
