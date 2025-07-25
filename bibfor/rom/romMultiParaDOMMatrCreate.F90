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
subroutine romMultiParaDOMMatrCreate(ds_multipara, i_coef, ds_solve)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/utmess.h"
#include "asterfort/infniv.h"
!
    type(ROM_DS_MultiPara), intent(in) :: ds_multipara
    integer(kind=8), intent(in) :: i_coef
    type(ROM_DS_Solve), intent(in) :: ds_solve
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Create matrix for multiparametric problems (complete problem)
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_multipara     : datastructure for multiparametric problems
! In  i_coef           : index of coefficient
! In  ds_solve         : datastructure to solve systems
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: nb_matr_maxi = 8
    character(len=1) :: type_comb(nb_matr_maxi)
    real(kind=8) :: coef_comb(2*nb_matr_maxi)
    character(len=24) :: matr_comb(nb_matr_maxi)
    integer(kind=8) :: i_coef_comb, i_matr, nb_matr
    aster_logical :: l_coefm_cplx
    character(len=19) :: syst_matr
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM5_63')
    end if
!
! - Initializations
!
    syst_matr = ds_solve%syst_matr
    nb_matr = ds_multipara%nb_matr
    ASSERT(nb_matr .le. nb_matr_maxi)
    type_comb(:) = ' '
    matr_comb(:) = ' '
    coef_comb(:) = 0.d0
!
! - Compute matrix
!
    i_coef_comb = 0
    do i_matr = 1, nb_matr
        matr_comb(i_matr) = ds_multipara%matr_name(i_matr)
        l_coefm_cplx = ds_multipara%matr_coef(i_matr)%l_cplx
        if (l_coefm_cplx) then
            type_comb(i_matr) = 'C'
            i_coef_comb = i_coef_comb+1
            coef_comb(i_coef_comb) = real(ds_multipara%matr_coef(i_matr)%coef_cplx(i_coef))
            i_coef_comb = i_coef_comb+1
            coef_comb(i_coef_comb) = dimag(ds_multipara%matr_coef(i_matr)%coef_cplx(i_coef))
        else
            type_comb(i_matr) = 'R'
            i_coef_comb = i_coef_comb+1
            coef_comb(i_coef_comb) = ds_multipara%matr_coef(i_matr)%coef_real(i_coef)
        end if
    end do
    call mtcmbl(nb_matr, type_comb, coef_comb, matr_comb, syst_matr, &
                ' ', ' ', 'ELIM=')
!
end subroutine
