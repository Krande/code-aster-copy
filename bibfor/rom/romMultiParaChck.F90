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
subroutine romMultiParaChck(ds_multipara, l_stab_fsi)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/utmess.h"
#include "asterfort/romFieldChck.h"
!
    type(ROM_DS_MultiPara), intent(in) :: ds_multipara
    aster_logical, intent(in) :: l_stab_fsi
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Check data for multiparametric problems
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_multipara     : datastructure for multiparametric problems
! In  l_stab_fsi       : flag for IFS stabilization
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: syme
    integer(kind=8) :: nb_matr, nb_vect, nb_vari_para, nb_vale_para, ndim
    integer(kind=8) :: i_matr, i_vect, i_vari_para
    character(len=24) :: nume_dof_ref, nume_dof
    character(len=24) :: numeq_ref, numeq
!
! --------------------------------------------------------------------------------------------------
!
    nb_matr = ds_multipara%nb_matr
    nb_vect = ds_multipara%nb_vect
    nb_vari_para = ds_multipara%nb_vari_para
!
! - Check numbering in matrix
!
    if (ds_multipara%matr_name(1) .ne. ' ') then
        call dismoi('NOM_NUME_DDL', ds_multipara%matr_name(1), 'MATR_ASSE', repk=nume_dof_ref)
        do i_matr = 2, nb_matr
            call dismoi('NOM_NUME_DDL', ds_multipara%matr_name(i_matr), 'MATR_ASSE', repk=nume_dof)
            if (nume_dof .ne. nume_dof_ref) then
                call utmess('F', 'ROM2_21')
            end if
        end do
    end if
!
! - Check numbering in vector
!
    if (ds_multipara%vect_name(1) .ne. ' ') then
        call dismoi('NUME_EQUA', ds_multipara%vect_name(1), 'CHAM_NO', repk=numeq_ref)
        do i_vect = 2, nb_vect
            call dismoi('NUME_EQUA', ds_multipara%vect_name(i_vect), 'CHAM_NO', repk=numeq)
            if (numeq .ne. numeq_ref) then
                call utmess('F', 'ROM2_21')
            end if
        end do
    end if
!
! - Only symmetric matrix
!
    if (ds_multipara%matr_name(1) .ne. ' ') then
        do i_matr = 1, nb_matr
            call dismoi('TYPE_MATRICE', ds_multipara%matr_name(i_matr), 'MATR_ASSE', repk=syme)
            if (syme .eq. 'NON_SYM') then
                call utmess('F', 'ROM2_22')
            end if
        end do
    end if
!
! - If no definition of variation => no functions !
!
    if (nb_vari_para .eq. 0) then
        do i_matr = 1, nb_matr
            if (ds_multipara%matr_coef(i_matr)%l_func) then
                call utmess('F', 'ROM2_25', sk=ds_multipara%matr_name(i_matr))
            end if
        end do
        do i_vect = 1, nb_vect
            if (ds_multipara%vect_coef(i_vect)%l_func) then
                call utmess('F', 'ROM2_31')
            end if
        end do
    end if
!
! - Same number of values for each parameter
!
    if (nb_vari_para .ne. 0) then
        nb_vale_para = ds_multipara%vari_para(1)%nb_vale_para
        do i_vari_para = 2, nb_vari_para
            if (ds_multipara%vari_para(i_vari_para)%nb_vale_para .ne. nb_vale_para) then
                call utmess('F', 'ROM2_29')
            end if
        end do
    end if
!
! - Check components PRES and PHI are in the model if we active l_stab_fsi
!
    if (l_stab_fsi) then
        call dismoi('DIM_GEOM', ds_multipara%field%model, 'MODELE', repi=ndim)
        if (ndim .eq. 2) then
            call romFieldChck(ds_multipara%field, fieldName_='UPPHI_2D')
        elseif (ndim .eq. 3) then
            call romFieldChck(ds_multipara%field, fieldName_='UPPHI_3D')
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
end subroutine
