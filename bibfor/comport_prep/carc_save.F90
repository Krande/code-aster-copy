! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine carc_save(mesh, carcri, ds_compor_para)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lcdiscard.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/setBehaviourParaValue.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: carcri
    type(Behaviour_PrepCrit), intent(in) :: ds_compor_para
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Save informations in <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  carcri           : name of <CARTE> CARCRI
! In  ds_compor_para   : datastructure to prepare parameters for constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer, parameter :: nbCmp = CARCRI_SIZE
    character(len=24), parameter :: list_elem_affe = '&&CARCSAVE.LIST'
    aster_logical :: l_affe_all, l_parallel_mesh
    integer :: nb_elem_affe
    integer, pointer :: v_elem_affe(:) => null()
    integer :: iFactorKeyword, nbFactorKeyword
    real(kind=8), pointer :: carcriValv(:) => null()
    real(kind=8) :: parm_theta_thm, parm_alpha_thm
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = ds_compor_para%nb_comp
    l_parallel_mesh = isParallelMesh(mesh)

! - Access to MAP
    call jeveuo(carcri//'.VALV', 'E', vr=carcriValv)

! - Get parameters from SCHEMA_THM
    parm_theta_thm = ds_compor_para%parm_theta_thm
    parm_alpha_thm = ds_compor_para%parm_alpha_thm

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get list of elements where comportment is defined
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            list_elem_affe, l_affe_all, nb_elem_affe)

! ----- Set in <CARTE>
        call setBehaviourParaValue(ds_compor_para%v_crit, &
                                   parm_theta_thm, parm_alpha_thm, &
                                   iFactorKeyword, carcriMap_=carcriValv)

! ----- Affect in <CARTE>
        if (l_affe_all) then
            call nocart(carcri, 1, nbCmp)
        else
            if (nb_elem_affe > 0 .or. .not. l_parallel_mesh) then
                call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
                call nocart(carcri, 3, nbCmp, mode='NUM', nma=nb_elem_affe, &
                            limanu=v_elem_affe)
            end if
            call jedetr(list_elem_affe)
        end if
    end do
!
    call jedetr(carcri//'.NCMP')
!
end subroutine
