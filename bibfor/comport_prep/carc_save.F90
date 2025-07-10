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
subroutine carc_save(mesh, carcri, prepMapCarcri)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/setBehaviourParaValue.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: carcri
    type(BehaviourPrep_MapCarcri), intent(in) :: prepMapCarcri
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
! In  prepMapCarcri    : datastructure to construct CARCRI map
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer(kind=8), parameter :: nbCmp = CARCRI_SIZE
    character(len=24), parameter :: list_elem_affe = '&&CARCSAVE.LIST'
    aster_logical :: l_affe_all, l_parallel_mesh
    integer(kind=8) :: nb_elem_affe
    integer(kind=8), pointer :: v_elem_affe(:) => null()
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword
    real(kind=8), pointer :: carcriValv(:) => null()
    real(kind=8) :: parm_theta_thm, parm_alpha_thm
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCarcri%nb_comp
    l_parallel_mesh = isParallelMesh(mesh)

! - Access to MAP
    call jeveuo(carcri//'.VALV', 'E', vr=carcriValv)

! - Get parameters from SCHEMA_THM
    parm_theta_thm = prepMapCarcri%parm_theta_thm
    parm_alpha_thm = prepMapCarcri%parm_alpha_thm

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get list of elements where comportment is defined
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            list_elem_affe, l_affe_all, nb_elem_affe)

! ----- Set in <CARTE>
        call setBehaviourParaValue(prepMapCarcri%prepCrit, &
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
