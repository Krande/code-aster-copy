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
subroutine setBehaviourParaValue(behaviourCrit, parm_theta_thm, parm_alpha_thm, &
                                 iFactorKeyword_, carcriList_, carcriMap_)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/setMFrontPara.h"
!
    type(Behaviour_Crit), pointer :: behaviourCrit(:)
    real(kind=8), intent(in) :: parm_theta_thm, parm_alpha_thm
    integer, optional, intent(in) :: iFactorKeyword_
    real(kind=8), intent(out), optional :: carcriList_(:)
    real(kind=8), pointer, optional :: carcriMap_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Set values in the map or in list
!
! --------------------------------------------------------------------------------------------------
!
! In  behaviourCrit    : parameters for integration of constitutive law
! In  iFactorKeyword   : index of factor keyword (for map)
! In  carcriList       : list for parameters for integration of constitutive law
! In  carcriMap        : map for parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
    iFactorKeyword = 1
    if (present(iFactorKeyword_)) then
        iFactorKeyword = iFactorKeyword_
    end if
!
    if (present(carcriMap_)) then
        carcriMap_(ITER_INTE_MAXI) = behaviourCrit(iFactorKeyword)%iter_inte_maxi
        carcriMap_(TYPE_MATR_T) = behaviourCrit(iFactorKeyword)%type_matr_t
        carcriMap_(RESI_INTE_RELA) = behaviourCrit(iFactorKeyword)%resi_inte_rela
        carcriMap_(PARM_THETA) = behaviourCrit(iFactorKeyword)%parm_theta
        carcriMap_(ITER_INTE_PAS) = behaviourCrit(iFactorKeyword)%iter_inte_pas
        carcriMap_(ALGO_INTE_R) = behaviourCrit(iFactorKeyword)%algo_inte_r
        carcriMap_(VALE_PERT_RELA) = behaviourCrit(iFactorKeyword)%vale_pert_rela
        carcriMap_(RESI_DEBORST_MAX) = behaviourCrit(iFactorKeyword)%resi_deborst_max
        carcriMap_(ITER_DEBORST_MAX) = behaviourCrit(iFactorKeyword)%iter_deborst_max
        carcriMap_(RESI_RADI_RELA) = behaviourCrit(iFactorKeyword)%resi_radi_rela
        carcriMap_(IVARIEXT1) = behaviourCrit(iFactorKeyword)%jvariext1
        carcriMap_(IVARIEXT2) = behaviourCrit(iFactorKeyword)%jvariext2
        carcriMap_(PARM_THETA_THM) = parm_theta_thm
        carcriMap_(PARM_ALPHA_THM) = parm_alpha_thm
        carcriMap_(IPOSTITER) = behaviourCrit(iFactorKeyword)%ipostiter
        if (behaviourCrit(iFactorKeyword)%l_matr_unsymm) then
            carcriMap_(CARCRI_MATRSYME) = 1
        else
            carcriMap_(CARCRI_MATRSYME) = 0
        end if
        carcriMap_(IPOSTINCR) = behaviourCrit(iFactorKeyword)%ipostincr
! ----- For external solvers (UMAT / MFRONT)
        carcriMap_(EXTE_PTR) = behaviourCrit(iFactorKeyword)%extern_ptr
        carcriMap_(EXTE_TYPE) = behaviourCrit(iFactorKeyword)%extern_type
        carcriMap_(EXTE_STRAIN) = behaviourCrit(iFactorKeyword)%exte_strain
    end if
    if (present(carcriList_)) then
        carcriList_(ITER_INTE_MAXI) = behaviourCrit(iFactorKeyword)%iter_inte_maxi
        carcriList_(TYPE_MATR_T) = behaviourCrit(iFactorKeyword)%type_matr_t
        carcriList_(RESI_INTE_RELA) = behaviourCrit(iFactorKeyword)%resi_inte_rela
        carcriList_(PARM_THETA) = behaviourCrit(iFactorKeyword)%parm_theta
        carcriList_(ITER_INTE_PAS) = behaviourCrit(iFactorKeyword)%iter_inte_pas
        carcriList_(ALGO_INTE_R) = behaviourCrit(iFactorKeyword)%algo_inte_r
        carcriList_(VALE_PERT_RELA) = behaviourCrit(iFactorKeyword)%vale_pert_rela
        carcriList_(RESI_DEBORST_MAX) = behaviourCrit(iFactorKeyword)%resi_deborst_max
        carcriList_(ITER_DEBORST_MAX) = behaviourCrit(iFactorKeyword)%iter_deborst_max
        carcriList_(RESI_RADI_RELA) = behaviourCrit(iFactorKeyword)%resi_radi_rela
        carcriList_(IVARIEXT1) = behaviourCrit(iFactorKeyword)%jvariext1
        carcriList_(IVARIEXT2) = behaviourCrit(iFactorKeyword)%jvariext2
        carcriList_(PARM_THETA_THM) = parm_theta_thm
        carcriList_(PARM_ALPHA_THM) = parm_alpha_thm
        carcriList_(IPOSTITER) = behaviourCrit(iFactorKeyword)%ipostiter
        if (behaviourCrit(iFactorKeyword)%l_matr_unsymm) then
            carcriList_(CARCRI_MATRSYME) = 1
        else
            carcriList_(CARCRI_MATRSYME) = 0
        end if
        carcriList_(IPOSTINCR) = behaviourCrit(iFactorKeyword)%ipostincr
! ----- For external solvers (UMAT / MFRONT)
        carcriList_(EXTE_PTR) = behaviourCrit(iFactorKeyword)%extern_ptr
        carcriList_(EXTE_TYPE) = behaviourCrit(iFactorKeyword)%extern_type
        carcriList_(EXTE_STRAIN) = behaviourCrit(iFactorKeyword)%exte_strain
    end if

! - Set values for MFRONT
    call setMFrontPara(behaviourCrit, iFactorKeyword)
!
end subroutine
