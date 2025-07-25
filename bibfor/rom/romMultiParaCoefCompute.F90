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
subroutine romMultiParaCoefCompute(ds_empi, ds_multipara, ds_algoGreedy, &
                                   i_mode_until, i_mode_coef, i_coef_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/romMultiParaCoefOneCompute.h"
!
    type(ROM_DS_Empi), intent(in) :: ds_empi
    type(ROM_DS_MultiPara), intent(inout) :: ds_multipara
    type(ROM_DS_AlgoGreedy), intent(in) :: ds_algoGreedy
    integer(kind=8), intent(in) :: i_mode_until, i_mode_coef
    integer(kind=8), optional, intent(in) :: i_coef_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Compute reduced coefficients
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_empi             : datastructure for empiric modes
! IO  ds_multipara        : datastructure for multiparametric problems
! In  ds_algoGreedy       : datastructure for Greedy algorithm
! In  i_mode_until        : last mode until compute
! In  i_mode_coef         : index of mode to compute coefficients
! In  i_coef_             : index of coefficient
!
! NB: if not(i_coef) : all coefficients for the mode
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_coef_b, i_coef_e, nb_mode, i_mode, i_coef, nb_vari_coef
    character(len=1) :: syst_2mbr_type
    character(len=24) :: coef_redu, syst_solu
    complex(kind=8), pointer :: vc_coef_redu(:) => null()
    complex(kind=8), pointer :: vc_syst_solu(:) => null()
    real(kind=8), pointer :: vr_coef_redu(:) => null()
    real(kind=8), pointer :: vr_syst_solu(:) => null()
    type(ROM_DS_Solve) :: ds_solveROM
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_45')
    end if
!
! - Get parameters
!
    coef_redu = ds_algoGreedy%coef_redu
    ds_solveROM = ds_algoGreedy%solveROM
    nb_vari_coef = ds_multipara%nb_vari_coef
    nb_mode = ds_solveROM%syst_size
    syst_2mbr_type = ds_solveROM%syst_2mbr_type
    syst_solu = ds_solveROM%syst_solu
!
! - Access to objects
!
    if (syst_2mbr_type .eq. 'R') then
        call jeveuo(syst_solu, 'L', vr=vr_syst_solu)
        call jeveuo(coef_redu, 'E', vr=vr_coef_redu)
    else if (syst_2mbr_type .eq. 'C') then
        call jeveuo(syst_solu, 'L', vc=vc_syst_solu)
        call jeveuo(coef_redu, 'E', vc=vc_coef_redu)
    else
        ASSERT(.false.)
    end if
    ASSERT(i_mode_until .le. nb_mode)
    ASSERT(i_mode_coef .le. nb_mode)
!
! - Select coefficient
!
    if (present(i_coef_)) then
        i_coef_b = i_coef_
        i_coef_e = i_coef_
    else
        i_coef_b = 1
        i_coef_e = nb_vari_coef
    end if
!
! - Compute reduced coefficients
!
    do i_coef = i_coef_b, i_coef_e
        if (niv .ge. 2) then
            call utmess('I', 'ROM2_46', si=i_coef)
        end if
! ----- Compute reduced coefficients for one evaluation of coefficient
        call romMultiParaCoefOneCompute(ds_empi, ds_multipara, ds_solveROM, &
                                        i_mode_until, i_mode_coef, i_coef)
! ----- Copy coefficients
        if (syst_2mbr_type .eq. 'R') then
            do i_mode = 1, i_mode_until
                vr_coef_redu(nb_vari_coef*(i_mode-1)+i_coef) = vr_syst_solu(i_mode)
            end do
        else if (syst_2mbr_type .eq. 'C') then
            do i_mode = 1, i_mode_until
                vc_coef_redu(nb_vari_coef*(i_mode-1)+i_coef) = vc_syst_solu(i_mode)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
end subroutine
