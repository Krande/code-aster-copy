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
subroutine romMultiParaCoefOneCompute(ds_empi, ds_multipara, ds_solveROM, &
                                      i_mode_until, i_mode_coef, i_coef)
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
#include "asterfort/romEvalCoef.h"
#include "asterfort/romMultiParaROMMatrCreate.h"
#include "asterfort/romMultiParaROM2mbrCreate.h"
#include "asterfort/romSolveROMSystSolve.h"
!
    type(ROM_DS_Empi), intent(in) :: ds_empi
    type(ROM_DS_MultiPara), intent(inout) :: ds_multipara
    type(ROM_DS_Solve), intent(in) :: ds_solveROM
    integer(kind=8), intent(in) :: i_mode_until, i_mode_coef, i_coef
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Compute reduced coefficients for one evaluation of coefficient
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_empi             : datastructure for empiric modes
! IO  ds_multipara        : datastructure for multiparametric problems
! In  ds_solveROM         : datastructure to solve systems (ROM)
! In  i_mode_until        : last mode until compute
! In  i_mode_coef         : index of mode to compute coefficients
! In  i_coef              : index of coefficient
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_mode, i_mode, vali(2)
    character(len=1) :: syst_2mbr_type
    character(len=24) :: syst_solu, syst_matr, syst_2mbr
    complex(kind=8), pointer :: vc_syst_solu(:) => null()
    complex(kind=8), pointer :: vc_syst_matr(:) => null()
    complex(kind=8), pointer :: vc_syst_2mbr(:) => null()
    real(kind=8), pointer :: vr_syst_solu(:) => null()
    real(kind=8), pointer :: vr_syst_matr(:) => null()
    real(kind=8), pointer :: vr_syst_2mbr(:) => null()
    real(kind=8) :: valr(2)
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_51', si=i_coef)
    end if
!
! - Get parameters
!
    nb_mode = ds_solveROM%syst_size
    syst_2mbr_type = ds_solveROM%syst_2mbr_type
    syst_matr = ds_solveROM%syst_matr
    syst_solu = ds_solveROM%syst_solu
    syst_2mbr = ds_solveROM%syst_2mbr
!
! - Access to objects
!
    if (syst_2mbr_type .eq. 'R') then
        call jeveuo(syst_solu, 'L', vr=vr_syst_solu)
    else if (syst_2mbr_type .eq. 'C') then
        call jeveuo(syst_solu, 'L', vc=vc_syst_solu)
    else
        ASSERT(ASTER_FALSE)
    end if
    ASSERT(i_mode_until .le. nb_mode)
    ASSERT(i_mode_coef .le. nb_mode)
!
! - Initialization of matrix
!
    if (syst_2mbr_type .eq. 'R') then
        call jeveuo(syst_matr, 'E', vr=vr_syst_matr)
        vr_syst_matr(1:nb_mode*nb_mode) = 0.d0
    else if (syst_2mbr_type .eq. 'C') then
        call jeveuo(syst_matr, 'E', vc=vc_syst_matr)
        vc_syst_matr(1:nb_mode*nb_mode) = dcmplx(0.d0, 0.d0)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Initialization of vector
!
    if (syst_2mbr_type .eq. 'R') then
        call jeveuo(syst_2mbr, 'E', vr=vr_syst_2mbr)
        vr_syst_2mbr(1:nb_mode) = 0.d0
    else if (syst_2mbr_type .eq. 'C') then
        call jeveuo(syst_2mbr, 'E', vc=vc_syst_2mbr)
        vc_syst_2mbr(1:nb_mode) = dcmplx(0.d0, 0.d0)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Evaluate coefficients
!
    call romEvalCoef(ds_multipara, l_init=ASTER_FALSE, &
                     i_mode_coef_=i_mode_coef, i_coef_=i_coef)
!
! - Compute reduced second member
!
    call romMultiParaROM2mbrCreate(ds_empi, ds_multipara, i_coef, syst_2mbr)
!
! - Compute reduced matrix
!
    call romMultiParaROMMatrCreate(ds_empi, ds_multipara, i_coef, syst_matr)
!
! - Solve reduced system
!
    call romSolveROMSystSolve(ds_solveROM, size_to_solve_=i_mode_until)
!
! - Debug print
!
    if (niv .ge. 2) then
        if (syst_2mbr_type .eq. 'R') then
            do i_mode = 1, i_mode_until
                valr(1) = vr_syst_solu(i_mode)
                valr(2) = 0.d0
                vali(1) = i_mode
                vali(2) = i_coef
                call utmess('I', 'ROM2_52', ni=2, vali=vali, nr=2, valr=valr)
            end do
        else if (syst_2mbr_type .eq. 'C') then
            do i_mode = 1, i_mode_until
                valr(1) = real(vc_syst_solu(i_mode))
                valr(2) = dimag(vc_syst_solu(i_mode))
                vali(1) = i_mode
                vali(2) = i_coef
                call utmess('I', 'ROM2_52', ni=2, vali=vali, nr=2, valr=valr)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
end subroutine
