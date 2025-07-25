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

subroutine romSolveROMSystCreate(syst_matr_type, syst_2mbr_type, syst_type, &
                                 nb_mode, ds_solve)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=1), intent(in) :: syst_matr_type
    character(len=1), intent(in) :: syst_2mbr_type
    character(len=1), intent(in) :: syst_type
    integer(kind=8), intent(in) :: nb_mode
    type(ROM_DS_Solve), intent(inout) :: ds_solve
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Create objects to solve system (ROM)
!
! --------------------------------------------------------------------------------------------------
!
! In  syst_matr_type   : type of matrix (real or complex)
! In  syst_2mbr_type   : type of second member (real or complex)
! In  syst_type        : global type of system (real or complex)
! In  nb_mode          : number of empiric modes
! IO  ds_solve         : datastructure to solve systems
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: syst_matr, syst_2mbr, vect_zero, syst_solu
    integer(kind=8) :: jv_dummy
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM2_33')
    end if
!
! - Get parameters
!
    syst_matr = ds_solve%syst_matr
    syst_2mbr = ds_solve%syst_2mbr
    vect_zero = ds_solve%vect_zero
    syst_solu = ds_solve%syst_solu
!
! - Create objects
!
    call wkvect(syst_matr, 'V V '//syst_matr_type, nb_mode*nb_mode, jv_dummy)
    call wkvect(syst_2mbr, 'V V '//syst_2mbr_type, nb_mode, jv_dummy)
    call wkvect(syst_solu, 'V V '//syst_type, nb_mode, jv_dummy)
    call wkvect(vect_zero, 'V V '//syst_type, nb_mode, jv_dummy)
!
! - Save parameters
!
    ds_solve%syst_size = nb_mode
    ds_solve%syst_type = syst_type
    ds_solve%syst_matr_type = syst_matr_type
    ds_solve%syst_2mbr_type = syst_2mbr_type
!
end subroutine
