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

!
!
interface 
    subroutine romAlgoNLCorrEFMecaResidual(v_cnequi, ds_algorom, l_cine, v_ccid,&
                                           resi)
        use Rom_Datastructure_type
        real(kind=8), pointer :: v_cnequi(:)
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        aster_logical, intent(in) :: l_cine
        integer(kind=8), pointer :: v_ccid(:)
        real(kind=8), intent(out) :: resi
    end subroutine romAlgoNLCorrEFMecaResidual
end interface
