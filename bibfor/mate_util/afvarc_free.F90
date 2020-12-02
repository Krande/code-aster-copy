! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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

subroutine afvarc_free(varc_cata, varc_affe)
!
use Material_Datastructure_type
!
implicit none
!
    type(Mat_DS_VarcListCata), intent(inout) :: varc_cata
    type(Mat_DS_VarcListAffe), intent(inout) :: varc_affe
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Free allocated object
!
! --------------------------------------------------------------------------------------------------
!
! InOut varc_cata      : datastructure for catalog of external state variables
! InOut varc_affe      : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i
!
    do i = 1, varc_cata%nb_varc
        deallocate(varc_cata%list_cata_varc(i)%list_cmp)
    end do
!
    deallocate(varc_cata%list_cata_varc)
!
    if (varc_affe%nb_varc_cmp .ne. 0) then
        deallocate(varc_affe%list_affe_varc)
    endif
!
end subroutine
