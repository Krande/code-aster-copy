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
interface
    subroutine ace_affe_verif_elem(noma, jdme, lesmailles, ng, grpmail, grpnbma, &
                                   ace_nu, mclf, coderet, TFVCode)
    !
        character(len= 8), intent(in)         :: noma
        integer(kind=8),   intent(in)         :: ng, jdme, ace_nu
        integer(kind=8),   intent(inout)      :: lesmailles(*)
        character(len=24), intent(inout)      :: grpmail(*)
        integer(kind=8),   intent(inout)      :: grpnbma(*)
        character(len= *), intent(in)         :: mclf
        integer(kind=8),   intent(out)        :: coderet
        integer(kind=8), optional, intent(in) :: TFVCode
!
    end subroutine ace_affe_verif_elem
end interface
