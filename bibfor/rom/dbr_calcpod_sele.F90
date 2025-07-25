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
subroutine dbr_calcpod_sele(nb_mode_maxi, tole_svd, s, nb_sing, nb_mode)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: nb_mode_maxi
    real(kind=8), intent(in) :: tole_svd
    real(kind=8), pointer :: s(:)
    integer(kind=8), intent(in) :: nb_sing
    integer(kind=8), intent(out) :: nb_mode
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE - Compute
!
! Select singular vectors
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_mode_maxi     : maximum number of emprical modes
! In  tole_svd         : tolerance for SVD
! In  nb_sing          : total number of singular values
! In  s                : singular values
! Out nb_mode          : number of modes selected
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: vale_mini, vale_maxi, vale_tole, valr(2)
    integer(kind=8) :: i_sing, vali(2)
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM5_4')
    end if
!
! - Init
!
    nb_mode = 0
!
! - Get parameters
!
    vale_mini = s(nb_sing)
    vale_maxi = s(1)
    vale_tole = tole_svd*vale_maxi
!
! - Select singular values
!
    if (nb_mode_maxi .eq. 0) then
        do i_sing = 1, nb_sing
            if (s(i_sing) .ge. vale_tole) then
                nb_mode = nb_mode+1
            end if
        end do
    else
        if (nb_sing .le. nb_mode_maxi) then
            nb_mode = nb_sing
        else
            nb_mode = nb_mode_maxi
        end if
    end if
    valr(1) = vale_mini
    valr(2) = vale_maxi
    vali(1) = nb_sing
    vali(2) = nb_mode
    call utmess('I', 'ROM5_5', ni=2, vali=vali, nr=2, valr=valr)
!
    if (nb_mode .lt. 1) then
        call utmess('F', 'ROM5_6')
    end if
!
end subroutine
