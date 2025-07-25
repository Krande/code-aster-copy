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
subroutine getBehaviourPara(l_mfront_proto, l_kit_thm, &
                            keywf, i_comp, algo_inte, &
                            iter_inte_maxi, resi_inte)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmdocv.h"
!
    aster_logical, intent(in) :: l_mfront_proto
    aster_logical, intent(in) :: l_kit_thm
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: i_comp
    character(len=16), intent(in) :: algo_inte
    integer(kind=8), pointer :: iter_inte_maxi
    real(kind=8), pointer :: resi_inte
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  l_mfront_proto     : .true. if MFront prototype
! In  l_kit_thm          : .true. if kit THM
! In  keywf              : factor keyword to read (COMPORTEMENT)
! In  i_comp             : factor keyword index
! In  algo_inte          : algorithm for integration of behaviour
! Out iter_inte_maxi     : value for ITER_INTE_MAXI
! Out resi_inte     : value for RESI_INTE
!
! --------------------------------------------------------------------------------------------------
!
    call nmdocv(keywf, i_comp, algo_inte, 'ITER_INTE_MAXI', &
                l_mfront_proto, l_kit_thm, vali=iter_inte_maxi)
    call nmdocv(keywf, i_comp, algo_inte, 'RESI_INTE     ', &
                l_mfront_proto, l_kit_thm, valr=resi_inte)
!
end subroutine
