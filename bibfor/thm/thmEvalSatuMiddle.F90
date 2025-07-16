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
subroutine thmEvalSatuMiddle(ds_thm, j_mater, p1, temp, &
                             satur, dsatur, retcom)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), intent(in) :: p1, temp
    real(kind=8), intent(out) :: satur, dsatur
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Evaluation of "middle" saturation (only LIQU_VAPE)
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  j_mater          : coded material address
! In  p1               : capillary pressure - At end of current step
! Out satur            : saturation
! Out dsatur           : derivative of saturation (/pc)
! Out retcom           : return code for error
!                         2 - If saturation doesn't belon to ]0,1[
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_para = 2
    real(kind=8) :: para_vale(nb_para)
    integer(kind=8) :: icodre(nb_para)
    character(len=16), parameter :: para_name(nb_para) = (/'SATU_PRES  ', 'D_SATU_PRES'/)
    character(len=16), parameter :: npar(2) = (/'PCAP', 'TEMP'/)
!
! --------------------------------------------------------------------------------------------------
!
    retcom = 0
    para_vale(:) = 0.d0
    if (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_UTIL' .or. &
        ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_ENDO' .or. &
        ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_TABBAL') then
        call rcvala(j_mater, ' ', 'THM_DIFFU', &
                    2, npar, [p1, temp], &
                    nb_para, para_name, para_vale, icodre, &
                    1)
        satur = para_vale(1)
        dsatur = para_vale(2)
        ASSERT(ds_thm%ds_behaviour%satur_type .eq. SATURATED_SPEC)
    else
        ASSERT(.false.)
    end if
    if (satur .gt. 1.d0 .or. satur .lt. 0.d0) then
        retcom = 2
    end if
!
end subroutine
