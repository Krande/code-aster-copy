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
subroutine thmGetParaBJH(ds_thm, j_mater, p1)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
#include "asterfort/THM_type.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), intent(in) :: p1
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Evaluation of BJH Parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  j_mater      : coded material address
! In  p1           : capillary pressure - At end of current step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_para_bjh = 5
    real(kind=8) :: para_vale_bjh(nb_para_bjh)
    integer(kind=8) :: icodre_bjh(nb_para_bjh)
    character(len=16), parameter :: para_name_bjh(nb_para_bjh) = (/'A0     ', &
                                                                   'SHUTTLE', &
                                                                   'EPAI   ', &
                                                                   'S_BJH  ', &
                                                                   'W_BJH  '/)

    real(kind=8) :: ep, surf, sbjh, wbjh
!
! --------------------------------------------------------------------------------------------------
    para_vale_bjh(:) = 0.d0
    ep = 0.d0
    surf = 0.d0

    if ((ds_thm%ds_behaviour%rela_hydr) .eq. 'HYDR_TABBAL') then

        call rcvala(j_mater, ' ', 'THM_DIFFU', &
                    1, 'PCAP', [p1], &
                    nb_para_bjh, para_name_bjh, para_vale_bjh, icodre_bjh, &
                    1)

        surf = para_vale_bjh(1)
        ds_thm%ds_material%bjh%shuttle = para_vale_bjh(2)

        ep = para_vale_bjh(3)
        sbjh = para_vale_bjh(4)
        wbjh = para_vale_bjh(5)
!~         write (6,*) 'thmgetBJH1',ds_thm%ds_behaviour%rela_hydr

        if (surf .lt. 0.d0) then

            call utmess('F', 'THM1_95')

        else

            ds_thm%ds_material%bjh%A0 = surf

        end if

        if (ep .lt. 0.d0) then

            call utmess('F', 'THM1_96')

        else
            ds_thm%ds_material%bjh%epai = ep

        end if

        if (sbjh .lt. 0.d0 .or. sbjh .gt. 1.d0 .or. wbjh .lt. 0.d0 .or. wbjh .gt. 1.d0) then

            call utmess('F', 'THM1_97')

        else

            ds_thm%ds_material%bjh%SBJH = sbjh
            ds_thm%ds_material%bjh%wBJH = wbjh

        end if

    end if

end subroutine
