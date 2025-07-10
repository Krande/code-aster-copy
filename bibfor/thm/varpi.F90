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
subroutine varpi(ds_thm, j_mater, p1, p1m, dp1, dp2, &
                 ep, surf, shut, &
                 phi0, dpi, sbjhm, &
                 wbjhm, epm, sbjh, wbjh)
!
    use THM_type
!
    implicit none
!
#include "asterfort/THM_type.h"
#include "asterfort/rcvala.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8), intent(in) :: j_mater
    real(kind=8), intent(in) :: p1, p1m, dp1, dp2
    real(kind=8), intent(in) :: phi0
    real(kind=8), intent(in) :: ep, surf, shut, sbjh, wbjh
    real(kind=8), intent(out) :: dpi
    real(kind=8), intent(out) :: sbjhm, wbjhm, epm

!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute the variation of the hydraulic pressure
!
! --------------------------------------------------------------------------------------------------
! In  j_mater          : coded material address
! In  p1m              : capillary pressure - At beginning of step
! In  p1               : capillary pressure - At end of current step!
! In  dp1              : increment of capillary pressure
! In  dp2              : increment of gaz pressure
! In  phi0             : initial porosity (THM_INIT)
! In  ep               : thickness of the adsorbed water layer
! In  surf             : specific surface of the material
! In  shut             : shuttleworth parameter
! In  sbjh             : saturated pores volume fraction from BJH - At end of current step
! In  wbjh             : unsaturatred pores surface fraction from BJH - At end of current step
! In  sbjhm            : saturated pores volume fraction from BJH - At beginning of step
! In  wbjhm            : unsaturatred pores surface fraction from BJH - At beginning of step
! In  epm              : thickness of the adsorbed water layer  - At beginning of step
! Out dpi              : variation of the hydraulic pressure at end of current time
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_para_bjh = 5
    real(kind=8) :: para_vale_bjh(nb_para_bjh)
    integer(kind=8) :: icodre_bjh(nb_para_bjh)
    character(len=16), parameter :: para_name_bjh(nb_para_bjh) = (/'A0     ', &
                                                                   'SHUTTLE', &
                                                                   'EPAI   ', &
                                                                   'S_BJH  ', &
                                                                   'W_BJH  '/)
! --------------------------------------------------------------------------------------------------

    dpi = 0.d0

! Value of sBJH and wbjh at beginning of step
    call rcvala(j_mater, ' ', 'THM_DIFFU', &
                1, 'PCAP', [p1m], &
                nb_para_bjh, para_name_bjh, para_vale_bjh, icodre_bjh, &
                1)
    sbjhm = para_vale_bjh(4)
    wbjhm = para_vale_bjh(5)
    epm = para_vale_bjh(3)

    dpi = dp2-(sbjh*p1)+(sbjhm*p1m)+((2./3.)*(0.5*(p1+p1m)*(sbjh-sbjhm))) &
          +((2./3.)*(surf/phi0)*(((wbjh*ep)+(wbjhm*epm))*0.5)*(-dp1)) &
          -((2./3.)*(surf/phi0)*shut*(wbjh-wbjhm))

end subroutine
