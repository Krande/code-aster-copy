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

subroutine vechme_wrap(stop, modelz, lload_namez, lload_infoz, &
                       time_value, time_delta, time_theta, &
                       cara_elem, mate, mateco, vect_elemz, varc_currz)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/vechme.h"
!
! person_in_charge: jacques.pellet at edf.fr
!
    character(len=1), intent(in) :: stop
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: lload_namez
    character(len=*), intent(in) :: lload_infoz
    real(kind=8), intent(in) :: time_value, time_delta, time_theta
    character(len=*), intent(in) :: cara_elem
    character(len=*), intent(in) :: mate, mateco
    character(len=*), intent(inout) :: vect_elemz
    character(len=*), intent(in) :: varc_currz
!
! --------------------------------------------------------------------------------------------------
!
! Compute Neumann loads
!
! Dead and fixed loads
!
! --------------------------------------------------------------------------------------------------
!
! In  stop           : continue or stop computation if no loads on elements
! In  model          : name of model
! In  mate           : name of material characteristics (field)
! In  cara_elem      : name of elementary characteristics (field)
! In  lload_name     : name of object for list of loads name
! In  lload_info     : name of object for list of loads info
! In  inst           : times informations
! In  varc_curr      : command variable for current time
! IO  vect_elem      : name of vect_elem result
!
! ATTENTION :
!   LE VECT_ELEM (VECELZ) RESULTAT A 1 PARTICULARITE :
!   CERTAINS RESUELEM NE SONT PAS DES RESUELEM MAIS DES CHAM_NO (.VEASS)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: inst(3)
    inst(1) = time_value
    inst(2) = time_delta
    inst(3) = time_theta
    call vechme(stop, &
                modelz, cara_elem, mate, mateco, &
                lload_namez, lload_infoz, &
                inst, &
                vect_elemz, &
                varcCurrZ_=varc_currz, jvBase_='G')
!
end subroutine
