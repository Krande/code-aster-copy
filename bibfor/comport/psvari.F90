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

subroutine psvari(rela_comp, nbvari, ipop1, ipop2)
    implicit none
#include "asterfort/assert.h"
    character(len=16), intent(in) :: rela_comp
    integer(kind=8), intent(in) :: nbvari
    integer(kind=8), intent(out) :: ipop1, ipop2
!
!     FONCTION REALISEE :
!
!     PERMET DE CONNAITRE EN FONCTION DE LA RELATION DE COMPORTEMENT
!     PARMI LES VARIABLES INTERNES LA POSITION DE :
!
!         - LA DEFORMATION PLASTIQUE CUMULEE
!         - L'INDICATEUR DE PLASTICITE
!
! ENTREE  --->  COMPOR : NOM DE LA RELATION DE COMPORTEMENT
!         --->  NBVARI : NOMBRE DE VARIABLES INTERNES
!
! SORTIE
!         --->  IPOS1  : POSITION DE LA DEFORMATION PLASTIQUE CUMULEE
!         --->  IPOS2  : POSITION DE L'INDICATEUR DE PLASTICITE
!
!     ------------------------------------------------------------------
!
    ipop1 = 0
    ipop2 = 0
    if ((rela_comp .eq. 'LEMAITRE') .or. (rela_comp .eq. 'VMIS_ECMI_TRAC') .or. &
        (rela_comp .eq. 'VMIS_ECMI_LINE') .or. (rela_comp .eq. 'VMIS_CIN1_CHAB') .or. &
        (rela_comp .eq. 'VMIS_CIN2_CHAB') .or. (rela_comp .eq. 'VISC_CIN1_CHAB') .or. &
        (rela_comp .eq. 'VISC_CIN2_CHAB') .or. (rela_comp .eq. 'VMIS_ISOT_TRAC') .or. &
        (rela_comp .eq. 'VMIS_ISOT_LINE') .or. (rela_comp .eq. 'VISC_ISOT_TRAC') .or. &
        (rela_comp .eq. 'VISC_ISOT_LINE')) then
        ipop1 = 1
        ipop2 = 2
    else if ((rela_comp .eq. 'ROUSSELIER')) then
        ipop1 = 1
        ipop2 = 3
    else if ((rela_comp .eq. 'ROUSS_PR') .or. (rela_comp .eq. 'ROUSS_VISC')) &
        then
        ipop1 = 1
        ipop2 = nbvari
    else if (rela_comp .eq. 'MONOCRISTAL') then
        ipop1 = nbvari-1
        ipop2 = nbvari
    else if (rela_comp .eq. 'POLYCRISTAL') then
        ipop1 = 7
        ipop2 = nbvari
    else
        ASSERT(ASTER_FALSE)
!
    end if
!
end subroutine
