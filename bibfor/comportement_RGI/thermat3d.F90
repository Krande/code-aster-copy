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
subroutine thermat3d(teta1, nrjm, tetas, tetar, dt80, &
                     dth0, dth, cthp, cthv)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!   influence de la temperature sur les parametres de fluage
!
!   declarations externes
    implicit none
#include "rgi_module.h"
!
!   variables externes
    real(kind=8), intent(in) :: teta1, nrjm, tetas, tetar, dt80, dth0
    real(kind=8), intent(out) :: dth, cthp, cthv
!
!   variables locales
    real(kind=8) :: easurrm, unsurtr, unsurts, unsurt, unsurt80
    real(kind=8) :: xxx180, xxx1, xxx2, xxx280, ath, cth80
!
!**********************************************************************
!   reglage des activation thermiques pour le fluage
!   +++ les teta sont en degres Celsius +++
!**********************************************************************
!
!   calcul du terme d activation d Arrhenius pour le potentiel
    easurrm = nrjm/KGAZP
!   la temperature de reference est 20°C
    unsurtr = 1.d0/(tetar+ZEROKLV)
!   unsurts la temperature de seuil pour la modif du potentiel
    unsurts = 1.d0/(tetas+ZEROKLV)
!   calcul des coeff d activation thermique
    unsurt = (1.d0/(teta1+ZEROKLV))
!   cas de l eau
    cthv = exp(-EASURRW*(unsurt-unsurtr))
!   cas de l endommagement thermique
    xxx1 = unsurts-unsurt
    xxx2 = 0.5d0*(xxx1+dabs(xxx1))
    cthp = exp(easurrm*xxx2)
!
!
!************************************************************************
!   endommagement thermique
!************************************************************************
    if (tetas .lt. T80DEG) then
        unsurt80 = 1.d0/(ZEROKLV+T80DEG)
        xxx180 = unsurts-unsurt80
        xxx280 = 0.5d0*(xxx180+dabs(xxx180))
        cth80 = exp(easurrm*xxx280)
        ath = 1.d0/(cth80-1.d0)*(dt80/(1.d0-dt80))
        dth = 1.d0-1.d0/(1.d0+ath*(cthp-1.d0))
    else
        dth = 0.d0
    end if
    dth = dmax1(dth0, dth)
!***********************************************************************!
end subroutine
