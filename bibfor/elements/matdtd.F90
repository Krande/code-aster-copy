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

subroutine matdtd(nomte, testl1, testl2, dsidep, cisail, &
                  x3, cour, r, cosa, kappa, &
                  dtildi)
    implicit none
!
#include "asterf_types.h"
#include "asterfort/r8inir.h"
    character(len=16) :: nomte
    real(kind=8) :: cisail, x3, cour, r, cosa, rhos, rhot, kappa
    real(kind=8) :: dsidep(6, 6), mata(3, 3), dtildi(5, 5)
    aster_logical :: testl1, testl2
!
!-----------------------------------------------------------------------
    real(kind=8) :: rhos2, rhost, rhot2, x32
!-----------------------------------------------------------------------
    call r8inir(25, 0.d0, dtildi, 1)
!
!     CALCULS DE LA MATRICE DTILDI
!
!
!   EXTRACTION DE DSIDEP LA SOUS MATRICE MATA CONCERNEE
!   ON NOTE QUE LA MATRICE DSIDEP OBTENUE EST DE TYPE CONTRAINTES PLANES
!   LES COMPOSANTES CORRESPONDANT AU TERME DE CISAILLEMENT SONT REMPLIES
!   "ARTIFICIELLEMENT"
!
    if (nomte .eq. 'MECXSE3') then
!
        mata(1, 1) = dsidep(1, 1)
        mata(1, 2) = dsidep(1, 2)
        mata(1, 3) = 0.d0
        mata(2, 1) = dsidep(2, 1)
        mata(2, 2) = dsidep(2, 2)
        mata(2, 3) = 0.d0
        mata(3, 1) = 0.d0
        mata(3, 2) = 0.d0
        mata(3, 3) = cisail*kappa/2.d0
!
        if (testl1) then
            rhos = 1.d0
        else
            rhos = 1.d0+x3*cour
        end if
        if (testl2) then
            rhot = 1.d0
        else
            rhot = 1.d0+x3*cosa/r
        end if
!
        x32 = x3*x3
        rhos2 = rhos*rhos
        rhot2 = rhot*rhot
        rhost = rhos*rhot
!
        dtildi(1, 1) = mata(1, 1)/rhos2
        dtildi(1, 2) = mata(1, 1)*x3/rhos2
        dtildi(1, 3) = mata(1, 2)/rhost
        dtildi(1, 4) = mata(1, 2)*x3/rhost
        dtildi(1, 5) = mata(1, 3)/rhos2
!
        dtildi(2, 2) = mata(1, 1)*x32/rhos2
        dtildi(2, 3) = mata(1, 2)*x3/rhost
        dtildi(2, 4) = mata(1, 2)*x32/rhost
        dtildi(2, 5) = mata(1, 3)*x3/rhos2
!
        dtildi(3, 3) = mata(2, 2)/rhot2
        dtildi(3, 4) = mata(2, 2)*x3/rhot2
        dtildi(3, 5) = mata(2, 3)/rhost
!
        dtildi(4, 4) = mata(2, 2)*x32/rhot2
        dtildi(4, 5) = mata(2, 3)*x3/rhost
!
        dtildi(5, 5) = mata(3, 3)/rhos2
!
        dtildi(2, 1) = dtildi(1, 2)
        dtildi(3, 1) = dtildi(1, 3)
        dtildi(3, 2) = dtildi(2, 3)
        dtildi(4, 1) = dtildi(1, 4)
        dtildi(4, 2) = dtildi(2, 4)
        dtildi(4, 3) = dtildi(3, 4)
        dtildi(5, 1) = dtildi(1, 5)
        dtildi(5, 2) = dtildi(2, 5)
        dtildi(5, 3) = dtildi(3, 5)
        dtildi(5, 4) = dtildi(4, 5)
!
    else
!
        if (testl1) then
            rhos = 1.d0
        else
            rhos = 1.d0+x3*cour
        end if
!
        x32 = x3*x3
        rhos2 = rhos*rhos
!
        dtildi(1, 1) = mata(1, 1)/rhos2
        dtildi(1, 2) = mata(1, 1)*x3/rhos2
        dtildi(1, 3) = mata(1, 2)/rhos2
        dtildi(2, 2) = mata(1, 1)*x32/rhos2
        dtildi(2, 3) = mata(1, 2)*x3/rhos2
        dtildi(3, 3) = mata(2, 2)/rhos2
!
        dtildi(2, 1) = dtildi(1, 2)
        dtildi(3, 1) = dtildi(1, 3)
        dtildi(3, 2) = dtildi(2, 3)
    end if
!
end subroutine
