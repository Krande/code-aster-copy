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

subroutine lcjpla(fami, kpg, ksp, loi, mod, &
                  nr, imat, nmat, mater, nvi, &
                  deps, sigf, vin, dsde, sigd, &
                  vind, vp, vecp, theta, dt, &
                  devg, devgii, codret)
! aslint: disable=W1504
    implicit none
!       MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT ELASTO-PLASTIQUE
!       VISCO-PLASTIQUE EN VITESSE A T+DT OU T
!       IN  FAMI   :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!           KPG,KSP:  NUMERO DU (SOUS)POINT DE GAUSS
!           LOI    :  MODELE DE COMPORTEMENT
!           MOD    :  TYPE DE MODELISATION
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           NVI    :  NB VARIABLES INTERNES
!           NR     :  NB DE TERME DANS LE SYSTEME DE NEWTOW
!           DEPS   :  INCREMENT DE DEFORMATION
!           SIGF   :  CONTRAINTE A L INSTANT +
!           VIN    :  VARIABLES INTERNES A L INSTANT +
!           SIGD   :  CONTRAINTE A L INSTANT -
!           VIND   :  VARIABLES INTERNES A L INSTANT -
!           THETA  :  ?? COMP_INCR/PARM_THETA
!           DT     :  ??
!           DEVG   :  ??
!           DEVGII :  ??
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!           VP     : VALEURS PROPRES DU DEVIATEUR ELASTIQUE (HOEK-BROWN)
!           VECP   : VECTEURS PROPRES DU DEVIATEUR ELASTIQUE(HOEK-BROWN)
!           CODRET : CODE RETOUR
!                    = 0, TOUT VA BIEN PAS DE REDECOUPAGE
!                    = 1 ou 2, CORRESPOND AU CODE RETOUR DE PLASTI.F
!       ----------------------------------------------------------------
#include "asterfort/hbrjpl.h"
#include "asterfort/irrjpl.h"
#include "asterfort/lgljpl.h"
#include "asterfort/rsljpl.h"
    integer(kind=8) :: imat, nmat, nvi, nr, kpg, ksp, codret
    real(kind=8) :: dsde(6, 6), devg(*), devgii, sigf(6), deps(6)
    real(kind=8) :: vin(*), vind(*), theta, dt, mater(nmat, 2)
    real(kind=8) :: vp(3), vecp(3, 3), sigd(6)
    character(len=8) :: mod
    character(len=16) :: loi
    character(len=*) :: fami
!       ----------------------------------------------------------------
!
    codret = 0
!
    if (loi(1:8) .eq. 'ROUSS_PR' .or. loi(1:10) .eq. 'ROUSS_VISC') then
        call rsljpl(fami, kpg, ksp, loi, imat, &
                    nmat, mater, sigf, vin, vind, &
                    deps, theta, dt, dsde)
!
    else if (loi(1:6) .eq. 'LAIGLE') then
        call lgljpl(mod, nmat, mater, sigf, devg, &
                    devgii, vin, dsde, codret)
!
    elseif ((loi(1:10) .eq. 'HOEK_BROWN') .or. (loi(1:14) .eq. &
                                                'HOEK_BROWN_EFF')) then
        call hbrjpl(mod, nmat, mater, sigf, vin, &
                    vind, vp, vecp, dsde)
!
    else if (loi(1:7) .eq. 'IRRAD3M') then
        call irrjpl(mod, nmat, mater, sigf, vind, &
                    vin, dsde)
    end if
!
end subroutine
