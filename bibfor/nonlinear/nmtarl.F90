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
subroutine nmtarl(mode, ndimsi, mat, sigel, vim, &
                  epm, dp, sp, xi, dirdp, &
                  dirsp, dirxi, min, rho, ener)
!
    implicit none
#include "asterfort/nmtacr.h"
#include "asterfort/zeroco.h"
    integer(kind=8) :: mode, ndimsi
    real(kind=8) :: mat(*), sigel(*), vim(*), epm(*), dp, sp, xi
    real(kind=8) :: dirdp, dirsp, dirxi, min
    real(kind=8) :: rho, ener
! ----------------------------------------------------------------------
! TAHERI: MINIMISATION DE F**2+G**2 RT A (P,SP,XI)+RHO(DIRP,DIRSP,DIRXI)
! ----------------------------------------------------------------------
! IN  MODE   2: (P,SP)     3: (XI,SP)
! IN  NDIMSI DIMENSION DES TENSEURS
! IN  MAT    TABLEAU DES CONSTANTES MATERIAUX
! IN  SIGEL  DEVIATEUR DES CONTRAINTES ELASTIQUES
! IN  VIM    VARIABLES INTERNES EN T-
! IN  EPM    DEFORMATION PLASTIQUE EN T-
! IN  DP     INCREMENT DE DEFORMATION PLASTIQUE CUMULEE
! IN  SP     CONTRAINTE DE PIC
! IN  XI     PILOTAGE DE EPN
! IN  DIRDP  DIRECTION POUR DP
! IN  DIRSP  DIRECTION POUR SP
! IN  DIRXI  DIRECTION POUR XI
! IN  MIN    VALEUR DE LA DERIVEE EN RHO=0
! VAR RHO    DISTANCE PARCOURUE  IN: RHOMAX  OUT: RHO
! VAR ENER   VALEUR DE (F**2+G**2)/2   IN: EN RHO=0   OUT: EN RHO
! ----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: niter, itelin
    real(kind=8) :: f, g, fds, gds, fdp, gdp, fdx, gdx, dpmax, sig(6)
    real(kind=8) :: tang(6, 6)
    real(kind=8) :: x(4), y(4), energ(4)
    real(kind=8) :: rhomax, refe, prelin
!
    parameter(prelin=1.d-2, itelin=3)
!
!
!    INITIALISATION DE LA SOLUTION RHO = 0
    x(1) = 0.d0
    y(1) = min
    refe = min
    if (refe .ge. 0.d0) then
        rho = 0.d0
        goto 999
    end if
    energ(1) = ener
!
!
!    EXAMEN DE LA SOLUTION RHO = RHOMAX
    rhomax = rho
    x(2) = rhomax
    call nmtacr(mode, ndimsi, mat, sigel, vim, &
                epm, dp+rhomax*dirdp, sp+rhomax*dirsp, xi+rhomax*dirxi, f, &
                g, fds, gds, fdp, gdp, &
                fdx, gdx, dpmax, sig, tang)
    if (mode .eq. 2) then
        y(2) = (f*fdp+g*gdp)*dirdp+(f*fds+g*gds)*dirsp
    else
        y(2) = (f*fdx+g*gdx)*dirxi+(f*fds+g*gds)*dirsp
    end if
    ener = (f**2+g**2)/2.d0
    if (y(2) .le. 0.d0) then
        rho = rhomax
        goto 999
    end if
    energ(2) = ener
!
!
!    CALCUL DE RHO
!
    x(3) = x(1)
    y(3) = y(1)
    energ(3) = energ(1)
    x(4) = x(2)
    y(4) = y(2)
    energ(4) = energ(2)
!
    do niter = 1, itelin
        if (abs(y(4)/refe) .lt. prelin) goto 110
        call zeroco(x, y)
        call nmtacr(mode, ndimsi, mat, sigel, vim, &
                    epm, dp+x(4)*dirdp, sp+x(4)*dirsp, xi+x(4)*dirxi, f, &
                    g, fds, gds, fdp, gdp, &
                    fdx, gdx, dpmax, sig, tang)
        if (mode .eq. 2) then
            y(4) = (f*fdp+g*gdp)*dirdp+(f*fds+g*gds)*dirsp
        else
            y(4) = (f*fdx+g*gdx)*dirxi+(f*fds+g*gds)*dirsp
        end if
        energ(4) = (f**2+g**2)/2.d0
    end do
110 continue
    rho = x(4)
    ener = energ(4)
!
!
999 continue
end subroutine
