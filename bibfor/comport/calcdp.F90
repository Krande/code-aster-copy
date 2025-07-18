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
subroutine calcdp(crit, seuil, dt, rprim, mutrbe, &
                  sigm0, epsi0, coefm, dp, iret)
    implicit none
#include "asterfort/calcfp.h"
#include "asterfort/zeroco.h"
    real(kind=8) :: seuil, dt, rprim, mutrbe, crit(3)
    real(kind=8) :: sigm0, epsi0, coefm, dp
!-------------------------------------------------------------
! CALCUL DE DP PAR RESOLUTION DE L'EQUATION SCALAIRE FPRIM=0
! FPRIM=F1-F(I) TEL QUE F1 : DROITE DECROISSANTE
!                       F(I): SIGM0*ARGSH((DP/DT/EPSI0)^(1/M))
! ------------------------------------------------------------
! IN CRIT  : CRITERES DE CONVERGENCE LOCAUX
! IN SEUIL : VALEUR DE LA FONCTION SEUIL F
! IN DT    : PAS DE TEMPS
! IN RPRIM : ECROUISSAGE
! IN MUTRBE : MU * TR (Bel)
! IN SIGM0,EPSI0,COEFM : PARAMETRE LOI VISQUEUSE ROUSS_VISC
! OUT DP    : INCREMENT VISCOPLASTIQUE
! OUT IRET    : CODE RETOUR DE LA DETERMINATION DE DP
!               IRET=0 => PAS DE PROBLEME
!               IRET=1 => ABSENCE DE CONVERGENCE
!-------------------------------------------------------------
    real(kind=8) :: fprim, fplas, dfprim, dpnew, fpnew, test, arg
    real(kind=8) :: r0, dpplas, dpmin, dpmax, x(4), y(4)
    integer(kind=8) :: iter, iret
!
! INITIALISATION
    r0 = mutrbe+rprim
!
! ESSAI PLASTIQUE PUR (DP = DPMAX)
!
    dpplas = seuil/r0
    dp = dpplas
!    CALCUL DU SEUIL
    call calcfp(mutrbe, rprim, seuil, dt, dp, &
                sigm0, epsi0, coefm, fplas, fprim, &
                dfprim)
    if (abs(fprim)/(1.d0+seuil) .lt. crit(3)) goto 999
!
    dpmax = dpplas
    x(1) = dp
    y(1) = fprim
!
!  RESOLUTION DE FPRIM=0 - METHODE SECANTE AVEC NEWTON SI BESOIN
!
!------EXAMEN DE LA SOLUTION DP=DPMIN
    dpmin = 0.d0
    dp = dpmin
    call calcfp(mutrbe, rprim, seuil, dt, dp, &
                sigm0, epsi0, coefm, fplas, fprim, &
                dfprim)
    x(2) = dp
    y(2) = fprim
!
!     AMELIORATION DES BORNES PAR ESTMATION DE DP^VP
    arg = (dpmax/dt/epsi0)**(1.d0/coefm)
    test = sigm0*log(arg+sqrt(arg**2+1))
!
    if (fplas .lt. test) then
        dp = dt*epsi0*(sinh(fplas/sigm0))**coefm
        if ((dp .ge. dpmin) .and. (dp .le. dpmax)) then
            call calcfp(mutrbe, rprim, seuil, dt, dp, &
                        sigm0, epsi0, coefm, fplas, fprim, &
                        dfprim)
            if (abs(fprim)/(1.d0+seuil) .lt. crit(3)) goto 999
            if (fprim .lt. 0) then
                x(1) = dp
                y(1) = fprim
            else
                x(2) = dp
                y(2) = fprim
            end if
        end if
    end if
!
!------CALCUL DE DP : EQUATION SCALAIRE FPRIM = 0 AVEC DPMIN< DP < DPMAX
    x(3) = x(1)
    y(3) = y(1)
    x(4) = x(2)
    y(4) = y(2)
    do iter = 1, int(crit(1))
!
        call zeroco(x, y)
        dp = x(4)
        call calcfp(mutrbe, rprim, seuil, dt, dp, &
                    sigm0, epsi0, coefm, fplas, fprim, &
                    dfprim)
        if (abs(fprim)/(1.d0+seuil) .lt. crit(3)) goto 999
!
        dpnew = dp-fprim/dfprim
        if ((dpnew .ge. dpmin) .and. (dpnew .le. dpmax)) then
            call calcfp(mutrbe, rprim, seuil, dt, dpnew, &
                        sigm0, epsi0, coefm, fplas, fpnew, &
                        dfprim)
!
            if (abs(fpnew)/(1.d0+seuil) .lt. crit(3)) then
                dp = dpnew
                goto 999
            end if
            if (abs(fpnew)/(1.d0+seuil) .lt. abs(fprim)/(1.d0+seuil)) then
                dp = dpnew
                fprim = fpnew
            end if
!
        end if
        y(4) = fprim
        x(4) = dp
    end do
!
    iret = 1
!
999 continue
!
end subroutine
