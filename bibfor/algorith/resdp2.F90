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

subroutine resdp2(materf, seq, i1e, pmoins, dp, &
                  plas)
    implicit none
#include "asterfort/schdp2.h"
#include "asterfort/utmess.h"
    real(kind=8) :: materf(5, 2), pmoins, dp, seq, i1e, plas
! =====================================================================
! --- RESOLUTION NUMERIQUE --------------------------------------------
! =====================================================================
    integer(kind=8) :: ndt, ndi
    real(kind=8) :: young, nu, troisk, deuxmu, alpha1, phi, c, pult, alpha
    real(kind=8) :: trois, deux, un, fcrit, valpro, gamapm, gamarp
    real(kind=8) :: neuf, douze, a1, b1, delta, quatre, valcoe, b2
    real(kind=8) :: fcrit0, pptest
    parameter(douze=12.0d0)
    parameter(neuf=9.0d0)
    parameter(quatre=4.0d0)
    parameter(trois=3.0d0)
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
! =====================================================================
    common/tdim/ndt, ndi
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    young = materf(1, 1)
    nu = materf(2, 1)
    troisk = young/(un-deux*nu)
    deuxmu = young/(un+nu)
    alpha1 = materf(1, 2)
    phi = materf(2, 2)
    c = materf(3, 2)
    pult = materf(4, 2)
    gamarp = sqrt(trois/deux)*pult
    gamapm = sqrt(trois/deux)*pmoins
    alpha = deux*sin(phi)/(trois-sin(phi))
! =====================================================================
! --- CALCUL ELASTIQUE ------------------------------------------------
! =====================================================================
    fcrit = schdp2(seq, i1e, phi, alpha1, c, pult, pmoins)
! =====================================================================
! --- CALCUL PLASTIQUE ------------------------------------------------
! =====================================================================
    if (fcrit .gt. 0.0d0) then
        plas = 1.0d0
        if (pmoins .lt. pult) then
            a1 = -neuf*c*cos(phi)*(un-alpha1)*(un-alpha1)/gamarp/gamarp/(trois-sin(phi))
            b1 = -( &
                 trois*deuxmu/deux+trois*troisk*alpha*alpha-sqrt(trois/deux)*douze*c*cos(phi)&
                 &/(trois-sin(phi))*(un-(un-alpha1)/gamarp*gamapm)*(un-alpha1)/gamarp &
                 )
            delta = b1*b1-quatre*a1*fcrit
            if (a1 .eq. 0.0d0) then
                call utmess('F', 'ALGORITH10_43')
            end if
            dp = -(b1+sqrt(delta))/deux/a1
            valcoe = sqrt(deux/trois)*(gamarp-gamapm)
            if (dp .gt. valcoe) then
                fcrit = schdp2(seq, i1e, phi, alpha1, c, pult, pult)
                b2 = -(trois*deuxmu/deux+trois*troisk*alpha*alpha)
                if (b2 .eq. 0.0d0) then
                    call utmess('F', 'ALGORITH10_42')
                end if
                dp = -fcrit/b2
            end if
        else
            b2 = -(trois*deuxmu/deux+trois*troisk*alpha*alpha)
            if (b2 .eq. 0.0d0) then
                call utmess('F', 'ALGORITH10_42')
            end if
            dp = -fcrit/b2
        end if
    else
        plas = 0.0d0
        dp = 0.0d0
    end if
! =====================================================================
! --- PROJECTION AU SOMMET --------------------------------------------
! =====================================================================
    pptest = pmoins+dp
    b2 = trois*troisk*alpha*alpha
    fcrit0 = schdp2(0.0d0, i1e, phi, alpha1, c, pult, pptest)
    valpro = fcrit0/b2
!
    a1 = -neuf*c*cos(phi)*(un-alpha1)*(un-alpha1)/gamarp/gamarp/(trois-sin(phi))
!
    if ((plas .eq. 1) .and. (dp .le. valpro)) then
        plas = 2.0d0
        fcrit = schdp2(0.0d0, i1e, phi, alpha1, c, pult, pmoins)
        if (pmoins .lt. pult) then
            b1 = -( &
                 trois*troisk*alpha*alpha-sqrt(trois/deux)*douze*c*cos(phi)/(trois-sin(phi))* &
                 &(un-(un-alpha1)/gamarp*gamapm)*(un-alpha1)/gamarp &
                 )
            delta = b1*b1-quatre*a1*fcrit
            dp = -(b1+sqrt(delta))/deux/a1
            valcoe = sqrt(deux/trois)*(gamarp-gamapm)
            if (dp .gt. valcoe) then
                fcrit = schdp2(0.0d0, i1e, phi, alpha1, c, pult, pult)
                dp = fcrit/b2
            end if
        else
            dp = fcrit/b2
        end if
    end if
! =====================================================================
end subroutine
