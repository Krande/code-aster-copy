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

subroutine dpmate(mod, imat, materf, ndt, ndi, &
                  nvi, typedp)
    implicit none
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nvi, imat, typedp
    real(kind=8) :: materf(5, 2), ltyped(1)
    character(len=8) :: mod
! ======================================================================
! --- RECUPERATION DES DONNEES MATERIAU POUR LA LOI DE DRUCKER PRAGER --
! ======================================================================
!       IN  IMAT   :  ADRESSE DU MATERIAU CODE
!       OUT MATERF :  COEFFICIENTS MATERIAU A T+DT
!                     MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                     MATER(*,2) = CARACTERISTIQUES   PLASTIQUES
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
! ======================================================================
    real(kind=8) :: trois, deux, un, six, alpha, sy, syult, c, a, phi
    real(kind=8) :: typed, tabtmp(4), coe, dilat(1), psi
    integer(kind=8) :: icodre(8)
    character(len=8) :: nomc(8)
! ======================================================================
    parameter(six=6.0d0)
    parameter(trois=3.0d0)
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
! ======================================================================
! --- DEFINITION PARAMETRES MATERIAU ELAS ------------------------------
! ======================================================================
    nomc(1) = 'E'
    nomc(2) = 'NU'
    nomc(3) = 'ALPHA'
! ======================================================================
! --- RECUPERATION DONNEES MATERIAU ELAS -------------------------------
! ======================================================================
    materf(3, 1) = 0.0d0
    call rcvala(imat, ' ', 'ELAS', 0, ' ', &
                [0.d0], 2, nomc(1), materf(1, 1), icodre, &
                1)
    call rcvala(imat, ' ', 'ELAS', 0, ' ', &
                [0.d0], 1, nomc(3), materf(3, 1), icodre, &
                0)
! ======================================================================
! --- DEFINITION PARAMETRES MATERIAU DRUCKER ---------------------------
! ======================================================================
    nomc(4) = 'ALPHA'
    nomc(5) = 'SY'
    nomc(6) = 'P_ULTM'
! ======================================================================
! --- RECUPERATION MATERIAU SUIVANT LE TYPE D ECROUISSAGE --------------
! ======================================================================
    typed = r8vide()
    call rcvala(imat, ' ', 'DRUCK_PRAGER', 0, ' ', &
                [0.d0], 1, 'TYPE_DP', ltyped(1), icodre, &
                1)
    typed = ltyped(1)
    if (nint(typed) .eq. 1) then
! ======================================================================
! --- CAS LINEAIRE -----------------------------------------------------
! ======================================================================
        typedp = 1
        nomc(7) = 'H'
        call rcvala(imat, ' ', 'DRUCK_PRAGER', 0, ' ', &
                    [0.d0], 4, nomc(4), tabtmp(1), icodre, &
                    1)
! ======================================================================
! --- POUR DES COMMODITES DE PROGRAMMATION ON DEFINIT LES PARAMETRES ---
! --- MATERIAU DE LA FACON SUIVANTE ------------------------------------
! ======================================================================
        materf(1, 2) = tabtmp(2)
        materf(2, 2) = tabtmp(4)
        materf(3, 2) = tabtmp(1)
        materf(4, 2) = tabtmp(3)
        coe = materf(2, 2)+trois*materf(1, 1)*(un/deux/(un+materf(2, 1))+materf(3, 2)*mater&
              &f(3, 2)/(un-deux*materf(2, 1)))
        if (coe .le. 0.0d0) then
            call utmess('F', 'ALGORITH3_37')
        end if
    else if (nint(typed) .eq. 2) then
! ======================================================================
! --- CAS PARABOLIQUE --------------------------------------------------
! ======================================================================
        typedp = 2
        nomc(7) = 'SY_ULTM'
        call rcvala(imat, ' ', 'DRUCK_PRAGER', 0, ' ', &
                    [0.d0], 4, nomc(4), tabtmp(1), icodre, &
                    1)
! ======================================================================
! --- POUR DES COMMODITES DE PROGRAMMATION ON DEFINIT LES PARAMETRES ---
! --- MATERIAU DE LA FACON SUIVANTE ------------------------------------
! ======================================================================
        alpha = tabtmp(1)
        sy = tabtmp(2)
        syult = tabtmp(4)
        phi = atan2((trois*alpha/deux/sqrt((deux*alpha+1.0d0)*(1.0d0-alpha))), 1.0d0)
!           PHI   = ASIN ( TROIS * ALPHA / ( DEUX + ALPHA ) )
        c = (trois-sin(phi))*sy/six/cos(phi)
        a = sqrt(syult/sy)
        materf(1, 2) = a
        materf(2, 2) = phi
        materf(3, 2) = c
        materf(4, 2) = tabtmp(3)
        nomc(8) = 'DILAT'
        call rcvala(imat, ' ', 'DRUCK_PRAGER', 0, ' ', &
                    [0.d0], 1, nomc(8), dilat(1), icodre, &
                    1)
        psi = atan2((trois*dilat(1)/deux/sqrt((deux*dilat(1)+1.0d0)*(1.0d0-dilat(1)))), 1.0d0)
        materf(5, 2) = psi
    else
        ASSERT(.false.)
    end if
! ======================================================================
! --- NOMBRE DE COMPOSANTES --------------------------------------------
! ======================================================================
    if (mod(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
    else if ((mod(1:6) .eq. 'D_PLAN') .or. (mod(1:4) .eq. 'AXIS')) then
        ndt = 4
        ndi = 3
    end if
! ======================================================================
! --- NOMBRE DE VARIABLES INTERNES -------------------------------------
! ======================================================================
    nvi = 3
! ======================================================================
end subroutine
