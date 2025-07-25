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

subroutine vdefro(np, matev, tensel, tenloc)
!.======================================================================
    implicit none
!
!      VDEFRO   -- PASSAGE DU VECTEUR DES EFFORTS GENERALISES
!                  OU DU VECTEUR DES DEFORMATIONS-COURBURES
!                  DU REPERE INTRINSEQUE AUX NOEUDS
!                  OU AUX POINTS D'INTEGRATION DE L'ELEMENT
!                  AU REPERE UTILISATEUR POUR LES ELEMENTS DE
!                  COQUE EPAISSE 3D .
!
!                 CETTE ROUTINE EST ANALOGUE A DXEFRO QUI EST
!                 OPERATIONELLE POUR LES ELEMENTS DE PLAQUE
!                 A L'EXCEPTION DES MATRICES DE PASSAGE QUI
!                 SONT DEFINIES EN DES POINTS DE L'ELEMENT.
!
!   ARGUMENT        E/S   TYPE         ROLE
!    NP             IN     I        NOMBRE DE POINTS OU SONT CALCULES
!                                   LES TENSEURS (I.E. IL S'AGIT DES
!                                   NOEUDS OU DES POINTS D'INTEGRATION
!                                   DE L'ELEMENT)
!    MATEV(2,2,10)  IN     R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX POINTS  DE
!                                   L'ELEMENT AU REPERE UTILISATEUR
!    TENSEL(1)      IN     R        VECTEUR DES EFFORTS GENERALISES
!                                   OU DES DEFORMATIONS-COURBURES
!                                   DANS LE REPERE INTRINSEQUE A
!                                   L'ELEMENT I.E.
!                                       NXX NYY NXY MXX MYY MXY VX VY
!                                   OU  EXX EYY EXY KXX KYY KXY GAX GAY
!    TENLOC(1)      OUT    R        VECTEUR DES EFFORTS GENERALISES
!                                   OU DES DEFORMATIONS-COURBURES
!                                   DANS LE REPERE UTILISATEUR
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/utbtab.h"
    real(kind=8) :: matev(2, 2, 1), tensel(*), tenloc(*)
! -----  VARIABLES LOCALES
    real(kind=8) :: nelem(4), melem(4), xab(2, 2)
    real(kind=8) :: nlocal(4), mlocal(4)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- BOUCLE SUR LES POINTS OU SONT CALCULES LES VECTEURS
! --- (I.E. LES NOEUDS OU LES POINTS D'INTEGRATION) :
!     ============================================
!-----------------------------------------------------------------------
    integer(kind=8) :: i, np
!-----------------------------------------------------------------------
    do i = 1, np
!
        nelem(1) = tensel(1+8*(i-1))
        nelem(2) = tensel(3+8*(i-1))
        nelem(3) = tensel(3+8*(i-1))
        nelem(4) = tensel(2+8*(i-1))
!
        melem(1) = tensel(4+8*(i-1))
        melem(2) = tensel(6+8*(i-1))
        melem(3) = tensel(6+8*(i-1))
        melem(4) = tensel(5+8*(i-1))
!
        call utbtab('ZERO', 2, 2, nelem, matev(1, 1, i), &
                    xab, nlocal)
        call utbtab('ZERO', 2, 2, melem, matev(1, 1, i), &
                    xab, mlocal)
!
        tenloc(1+8*(i-1)) = nlocal(1)
        tenloc(2+8*(i-1)) = nlocal(4)
        tenloc(3+8*(i-1)) = nlocal(2)
!
        tenloc(4+8*(i-1)) = mlocal(1)
        tenloc(5+8*(i-1)) = mlocal(4)
        tenloc(6+8*(i-1)) = mlocal(2)
!
        tenloc(7+8*(i-1)) = tensel( &
                            7+8*(i-1))*matev(1, 1, i)+tensel(8+8*(i-1))*matev(2, 1, i)
        tenloc(8+8*(i-1)) = tensel( &
                            7+8*(i-1))*matev(1, 2, i)+tensel(8+8*(i-1))*matev(2, 2, i)
!
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
