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

subroutine evala1(fami, kpg, ksp, mod, relcom, &
                  sig, vin, imat, module, icode)
! person_in_charge: alexandre.foucault at edf.fr
! =====================================================================
!    - FONCTION REALISEE:  CALCUL DU MODULE DE RIGIDITE
!                          DE MICRO-DILTATION CONTENU
!                          AUX POINTS DE GAUSS
!    - ARGUMENTS:
!        DONNEES:      RELCOM  -->  RELATION DE COMPORTEMENT LOCAL
!                      MOD     -->  TYPE DE MODELISATION
!                      SIG     --> ETAT DE CONTRAINTES
!                      VIN     --> VARIABLES INTERNES
!                      IMAT    --> INDICE DU MATERIAU
!        SORTIE:       MODULE  --> MODULE DE RIGIDITE
!                                  DE MICRO-DILATATION
!                      ICODE   -->
! =====================================================================
    implicit none
#include "asterc/r8prem.h"
#include "asterc/r8pi.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/hujtid.h"
#include "asterfort/redrpr.h"
    real(kind=8) :: module, sig(6), vin(50), dsde(6, 6), sigt(2, 2)
    real(kind=8) :: angle, rot(2, 2), degr, sigr1(2, 2), sigr2(2, 2), temp
    real(kind=8) :: sigf(6), ratio1, pi, value1, ratio2, value2, valeur
    real(kind=8) :: incang(3), angref, valmax
    integer(kind=8) :: i, j, k, iang, nbind(3), disc, kpg, ksp
    integer(kind=8) :: imat, icode, iangmx(3)
    character(len=*) :: fami
    character(len=8) :: mod
    character(len=16) :: relcom
!
! =====================================================================
! DEFINITION DES ELEMENTS NECESSAIRES A L'EVALUATION DU MODULE
! =====================================================================
! --- DISCRETISATION ANGULAIRE DE LA RECHERCHE ENTRE [0° ET 90°[
! =====================================================================
! =====================================================================
! --- DEFINITION DE 3 NIVEAUX DE RECHERCHE
! --- QUI FONT VARIER LA DISCRETISATION ANGULAIRE
! =====================================================================
    pi = r8pi()
    degr = r8dgrd()
    incang(1) = 5.d0
    nbind(1) = 36
    incang(2) = 1.0d0
    nbind(2) = 9
    incang(3) = 2.0d-1
    nbind(3) = 9
    iangmx(1) = 0
    iangmx(2) = 0
    iangmx(3) = 0
! =====================================================================
! --- INITIALISATION MATRICE DE ROTATION ANGULAIRE
! =====================================================================
    do i = 1, 2
        do j = 1, 2
            rot(i, j) = 0.d0
        end do
    end do
! =====================================================================
! --- EXPRESSION DU TENSEUR DES CONTRAINTES SOUS FORME MATRICIEL
! --- LIMITE AU 2D POUR L'INSTANT
! =====================================================================
    sigt(1, 1) = sig(1)
    sigt(1, 2) = sig(4)/sqrt(2.0d0)
    sigt(2, 1) = sigt(1, 2)
    sigt(2, 2) = sig(2)
! =====================================================================
! --- INITIALISATION A ZERO DE LA VALEUR DU MODULE
! --- UNIQUE VARIABLE DE SORTIE
! =====================================================================
    module = 0.d0
    valeur = 0.d0
! =====================================================================
! --- INITIALISATION A ZERO DE L'ANGLE DE DEPART DE RECHERCHE
! =====================================================================
    angref = 0.d0
! =====================================================================
! BOUCLE SUR LA VALEUR DE L'INCREMENT DE DISCRETISATION
! =====================================================================
    do disc = 1, 3
! =====================================================================
! BOUCLE SUR LES ANGLES POUR RECHERCHER LA VALEUR DU MODULE
! =====================================================================
        do iang = 1, nbind(disc)
! =====================================================================
! CONSTRUCTION MATRICE DE ROTATION
! =====================================================================
            angle = iang*incang(disc)*degr+angref
            rot(1, 1) = cos(angle)
            rot(1, 2) = -sin(angle)
            rot(2, 1) = sin(angle)
            rot(2, 2) = cos(angle)
! =====================================================================
! CALCUL DE L'ETAT DE CONTRAINTES SUITE A LA ROTATION
! --- SIGF = TRANSPOSE(ROT)*SIG*ROT
! =====================================================================
            do i = 1, 2
                do j = 1, 2
                    temp = 0.d0
                    do k = 1, 2
                        temp = temp+sigt(i, k)*rot(k, j)
                    end do
                    sigr1(i, j) = temp
                end do
            end do
!
! =====================================================================
! --- TRANSPOSE(ROT) = ROT
! =====================================================================
            temp = rot(2, 1)
            rot(2, 1) = rot(1, 2)
            rot(1, 2) = temp
!
            do i = 1, 2
                do j = 1, 2
                    temp = 0.d0
                    do k = 1, 2
                        temp = temp+rot(i, k)*sigr1(k, j)
                    end do
                    sigr2(i, j) = temp
                end do
            end do
!
! =====================================================================
! --- EXPRESSION DES CONTRAINTES SOUS FORME VECTORIELLE
! =====================================================================
            sigf(1) = sigr2(1, 1)
            sigf(2) = sigr2(2, 2)
            sigf(3) = sig(3)
            sigf(4) = sigr2(1, 2)*sqrt(2.0d0)
            sigf(5) = sig(5)
            sigf(6) = sig(6)
!
! =====================================================================
! --- CALCUL DE LA MATRICE TANGENTE -----------------------------------
! --- (FONCTION DE LA RELATION DE COMPORTEMENT) -----------------------
! =====================================================================
            if (relcom .eq. 'DRUCK_PRAGER') then
! =====================================================================
! --- LOI DE TYPE DRUCKER_PRAGER --------------------------------------
! =====================================================================
                call redrpr(mod, imat, sigf, vin, dsde, &
                            icode)
! =====================================================================
! ----------- LOI DE TYPE HUJEUX --------------------------------------
! =====================================================================
            else if (relcom .eq. 'HUJEUX') then
!
                call hujtid(fami, kpg, ksp, mod, imat, &
                            sigf, vin, dsde, icode)
!
            else
! =====================================================================
!C RELATION DE COMPORTEMENT INVALIDE
! =====================================================================
                ASSERT(.false.)
            end if
! =====================================================================
! ----------- EVALUATION DU MODULE ------------------------------------
! =====================================================================
            if (abs(dsde(4, 4)) .gt. r8prem()) then
                ratio1 = (dsde(4, 1)*dsde(1, 4)-dsde(4, 4)*dsde(1, 1))/(3*dsde(4, 4))
                value1 = (1.0d0/(2.0d0*pi))**2*ratio1
                ratio2 = (dsde(4, 2)*dsde(2, 4)-dsde(4, 4)*dsde(2, 2))/(3*dsde(4, 4))
                value2 = (1.0d0/(2.0d0*pi))**2*ratio2
!
                if (iang .eq. 1) then
                    iangmx(disc) = iang
                    if (disc .eq. 1) then
                        valeur = value1
                    else
                        valmax = max(value1, value2)
                        if (valmax .ge. valeur) then
                            valeur = valmax
                        end if
                    end if
                else
                    valmax = max(value1, value2)
                    if (valmax .ge. valeur) then
                        valeur = valmax
                        iangmx(disc) = iang
                    end if
                end if
            end if
!
        end do
! =====================================================================
! ------ FIN DE LA BOUCLE POUR CETTE VALEUR DE DISCRETISATION ---------
! =====================================================================
! =====================================================================
! ------ CHANGEMENT DE L'ANGLE DE REFERENCE POUR DEBUTER
! ------ LA RECHERCHE DU MAXIMUM
! =====================================================================
        angref = (iangmx(disc)-1)*incang(disc)*degr+angref
    end do
! =====================================================================
! ------ FIN DE LA BOUCLE POUR OBTENIR LE MAX DE VALUE1 ---------
! =====================================================================
!
    module = max(module, valeur)
!
end subroutine
