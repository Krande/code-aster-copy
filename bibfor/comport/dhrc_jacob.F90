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

subroutine dhrc_jacob(eps, vint, c, bp1, &
                      cp1, bp2, cp2, as1, bs1, &
                      cs1, as2, bs2, cs2, indi, &
                      neta1, neta2, cstseu, jacob)
!
! person_in_charge: sebastien.fayolle at edf.fr
!
    implicit none
!
    integer(kind=8), intent(in) :: indi(6)
    real(kind=8), intent(in) :: vint(*), eps(8)
    real(kind=8), intent(in) :: as1(6, 6), as2(6, 6)
    real(kind=8), intent(in) :: bp1(6, 2), bp2(6, 2), bs1(6, 2), bs2(6, 2)
    real(kind=8), intent(in) :: c(2, 2, 2), cp1(2, 2), cp2(2, 2), cs1(2, 2), cs2(2, 2)
    real(kind=8), intent(in) :: neta1(2), neta2(2), cstseu(6)
    real(kind=8), intent(out) :: jacob(6, 6)
!
! ----------------------------------------------------------------------
!
!      CALCUL DE LA MATRICE JACOBIENNE DES SEUILS
!      APPELE PAR "DHRC_LC"
!
! IN:
!       EPS     : TENSEUR DE DEFORMATIONS
!                 (EXX EYY 2EXY KXX KYY 2KXY)
!       VINT    : VECTEUR DES VARIABLES INTERNES
!                 VINT=(D1,D2,EPSP1X,EPSP1Y,EPSP2X,EPSP2Y)
!       AS1     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR ELASTIQUE PAR
!                 RAPPORT A D1
!       AS2     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR ELASTIQUE PAR
!                 RAPPORT A D2
!       BP1     : DERIVEE PREMIERE DU TENSEUR DE RAIDEUR COUPLE PAR
!                 RAPPORT A D1
!                 LA PREMIERE COMPOSANTE DE BP CORRESPOND AUX DEFORMATIONS M-F
!                 LA DEUXIEME COMPOSANTE DE BP CORRESPOND AUX GLISSEMENTS
!       BP2     : DERIVEE PREMIERE DU TENSEUR DE RAIDEUR COUPLE PAR
!                 RAPPORT A D2
!                 LA PREMIERE COMPOSANTE DE BP CORRESPOND AUX DEFORMATIONS M-F
!                 LA DEUXIEME COMPOSANTE DE BP CORRESPOND AUX GLISSEMENTS
!       BS1     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR COUPLE PAR
!                 RAPPORT A D1
!                 LA PREMIERE COMPOSANTE DE BS CORRESPOND AUX DEFORMATIONS M-F
!                 LA DEUXIEME COMPOSANTE DE BS CORRESPOND AUX GLISSEMENTS
!       BS2     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR COUPLE PAR
!                 RAPPORT A D2
!                 LA PREMIERE COMPOSANTE DE BS CORRESPOND AUX DEFORMATIONS M-F
!                 LA DEUXIEME COMPOSANTE DE BS CORRESPOND AUX GLISSEMENTS
!       C       : TENSEUR DE RAIDEUR D'ECROUISSAGE
!                 LA PREMIERE COMPOSANTE DE C CORRESPOND AUX GLISSEMENTS
!                 LA DEUXIEME COMPOSANTE DE C CORRESPOND AUX GLISSEMENTS
!                 LA TROISIEME COMPOSANTE DE C CORRESPOND A LA DISTINCTION
!                 ENTRE PARTIE SUPERIEURE ET INFERIEURE DE LA PLAQUE
!       CP1     : DERIVEE PREMIERE DU TENSEUR DE RAIDEUR ECROUISSAGE PAR
!                 RAPPORT A D1
!                 LES COMPOSANTES DE CP CORRESPONDENT AUX GLISSEMENTS
!       CP2     : DERIVEE PREMIERE DU TENSEUR DE RAIDEUR ECROUISSAGE PAR
!                 RAPPORT A D2
!                 LES COMPOSANTES DE CP CORRESPONDENT AUX GLISSEMENTS
!       CS1     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR ECROUISSAGE PAR
!                 RAPPORT A D1
!                 LES COMPOSANTES DE CS CORRESPONDENT AUX GLISSEMENTS
!       CS2     : DERIVEE SECONDE DU TENSEUR DE RAIDEUR ECROUISSAGE PAR
!                 RAPPORT A D2
!       NETA1   : CONTRAINTE ASSOCIEE AU GLISSEMENT SUR LA PARTIE 1 DE
!                 LA PLAQUE
!       NETA2   : CONTRAINTE ASSOCIEE AU GLISSEMENT SUR LA PARTIE 2 DE
!                 LA PLAQUE
!       CSTSEU  : PARAMETRES DE SEUILS
!            (1): POUR L'ENDOMMAGEMENT
!            (2): POUR LE GLISSEMENT
!                 LES COMPOSANTES DE CS CORRESPONDENT AUX GLISSEMENTS
!
! OUT:
!       JACOB   : MATRICE JACOBIENNE DES SEUILS ACTIVÉS
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l
    real(kind=8) :: jacobt(6, 6)
!
! ----------------------------------------------------------------------
! --  CALCUL DE LA JACOBIENNE TOTALE JACOBT
! ----------------------------------------------------------------------
!
    jacobt(:, :) = 0.0d0
!
    do k = 1, 2
        do i = 1, 2
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*as1(k, i)*eps(i)*0.5d0
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*as2(k, i)*eps(i)*0.5d0
!
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*bs1(k, i)*vint(i+2)
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*bs2(k, i)*vint(i+4)
!
            jacobt(1, 1) = jacobt(1, 1)-vint(k+2)*cs1(k, i)*vint(i+2)*0.5d0
            jacobt(2, 2) = jacobt(2, 2)-vint(k+4)*cs2(k, i)*vint(i+4)*0.5d0
        end do
        do i = 3, 6
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*as1(k, i)*eps(i)*0.5d0
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*as2(k, i)*eps(i)*0.5d0
        end do
    end do
    do k = 3, 6
        do i = 1, 2
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*as1(k, i)*eps(i)*0.5d0
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*as2(k, i)*eps(i)*0.5d0
!
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*bs1(k, i)*vint(i+2)
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*bs2(k, i)*vint(i+4)
        end do
        do i = 3, 6
            jacobt(1, 1) = jacobt(1, 1)-eps(k)*as1(k, i)*eps(i)*0.5d0
            jacobt(2, 2) = jacobt(2, 2)-eps(k)*as2(k, i)*eps(i)*0.5d0
        end do
    end do
    do k = 1, 6
!
!       jacobt(1,2)=0.d0 par construction du modèle
!
        jacobt(1, 3) = jacobt(1, 3)-eps(k)*bp1(k, 1)
        jacobt(1, 4) = jacobt(1, 4)-eps(k)*bp1(k, 2)
! Pour évolution
!        jacobt(1,5)=0.d0 - couplage endommagement haut glissement bas x
!        jacobt(1,6)=0.d0 - couplage endommagement haut glissement bas y
!
! Pour évolution
!        jacobt(2,3)=0.d0 - couplage endommagement bas glissement haut x
!        jacobt(2,4)=0.d0 - couplage endommagement bas glissement haut y
!
        jacobt(2, 5) = jacobt(2, 5)-eps(k)*bp2(k, 1)
        jacobt(2, 6) = jacobt(2, 6)-eps(k)*bp2(k, 2)
!
        jacobt(3, 1) = jacobt(3, 1)-eps(k)*bp1(k, 1)
        jacobt(4, 1) = jacobt(4, 1)-eps(k)*bp1(k, 2)
!
! Pour évolution
!        jacobt(5,1)=0.d0 - couplage endommagement haut glissement bas x
!        jacobt(6,1)=0.d0 - couplage endommagement haut glissement bas y
!
! Pour évolution
!        jacobt(3,2)=0.d0 - couplage endommagement bas glissement haut x
!        jacobt(4,2)=0.d0 - couplage endommagement bas glissement haut y
!
        jacobt(5, 2) = jacobt(5, 2)-eps(k)*bp2(k, 1)
        jacobt(6, 2) = jacobt(6, 2)-eps(k)*bp2(k, 2)
!
    end do
    do k = 1, 2
        jacobt(1, 3) = jacobt(1, 3)-vint(k+2)*cp1(k, 1)
        jacobt(1, 4) = jacobt(1, 4)-vint(k+2)*cp1(k, 2)
!
        jacobt(2, 5) = jacobt(2, 5)-vint(k+4)*cp2(k, 1)
        jacobt(2, 6) = jacobt(2, 6)-vint(k+4)*cp2(k, 2)
!
        jacobt(3, 1) = jacobt(3, 1)-vint(k+2)*cp1(k, 1)
!
        jacobt(4, 1) = jacobt(4, 1)-vint(k+2)*cp1(k, 2)
!
        jacobt(5, 2) = jacobt(5, 2)-vint(k+4)*cp2(k, 1)
!
        jacobt(6, 2) = jacobt(6, 2)-vint(k+4)*cp2(k, 2)

    end do
!
    jacobt(3, 3) = jacobt(3, 3)-c(1, 1, 1)
    jacobt(3, 4) = jacobt(3, 4)-c(2, 1, 1)
!
    jacobt(4, 3) = jacobt(4, 3)-c(1, 2, 1)
    jacobt(4, 4) = jacobt(4, 4)-c(2, 2, 1)
!
    jacobt(5, 5) = jacobt(5, 5)-c(1, 1, 2)
    jacobt(5, 6) = jacobt(5, 6)-c(2, 1, 2)
!
    jacobt(6, 5) = jacobt(6, 5)-c(1, 2, 2)
    jacobt(6, 6) = jacobt(6, 6)-c(2, 2, 2)
!
    do k = 1, 6
        jacobt(1, k) = jacobt(1, k)/cstseu(1)
        jacobt(2, k) = jacobt(2, k)/cstseu(2)
        jacobt(3, k) = jacobt(3, k)*2.0d0*neta1(1)/(cstseu(3)**2.0d0)
        jacobt(4, k) = jacobt(4, k)*2.0d0*neta1(2)/(cstseu(4)**2.0d0)
        jacobt(5, k) = jacobt(5, k)*2.0d0*neta2(1)/(cstseu(5)**2.0d0)
        jacobt(6, k) = jacobt(6, k)*2.0d0*neta2(2)/(cstseu(6)**2.0d0)
    end do
!
! ----------------------------------------------------------------------
! --  CALCUL DE LA JACOBIENNE ACTIVEE JACOB
! ----------------------------------------------------------------------
!
    l = 0
    do k = 1, 6
        if (k .eq. indi(k)) then
            l = l+1
            j = 0
            do i = 1, 6
                if (i .eq. indi(i)) then
                    j = j+1
                    jacob(l, j) = jacobt(k, i)
                end if
            end do
        end if
    end do
!
end subroutine
