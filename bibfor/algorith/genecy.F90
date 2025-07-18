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
subroutine genecy(cmod1, cmod2, neq, lmat, para, &
                  nbsec, beta1, beta2, ctrav)
!    C. VARE     DATE 20/01/94
!-----------------------------------------------------------------------
!  BUT: CALCULER LES PARAMETRES GENERALISES DES MODES CALCULES
    implicit none
!       PAR UNE METHODE CYCLIQUE
!
!       LES MODES ETANT A PRIORI DOUBLES, IL Y A DEUX PARAMETRES
!       EN SORTIE
!-----------------------------------------------------------------------
!
! CMOD1    /I/: VECTEUR DU PREMIER MODE COMPLEXE
! CMOD2    /I/: VECTEUR DU DEUXIEME MODE COMPLEXE
! NEQ      /I/: NOMBRE D'EQUATIONS ASSEMBLEES
! LMAT     /I/: ADRESSE DESCRIPTEUR MATRICE
! PARA     /O/: VECTEUR DES DEUX PARAMETRES GENERALISES
! NBSEC    /I/: NOMBRE DE SECTEURS
! BETA1    /I/: DEPHASAGE INTER-SECTEUR DU PREMIER MODE
! BETA2    /I/: DEPHASAGE INTER-SECTEUR DU DEUXIEME MODE
! CTRAV    /M/: VECTEUR DE TRAVAIL (NEQ)
!
!-----------------------------------------------------------------------
#include "asterfort/mcmult.h"
    integer(kind=8) :: i, j, lmat, nbsec, neq
    real(kind=8) :: beta1, beta2, xima, xrea
    real(kind=8) :: para(2), zero
    complex(kind=8) :: cmod1(neq), cmod2(neq), ctrav(neq), cfact1, cfact2
!-----------------------------------------------------------------------
    data zero/0.d+00/
!-----------------------------------------------------------------------
!
    para(1) = zero
    para(2) = zero
    do i = 1, neq
        ctrav(i) = dcmplx(0.d0, 0.d0)
    end do
!
!------CALCUL DU PRODUIT MATRICE ASSEMBLEE REELLE-MODE COMPLEXE---------
!
    call mcmult('ZERO', lmat, cmod2, ctrav, 1, &
                .true._1)
!
!-------------------BOUCLE SUR LES SECTEURS-----------------------------
!
    do i = 1, nbsec
!
!  CALCUL DU DEPHASAGE DU SECTEUR COURANT (ET DU CONJUGUE)
!
        cfact1 = dcmplx(cos((i-1)*beta1), sin((i-1)*beta1))
        cfact2 = dcmplx(cos((i-1)*beta2), sin((i-1)*beta2))
!
        xrea = zero
        xima = zero
!
!  BOUCLE SUR LES DDL ASSEMBLES POUR PRODUITS SCALAIRES
!
        do j = 1, neq
            xrea = xrea+dble(cfact1*cmod1(j))*dble(cfact2*ctrav(j))
            xima = xima+dimag(cfact1*cmod1(j))*dimag(cfact2*ctrav(j))
        end do
!
        para(1) = para(1)+xrea
        para(2) = para(2)+xima
!
    end do
!
end subroutine
