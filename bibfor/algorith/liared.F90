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
subroutine liared(nomres, fmli, iblo, liamod, nlilia, &
                  ncolia, promod, nlipro, ncopro, taille, &
                  indcol, nbcol)
    implicit none
!  O. NICOLAS     DATE 01/08/04
!
!  REFAIT INTEGRALEMENT
!
!    M. CORUS     DATE 05/02/10
!-----------------------------------------------------------------------
!  BUT:      < CALCUL DE LA MATRICE DE LIAISON REDUITE >
!
!  ON PROJETTE LA RESTICTION DE LA BASE MODALE A L'INTERFACE SUR
!  LA RESTRICTION DE LA BASE MODALE ESCLAVE. L'ORIENTATION DES
!  SOUS-STRUCTURES A DEJA ETE PRISE EN COMPTE
!  ON ELIMINE LES LIGNES CORRESPONDANT AUX VECTEURS DU PROJECTEUR
!  NE FAISANT PAS BOUGER L'INTERFACE
!
!-----------------------------------------------------------------------
!
! NOMRES  /I/ : NOM UTILISATEUR DU RESULTAT MODELE_GENE
! FMLI    /I/ : FAMILLE DES MATRICES DE LIAISON
! IBLO    /I/ : NUMERO DU BLOC DE LA LIAISON
! LIAMOD  /I/ : MATRICE A PROJETER
! NLILIA /I/ : NB DE LIGNES DE LA MATRICE A PROJETER
! NCOLIA /I/ : NB DE COLONNES DE LA MATRICE A PROJETER
! PROMOD  /I/ : MATRICE DU PROJECTEUR
! NLIPRO /I/ : NB DE LIGNES DU PROJECTEUR
! NCOPRO /I/ : NB DE COLONNES DU PROJECTEUR
! TAILLE  /O/ : TAILLE DE LA MATRICE DE LAISON
! INDCOL /I-O/ : VECTEUR DES INDICES DES COLONNES DE LIAISONS ACTIVES
! NBCOL   /I-O/ : NOMBRE DE COLONNES ACTIVES
!
!
!
!
!
!-- VARIABLES EN ENTREES / SORTIE
#include "jeveux.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iblo, nlilia, ncolia, nlipro, ncopro, taille(2), nbcol
    character(len=8) :: nomres
    character(len=24) :: fmli, liamod, promod, indcol
!
!-- VARIABLES DE LA ROUTINE
    integer(kind=8) :: i1, j1, k1, l1, lliamo, lpromo, ltemp, lpro
    real(kind=8) :: temp, eps, coeff
    parameter(eps=2.3d-16)
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
!-- TEST SUR LA TAILLE DES MATRICES POUR S'ASSURER DE POUVOIR
!-- FAIRE LA PROJECTION
!
    if (nlipro .eq. nlilia) then
!
        call jeveuo(liamod, 'L', lliamo)
        call jeveuo(promod, 'L', lpromo)
!
        if (liamod .eq. promod) then
            coeff = -1.d0
        else
            coeff = 1.d0
        end if
!
!-- TEST SUR LA PRESENCE DE COLONNES DE ZEROS DANS LE PROJECTEUR
!-- ET DETERMINATION DE LA TAILLE DE LA MATRICE DE LIAISON SI CA N'A
!-- PAS DEJA ETE FAIT
!
        if (indcol(1:5) .eq. 'BLANC') then
            indcol = '&&INDCOLONNES_NON_ZERO'
            call wkvect(indcol, 'V V I', ncopro, lpro)
            do i1 = 1, ncopro
                zi(lpro+i1-1) = 0
            end do
            nbcol = 0
            do j1 = 1, ncopro
                temp = 0.d0
                do i1 = 1, nlipro
                    temp = temp+zr(lpromo+(j1-1)*nlipro+i1-1)**2
                end do
                if (sqrt(temp)/nlipro .gt. eps) then
                    zi(lpro+nbcol) = j1
                    nbcol = nbcol+1
                end if
            end do
        else
            call jeveuo(indcol, 'L', lpro)
        end if
!
! --- CREATION DE LA NOUVELLE MATRICE DE LIAISON
!
        call jecroc(jexnum(fmli, iblo))
        call jeecra(jexnum(fmli, iblo), 'LONMAX', nbcol*ncolia)
        call jeveuo(jexnum(fmli, iblo), 'E', ltemp)
!
! --- INITIALISATION DES MATRICES ORIENTEES DE LIAISON----------
!
        do i1 = 1, nbcol*ncolia
            zr(ltemp+i1-1) = 0.d0
        end do
!
!
!        CALL DGEMM('T','N',NCOPRO,NCOLIA,NLIPRO,1.,ZR(LPROMO),
!     &            NLIPRO,ZR(LLIAMO),NLILIA,0.,ZR(LTEMP),NCOPRO)
!
!-- BOUCLE "A LA MAIN", TEL QUE DANS DGEMM :
!
!   C=A^T*B  => Cij = Sum(Aki*Bkj)
!
!    A : NLA*NCA
!    B : NLA*NCB
!    C : NCA*NCB
!
!   C(i,j) = C( (j-1)*NCA + i-1 )
!   A(i,k) = A( (k-1)*NLA + i-1 )
!    => A(k,i) = A( (k-1)*NLA + i-1 )
!   B(k,j) = B( (j-1)*NLB + k-1 )
!
!   pour j=1,NCB
!     pour k=1,NLA
!       pour i=1,NCA
!         C( (j-1)*NCA+(i-1) ) =
!              C( (j-1)*NCA+ i-1 ) +
!              A( (i-1)*NLA + k-1 ) *
!              B( (j-1)*NLB + k-1 )
!
!
!
!--
!-- PROJECTION INITIALE
!--
!
        do j1 = 1, ncolia
            do k1 = 1, nlipro
                temp = zr(lliamo+(j1-1)*nlilia+k1-1)
                do l1 = 1, nbcol
                    i1 = zi(lpro+l1-1)
                    zr(ltemp+(j1-1)*nbcol+l1-1) = zr(ltemp+(j1-1)* &
                                                   nbcol+l1-1)+coeff*temp*zr(lpromo+(i1-1)*nlipro+ &
                                                                               k1-1)
                end do
            end do
        end do
!
    end if
!
    taille(1) = nbcol
    taille(2) = ncolia
!
    call jedema()
!
end subroutine
