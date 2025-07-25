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
subroutine pcfull(n, icpl, icpc, icpd, icplp, &
                  icpcp, ind, lca, ier)
!                    S.P. PCFULL
!                    -----------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUT : CE SP CALCULE LES POINTEURS ICPLP,ICPCP
! ----  CORRESPONDANTS AU 'REMPLISSAGE' DE LA
!       MATRICE EN COURS DE FACTORISATION
!
!       VERSION GENERALE : MATRICE A STRUCTURE QUELCONQUE
!
!   PARAMETRES D'ENTREE:
!   -------------------
!
!   ICPL,ICPD,ICPC : LES POINTEURS ASSOCIES A LA MATRICE A FACTORISER
!
!   * ICPL(I) = ADRESSE DANS LE RANGEMENT DES COEFFICIENTS A(I,J)
!               DU DERNIER COEFFICIENT DE LA LIGNE I
!   * ICPC(K) = POUR K = ICPL(I-1)+1...ICPL(I) NUMEROS DES INDICES DE
!               COLONNE J, DES COEFFICIENTS A(I,J) DE LA LIGNE I
!               ( RANGES PAR ORDRE DE J CROISSANT)
!   * ICPD(I) = ADRESSE DANS LE RANGEMENT DES COEFFICIENTS A(I,J)
!               DU DERNIER COEFFICIENT DE LA LIGNE I AVEC A(I,J), J < I
!
!   IND         : TABLEAU UTILITAIRE
!   LCA         : TAILLE MAXIMALE ADMISE POUR LA MATRICE FACTORISEE
!
!   PARAMETRE DE SORTIE:
!   -------------------
!
!   ICPLP,ICPCP : LES POINTEURS ASSOCIES AU REMPLISSAGE
!
!   * ICPLP(I) = ADRESSE DANS LE RANGEMENT DES COEFFICIENTS A(I,J)
!               DU DERNIER COEFFICIENT DE REMPLISSAGE DE LA LIGNE I
!   * ICPCP(K) = POUR K = ICPL(I-1)+1...ICPL(I) NUMEROS DES INDICES DE
!               COLONNE J, DES COEFFICIENTS DE REMPLISSAGE DE LA LIGNE I
!               ( RANGES PAR ORDRE DE J CROISSANT)
!
!   ICPL,ICPC  : LES POINTEURS ASSOCIES A LA MATRICE FACTORISEE
!                ( REUNION DE ICPL ET ICPLP, ICPC ET ICPCP )
!   NZA        : NOMBRE DE COEFFICIENTS DE LA MATRICE FACTORISEE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    implicit none
#include "asterfort/pctrii.h"
    integer(kind=8) :: n
    integer(kind=4) :: icpc(*)
    integer(kind=8) :: icpl(0:n), icpd(n)
    integer(kind=8) :: icplp(0:n), icpcp(*), ind(n)
!
!=======================================================================
!
!     ON EFFECTUE UNE FACTORISATION LOGIQUE DE LA MATRICE A, LIGNE PAR
!     LIGNE, SANS ACTUALISATION DU REMPLISSAGE : ON NE CALCULE PAS LES
!      COEFFICIENTS L(I,J) ET U(I,J), MAIS LEUR POSITION (I,J).
!
!     LES FORMULES DE FACTORISATION UTILISEES :
!
! (1) L(I,JJ) = [ A(I,JJ) - SOM[J=1,JJ-1] L(I,J).U(J,JJ) ] / U(JJ,JJ)
!      POUR I=1,2,.....N ET JJ=1,2,....I
!
! (2) U(JJ,I) =   A(JJ,I) - SOM[J=1,JJ-1] L(JJ,J).U(JJ,J)
!      POUR I=1,2,.....N ET JJ=1,2,....I
!
!     ON PROCEDE AU CALCUL DES L(I,JJ) ET U(JJ,I) : IL Y A REMPLISSAGE
!     SI UN AU MOINS DES PRODUITS DE LA SOMME EST NON NUL
!
!     => POUR CHACUN DES L(I,J) NON NUL DE LA LIGNE I (BOUCLE 3: J < I)
!        ON RECHERCHE S'IL EXISTE UN U(J,JJ) NON NUL  (BOUCLE 4: J < JJ)
!        SI ON EN TROUVE UN, ON TESTE L'EXISTENCE DU COEFFICIENT A(I,JJ)
!        A L'AIDE DU TABLEAU IND. ON CREE UN COEFFICIENT DE REMPLISSAGE
!        SI NECESSAIRE.
!
!     => POUR CHACUN DES U(J,I) NON NUL DE LA LIGNE I (BOUCLE 3: J < I)
!        ON RECHERCHE S'IL EXISTE UN L(JJ,J) NON NUL  (BOUCLE 4: J < JJ)
!        SI ON EN TROUVE UN, ON TESTE L'EXISTENCE DU COEFFICIENT A(I,JJ)
!        A L'AIDE DU TABLEAU IND. ON CREE UN COEFFICIENT DE REMPLISSAGE
!        SI NECESSAIRE.
!
!    REMARQUE :  DANS LES FORMULES (1) ET (2) LES SOMMES SONT CALCULEES
!    --------    AVEC A(K,L) A LA PLACE DE L(K,L) OU U(K,L) CE QUI EST
!                COHERENT PUISQU'ON N'ACTUALISE PAS LE REMPLISSAGE
!
!=======================================================================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic1, ic2, ier, istop, j, jj
    integer(kind=8) :: k, k1, k2, kp1, kp2, l, lca
    integer(kind=8) :: ncremx, nzero
!-----------------------------------------------------------------------
    do i = 1, n
        ind(i) = 0
    end do
    ic1 = 0
    ic2 = 0
    k1 = 1
!
!
!     FACTORISATION LOGIQUE : LIGNE PAR LIGNE
!     ---------------------------------------
    do i = 1, n
!
!       TEST DE DEPASSEMENT DE DIMENSION (PROTECTION DES TABLEAUX)
!       REMARQUE JP : LA FICHE 12292 MONTRE QUE CE TEST EST INSUFFISANT
!       POUR NE PAS DEBORDER DU TABLEAU ICPCP, J'AI DU AJOUTER UN AUTRE
!       TEST AVANT ICPCP(IC1)=JJ (VOIR CI-DESSOUS)
        ncremx = i-icpc(k1)
        if (ic1+ncremx .gt. lca) then
            istop = i
            goto 90
        end if
!
!       MISE A JOUR DU TABLEAU IND
        k2 = icpl(i)
        do k = k1, k2
            j = icpc(k)
            ind(j) = i
        end do
        ind(i) = i
!
!       RECHERCHE DANS LA LIGNE I DES L(I,J) NON NULS
        do k = k1, icpd(i)
            j = icpc(k)
!         RECHERCHE DANS LA LIGNE J DES U(J,JJ) NON NULS
            do l = icpd(j)+1, icpl(j)
                jj = icpc(l)
!           LE COEFFICIENT L(I,JJ) OU U(JJ,I) EXISTE-T-IL ?
                if (ind(jj) .ne. i) then
!             NON ==> CREATION D'UN COEFFICIENT DE REMPLISSAGE
                    ic1 = ic1+1
!             STOCKAGE DE L'INDICE DE COLONNE DU COEFFICIENT LU(I,JJ)
                    if (ic1 .gt. lca) then
                        istop = i
                        goto 90
                    end if
                    icpcp(ic1) = jj
!             MISE A JOUR DU TABLEAU IND
                    ind(jj) = i
                end if
            end do
        end do
!
!       RECLASSEMENT DES INDICES DE COLONNE PAR ORDRE CROISSANT
        call pctrii(icpcp(ic2+1), ic1-ic2)
!
!       MISE A JOUR DU POINTEUR ICPLP
        icplp(i) = ic1
        ic2 = ic1
        k1 = k2+1
    end do
    icplp(0) = 0
!
!
!     TEST DE DEPASSEMENT DE DIMENSION (PROTECTION DES TABLEAUX)
    k1 = icpl(n)
    kp1 = icplp(n)
    nzero = k1+kp1
    if (nzero .gt. lca) then
        ier = nzero
        goto 140
    end if
!
!
!     CREATION DES TABLEAUX ICPL ET ICPC
!     POUR LA MATRICE FACTORISEE : REUNION DES TABLEAUX ICPC ET ICPCP
!     ---------------------------------------------------------------
    k = nzero
    do i = n, 1, -1
        icpl(i) = k
        kp2 = icplp(i-1)
        k2 = icpl(i-1)
60      continue
!
        if (k1 .gt. k2) then
!          LIGNE DE L EN COURS
            if (kp1 .gt. kp2) then
!           LIGNE DE COEF EN COURS
                if (icpc(k1) .lt. icpcp(kp1)) then
                    icpc(k) = int(icpcp(kp1), 4)
                    kp1 = kp1-1
                else
                    icpc(k) = icpc(k1)
                    k1 = k1-1
                end if
!
            else
!           LIGNE DE COEF FINIE
                icpc(k) = icpc(k1)
                k1 = k1-1
            end if
!
        else
!         LIGNE DE L FINIE
            if (kp1 .gt. kp2) then
!           LIGNE DE COEF EN COURS
                icpc(k) = int(icpcp(kp1), 4)
                kp1 = kp1-1
            else
                goto 70
            end if
        end if
!
        k = k-1
        goto 60
!
70      continue
    end do
    goto 140
!
!
90  continue
!     DEPASSEMENT DE DIMENSION ON CALCULE IC1= PLACE A AJOUTER
    do i = istop, n
        k2 = icpl(i)
        do k = k1, k2
            j = icpc(k)
            ind(j) = i
        end do
        ind(i) = i
        do k = k1, icpd(i)
            j = icpc(k)
            do l = icpd(j)+1, icpl(j)
                jj = icpc(l)
                if (jj .ge. i) goto 110
!
                if (ind(jj) .ne. i) then
!             NON ==> CREATION D'UN COEFFICIENT DE REMPLISSAGE
                    ic1 = ic1+1
                    ind(jj) = i
                end if
110             continue
            end do
        end do
        k1 = k2+1
    end do
!
!
!     NZERO=TAILLE MAT INI. + TAILLE MAT REMPLIE
    nzero = icpl(n)+ic1
    ier = nzero
140 continue
!
end subroutine
