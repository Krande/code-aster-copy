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
subroutine pcinfe(n, icpl, icpc, icpd, icplp, &
                  icpcp, ind, lca, ier)
!       S.P. PCINFE IDEM S-P PCFULL
!                    SAUF NON CREATION COEFS SUR U
!
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
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ic1, ic2, ier, istop, j, jj
    integer(kind=8) :: k, k1, k2, kp1, kp2, l, lca
    integer(kind=8) :: nzero
!-----------------------------------------------------------------------
    do i = 1, n
        ind(i) = 0
    end do
    ic1 = 0
    ic2 = 0
    k1 = 1
!
!     FACTORISATION LOGIQUE : LIGNE PAR LIGNE
!     ---------------------------------------
!
    do i = 1, n
        k2 = icpl(i)
!
!     MISE A JOUR DU TABLEAU IND
!
        do k = k1, k2
            j = icpc(k)
            ind(j) = i
        end do
        ind(i) = i
!
!     RECHERCHE DANS LA LIGNE I DES L(I,J) NON NULS
!
        do k = k1, icpd(i)
            j = icpc(k)
!
!     RECHERCHE DANS LA LIGNE J DES U(J,JJ) NON NULS
!
            do l = icpd(j)+1, icpl(j)
                jj = icpc(l)
!
!     LE COEFFICIENT L(I,JJ) EXISTE-T-IL ?
!                         ARRET AU PREMIER  U(JJ,I)
!
                if (jj .ge. i) goto 30
!
                if (ind(jj) .ne. i) then
!
!     NON ==> CREATION D'UN COEFFICIENT DE REMPLISSAGE
!
                    ic1 = ic1+1
!
!     TEST DE DEPASSEMENT DE DIMENSION (PROTECTION DES TABLEAUX)
!
                    if (ic1 .gt. lca) then
!               WRITE (6,107) NIV,LCA,I
                        istop = i
                        goto 100
                    end if
!
!     STOCKAGE DE L'INDICE DE COLONNE DU COEFFICIENT LU(I,JJ)
!
                    icpcp(ic1) = jj
!
!     MISE A JOUR DU TABLEAU IND
!
                    ind(jj) = i
                end if
30              continue
            end do
        end do
!
!     RECLASSEMENT DES INDICES DE COLONNE PAR ORDRE CROISSANT
!
        call pctrii(icpcp(ic2+1), ic1-ic2)
!
!     MISE A JOUR DU POINTEUR ICPLP
!
        icplp(i) = ic1
        ic2 = ic1
        k1 = k2+1
    end do
    icplp(0) = 0
!
!     AVANT FUSION DE ICPC ET ICPCP
!     TEST DE DEPASSEMENT DE DIMENSION (PROTECTION DES TABLEAUX)
!
    k1 = icpl(n)
    kp1 = icplp(n)
    nzero = k1+kp1
    if (nzero .gt. lca) then
!       WRITE (6,200) NIV,LCA,NZERO
        ier = nzero
        goto 150
    end if
!
!     CREATION DES TABLEAUX ICPL ET ICPC
!     POUR LA MATRICE FACTORISEE : REUNION DES TABLEAUX ICPC ET ICPCP
!     ---------------------------------------------------------------
!
    k = nzero
    do i = n, 1, -1
        icpl(i) = k
        kp2 = icplp(i-1)
        k2 = icpl(i-1)
60      continue
        if (k1 .gt. k2) then
            if (kp1 .gt. kp2) then
!       -------------------
                if (icpc(k1) .lt. icpcp(kp1)) then
                    icpc(k) = int(icpcp(kp1), 4)
                    k = k-1
                    kp1 = kp1-1
                else
                    icpc(k) = icpc(k1)
                    k = k-1
                    k1 = k1-1
                end if
            else
                icpc(k) = icpc(k1)
                k = k-1
                k1 = k1-1
            end if
        else
!     ---- LIGNE DE L EPUISEE ------
70          continue
            if (kp1 .gt. kp2) then
                icpc(k) = int(icpcp(kp1), 4)
                k = k-1
                kp1 = kp1-1
            else
                goto 80
            end if
            goto 70
        end if
!     ------
        goto 60
80      continue
    end do
!
!     LE NOMBRE DE COEFFICIENTS DE LA MATRICE FACTORISEE
!
!     NZCA = NZERO
!     WRITE (6,*) ' FIN DU S-P PCINFE  TAILLE FACTORISEE= ',NZCA
!
    goto 150
!
!  DEPASSEMENT DE DIMENSION ON CALCULE IC1= PLACE A AJOUTER
!
100 continue
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
                if (jj .ge. i) goto 120
!
                if (ind(jj) .ne. i) then
!     NON ==> CREATION D'UN COEFFICIENT DE REMPLISSAGE
                    ic1 = ic1+1
                    ind(jj) = i
                end if
120             continue
            end do
        end do
        k1 = k2+1
    end do
! NZERO=TAILLE MAT INI.+TAILLE MAT REMPLIE
    nzero = icpl(n)+ic1
!     WRITE (6,200) NIV,LCA,NZERO
    ier = nzero
150 continue
!
end subroutine
