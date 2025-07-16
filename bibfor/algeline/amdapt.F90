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
!> Brief description of the subroutine,
!> continuation line
!>
!
!> Detailed description
!>
!> The lines must start with '!>' to be extracted by doxygen.
!>
!> Paragraphs are separated by a blank line.
!> Without blank line, the following lines are concatenated in the same paragraph.
!>
!> End of line may be inserted to start\n a new line.
!>
!> Equations can be inserted in the document using LaTeX syntax.
!>
!> Compute \f$ \frac{d\lambda}{dt}, \frac{d\phi}{dt},  \frac{dz}{dt} \f$
!
!> @todo all todo comments will be visible in a single page.
!
!> @param[in]  neq      number of equations
!> @param[out] nnv      the output numbering
!
subroutine amdapt(neq, nbnd, nbsn, pe, nv, &
                  invp, parent, supnd, adress, lgind, &
                  fctnzs, fctops, llist, nnv)
! person_in_charge: olivier.boiteau at edf.fr
!
!     DONNEES
!     NEQ : NBRE TOTAL D'INCONNUES (LAGRANGES INCLUS)
!     NBND : NBRE DE DDL  RENUMEROTES PAR AMDBAR (NEQ - LAGRANGES)
!     INVP : INVP(I) : NOUVEAU NUMERO DU DDL  (RESULTATS DE AMDBAR)
!     PE  (RESULTATS DE AMDBAR)
!     NV (RESULTATS DE AMDBAR)
!     ON OUTPUT, PE HOLDS THE ASSEMBLY TREE/FOREST, WHICH IMPLICITLY
!     REPRESENTS A PIVOT ORDER WITH IDENTICAL FILL-IN AS THE ACTUAL
!     ORDER (VIA A DEPTH-FIRST SEARCH OF THE TREE).
!
!     ON OUTPUT:
!     IF NV (I) .GT. 0, THEN I REPRESENTS A NODE IN THE ASSEMBLY TREE,
!     AND THE PARENT OF I IS -PE (I), OR ZERO IF I IS A ROOT.
!     IF NV (I) = 0, THEN (I,-PE (I)) REPRESENTS AN EDGE IN A
!     SUBTREE, THE ROOT OF WHICH IS A NODE IN THE ASSEMBLY TREE.
!
!
!     RESULTATS
!     SUPND(1:NBSN+1) : DEFINITION DES SUPERNOEUDS, LE SN I EST DEFINI
!     PAR ( SUPND(I), SUPND(I+1)-1)
!     PARENT(NBSN) : PARENT(I) EST LE SN PARENT DU SN I
!     ADRESS : POINTEUR DANS LES TABLEAUX GLOBAL ET LOCAL
!     EN FAIT ADRESS SERA RECALCULE DANS FACSMB, IL SERT ICI
!     A EVALUER LGIND.
!
!     LGIND : LONGUEUR DE  GLOBAL (ET LOCAL),
!     FCTNZS      NBRE DE TERMES NON NULS DANS LA FACTORISEE,
!     FCTOPS   NBRE DE D OPERATIONS.
!     (CES VALEURS SONT INDICATIVES)
!     CAR NE PRENANT PAS EN COMPTE LES LAGRANGES.
!
!     TAB DE TRAVAIL
!     LLIST,NNV
!
    implicit none
!
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: neq, invp(neq), pe(neq+1), nv(neq)
    integer(kind=8) :: nbnd, nbsn, lgind, fctnzs
    real(kind=8) :: fctops
    integer(kind=8) :: parent(*), supnd(neq), adress(*)
!
    integer(kind=8) :: llist(neq), nnv(neq), ifm, niv
    integer(kind=8) :: i, j, k, ncol, nlig, deb, fin, snj, ndi
!
!     NNV EQUIVAUDRA A NV DANS LA NOUVELLE NUMEROTATION
!     PARENT SERA  LE PARENT PAR INCONNUE ET NON PAR SUPERNOEUD
    call infniv(ifm, niv)
    nbsn = 0
    do i = 1, nbnd
        if (nv(i) .ne. 0) nbsn = nbsn+1
        j = invp(i)
        if (pe(i) .ne. 0) then
            parent(j) = invp(-pe(i))
        else
            parent(j) = 0
        end if
        nnv(j) = nv(i)
    end do
!     NV CONTIENDRA  LA LARGEUR DES SN (LGSN AILLEURS)
    j = 1
    do i = 1, nbsn
        nv(i) = 1
111     continue
        if (nnv(j) .eq. 0) then
            nv(i) = nv(i)+1
            j = j+1
            goto 111
        end if
        j = j+1
    end do
!
    supnd(1) = 1
    do i = 1, nbsn
        supnd(i+1) = supnd(i)+nv(i)
    end do
!     LLIST SERA  L INVERSE DE SUPND, APPELE INVSUP AILLEURS
    k = 0
    adress(1) = 1
    do snj = 1, nbsn
        do i = 1, nv(snj)
            k = k+1
            llist(k) = snj
        end do
!     CALCUL DE ADRESS : ON CHERCHE DANS LE SN SNJ,
!     LE NOEUD I T.Q. NNV(I) =/= 0, CETTE VALEUR CONTIENT
!     LE "VRAI DEGRE" DE I,(SORTIE DE AMDBAR)
!     ON RAPPELLE VRAI DEGRE = LARGEUR DU SN + NBRE DE VOISINS
!     C A D ADRESS(SN+1) - ADRESS(SN)
        deb = supnd(snj)
        fin = supnd(snj+1)-1
        i = fin
128     continue
        if (i .lt. deb) then
            call utmess('F', 'ALGELINE5_6', si=snj)
        end if
        if (nnv(i) .ne. 0) goto 129
        i = i-1
        goto 128
129     continue
        adress(snj+1) = adress(snj)+nnv(i)
    end do
!
    lgind = adress(nbsn+1)-1
    fctnzs = 0
    fctops = 0.d0
    do i = 1, nbsn
        ncol = supnd(i+1)-supnd(i)
        nlig = adress(i+1)-adress(i)
        fctnzs = fctnzs+nlig*ncol-(ncol*(ncol+1))/2
        do j = 1, ncol
            nlig = nlig-1
            fctops = fctops+nlig*(nlig+3)
        end do
    end do
!     ON CALCUL PARENT EN SN A PARTIR DE PARENT/NOEUDS
!     NV CONTIENDRA PROVISOIREMENT PARENT/NOEUDS
    do i = 1, nbnd
        nv(i) = parent(i)
    end do
!
    do i = 1, nbnd
        parent(i) = 0
    end do
    do snj = 1, nbsn
        deb = supnd(snj)
        fin = supnd(snj+1)-1
        ndi = fin
175     continue
        if (ndi .lt. deb) then
            call utmess('F', 'ALGELINE5_6', si=snj)
        end if
        if (nnv(ndi) .ne. 0) goto 177
        ndi = ndi-1
        goto 175
177     continue
!        LE NOEUD NDI EST LE ND REPRESENTATIF DU SN
!        PARENT(SNJ) EST LE SN CONTENANT LE PARENT DE NDI
        if (nv(ndi) .ne. 0) then
            parent(snj) = llist(nv(ndi))
        else
            parent(snj) = 0
        end if
!
    end do
    if (niv .eq. 2) then
        write (ifm, *) 'AMDAPT  :  TRAITEMENT DE AMDBAR'
        write (ifm, *) ' NOMBRE DE SUPERNOEUDS: ', nbsn
        do i = 1, nbsn
            write (ifm, *) 'SN ', i, ' :NDS DE ', supnd(i), ' A ', supnd(i+1) &
                -1, ',PARENT ', parent(i), ',DEGRE :', adress(i+1)-adress(i)
        end do
    end if
!
end subroutine
