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
subroutine caladj(col, diag, xadj, adjncy, n, &
                  nnz, deb, tab, suiv, lmat, &
                  ladjn, nrl)
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
#include "asterfort/utmess.h"
!      CALCUL DES VOISINS DE TOUS LES NOEUDS ( VERSION ASTER )
!      DONNEES
!      -------
!             DIAG(1:N)        : POINTEURS DE LA MATRICE RANGEE
!             COL(1:LMAT )         SYMETRIQUE INFERIEURE PAR LIGNES
!      RESULTATS
!      ---------
!             XADJ(1:N+1)     POINTEURS ,LES VOISINS DE I SONT ENTRE
!             ADJNCY(1:NADJ)     LES ADRESSES XADJ(I) ET XADJ(I+1)-1
!                                PAR CONTINUITE ON XADJ(N+1) = NADJ+1
!      TAB DE TRAVAIL
!      --------------
!             TAB(1:LT),SUIV(1:LT)
!             NNZ(1:N)
! ATTENTION : XADJ SERT DE TAB DE TRAVAIL DANS LA 1ERE PARTIE (DO 1 )
!-----------
    integer(kind=8) :: lmat, n, col(lmat), diag(0:n), adjncy(*)
    integer(kind=8) :: xadj(n+1), nnz(n), suiv(*), tab(*)
    integer(kind=8) :: nrl
    integer(kind=8) :: deb(n)
!--------------------------------------------------------------
!     VAR. LOCALES
    integer(kind=8) :: i, j, k, ii, it, ladjn, iad
    integer(kind=8) :: vali(2)
    if (nrl .eq. 0) then
!     PAS DE RELATION LINEAIRE ENTRE PLUSIEURS DDL
        do j = 1, n
!
            nnz(j) = diag(j)-diag(j-1)-1
        end do
        do k = 1, diag(n)
!     PARTIE TRIANGULAIRE SUPERIEURE
            nnz(col(k)) = nnz(col(k))+1
        end do
!
        xadj(1) = 1
        do j = 1, n
!     ON DIMINUE DE 1 CAR ON NE VEUT PAS LE TERME
!     DIAGONAL DANS ADJNCY
            xadj(j+1) = xadj(j)+nnz(j)-1
            nnz(j) = 0
        end do
!
        do j = 1, n
            do ii = diag(j-1)+1, diag(j)-1
                i = col(ii)
                adjncy(xadj(j)+nnz(j)) = i
                nnz(j) = nnz(j)+1
                adjncy(xadj(i)+nnz(i)) = j
                nnz(i) = nnz(i)+1
            end do
        end do
!
!     ---------------------------
    else
!     AVEC RELATION LINEAIRE
!     CALCUL DES LISTES DE NOEUDS A AJOUTER ( FAIT DANS PREMLA)
!
        do i = 1, n
!            DEB(I) =0
            nnz(i) = 0
        end do
!     INITIALISATION DE NNZ : NBRE DE TERMES A AJOUTER
!     POUR CHAQUE LIGNE
        do j = 1, n
            it = deb(j)
219         continue
            if (it .gt. 0) then
                nnz(j) = nnz(j)+1
                it = suiv(it)
                goto 219
            end if
        end do
!     VERIFICATION
        do j = 1, n
!     TERMES A AJOUTER PARTIE INFERIEURE
            nnz(j) = nnz(j)+diag(j)-diag(j-1)-1
        end do
        do k = 1, diag(n)
!     PARTIE TRIANGULAIRE SUPERIEURE
            nnz(col(k)) = nnz(col(k))+1
        end do
        do j = 1, n
!     TERMES A AJOUTER PARTIE SUPERIEURE
            it = deb(j)
324         continue
            if (it .gt. 0) then
                nnz(tab(it)) = nnz(tab(it))+1
                it = suiv(it)
                goto 324
            end if
        end do
!
        xadj(1) = 1
        do j = 1, n
!     ON DIMINUE DE 1 CAR ON NE VEUT PAS LE TERME
!     DIAGONAL DANS ADJNCY
            xadj(j+1) = xadj(j)+nnz(j)-1
            nnz(j) = 0
        end do
        if ((xadj(n+1)-1) .gt. ladjn) then
!       TEST D'ESPACE SUFFISANT DANS ADJNCY
            vali(1) = ladjn
            vali(2) = xadj(n+1)-1
            call utmess('F', 'ALGELINE4_4', ni=2, vali=vali)
        end if
!
        iad = 0
        do j = 1, n
            do ii = diag(j-1)+1, diag(j)-1
                i = col(ii)
                adjncy(xadj(j)+nnz(j)) = i
                nnz(j) = nnz(j)+1
                adjncy(xadj(i)+nnz(i)) = j
                nnz(i) = nnz(i)+1
                iad = max(iad, (xadj(i)+nnz(i)))
                iad = max(iad, (xadj(j)+nnz(j)))
            end do
            it = deb(j)
344         continue
            if (it .gt. 0) then
                adjncy(xadj(j)+nnz(j)) = tab(it)
                nnz(j) = nnz(j)+1
                adjncy(xadj(tab(it))+nnz(tab(it))) = j
                nnz(tab(it)) = nnz(tab(it))+1
                it = suiv(it)
                goto 344
            end if
        end do
!
    end if
end subroutine
