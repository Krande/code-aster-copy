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
subroutine vpqlts(diag, surdia, neq, vecpro, mxcmp, &
                  mxiter, ier, nitqr)
    implicit none
#include "asterc/r8prem.h"
    integer(kind=8) :: neq, mxcmp, mxiter, ier, nitqr
    real(kind=8) :: diag(neq), surdia(neq), vecpro(mxcmp, neq)
!     CALCUL DE TOUTES LES VALEURS PROPRES ET DES VECTEURS PROPRES
!     ASSOCIES PAR LA METHODE QL-IMPLICITE POUR UNE MATRICE TRIDIAGONALE
!     SYMETRIQUE
!     ------------------------------------------------------------------
! VAR  DIAG  :    : TABLEAU DE R8 (1 .. NEQ)
!            EN ENTREE: VECTEUR CONTENANT LA DIAGONALE DE LA MATRICE
!            EN SORTIE: CONTIENT LES VALEURS PROPRES EN ORDRE QUELCONQUE
! VAR  SURDIA:    : TABLEAU DE R8 (1 .. NEQ)
!            EN ENTREE: CONTIENT LA  SUR-DIAGONALE     2,...,NEQ.
!            EN SORTIE: LA SUR-DIAGONALE EST PERDUE (ZONE DE TRAVAIL)
! IN   NEQ   : IS : ORDRE DE LA MATRICE
! OUT  VECPRO :   : TABLEAU DE R8 (1 .. NEQ)X(1 .. NEQ)
!            EN SORTIE:  CONTIENT LES VECTEURS PROPRES
!                        LA COLONNE J CONTIENT LE J-IEME VECTEUR PROPRE
!                        QUI CORRESPOND A LA VALEUR PROPRE DIAG(J).
! IN   MXCMP :  1-ERE DIMENSION (EXACTE) DE VECPRO (POUR LE FORTRAN)
!                   SI MXCMP < NEQ ALORS ON NE CALCULE PAS LES VECTEURS
!                                     ET VECPRO N'EST ALORS PAS UTILISE
! IN   MXITER:  NOMBRE MAXIMUM D'ITERATION (MXITER=30 EST UN BON CHOIX)
! OUT  IER   : IS : 0  PAS DE PROBLEME
!                   J  CONVERGENCE NON ATTEINTE A LA J VALEUR PROPRE
! OUT  NITQR : NOMBRE MAXIMAL ATTEINT D'ITERATIONS POUR AVOIR CONVERGE
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE ????
!     ------------------------------------------------------------------
    real(kind=8) :: b, c, f, g, h, p, r, s
    real(kind=8) :: zero, un, deux, epsmac
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ieq, j, jeq, jter, k, m
!
!-----------------------------------------------------------------------
    ier = 0
    zero = 0.0d0
    epsmac = r8prem()
    un = 1.0d0
    deux = 2.0d0
    nitqr = 0
!
    if (neq .eq. 1) goto 99999
!
!     --- TRANSFERT DES ELEMENTS SUR-DIAGONAUX ---
    do i = 2, neq
        surdia(i-1) = surdia(i)
    end do
!
!     --- SI NECESSAIRE CREATION DE LA MATRICE UNITE ---
    if (mxcmp .ge. neq) then
        do ieq = 1, neq
            do jeq = 1, neq
                vecpro(ieq, jeq) = zero
            end do
            vecpro(ieq, ieq) = un
        end do
    end if
!
    surdia(neq) = zero
    b = zero
    f = zero
    do j = 1, neq
        jter = 0
        h = epsmac*(abs(diag(j))+abs(surdia(j)))
        if (b .lt. h) b = h
!
!        --- RECHERCHE DU PLUS PETIT ELEMENT SUR-DIAGONAL ---
        do m = j, neq
            k = m
            if (abs(surdia(k)) .le. b) goto 15
        end do
15      continue
        m = k
        if (m .eq. j) goto 55
20      continue
        if (jter .eq. mxiter) goto 999
        jter = jter+1
        if (jter .gt. nitqr) then
            nitqr = jter
        end if
!
!        --- PREPARATION DU DECALAGE ---
        g = diag(j)
        p = (diag(j+1)-g)/(deux*surdia(j))
        r = abs(p)
        if (epsmac*abs(p) .lt. un) r = sqrt(p*p+un)
        diag(j) = surdia(j)/(p+sign(r, p))
        h = g-diag(j)
        do i = j+1, neq
            diag(i) = diag(i)-h
        end do
        f = f+h
!
!        --- ON APPLIQUE LA TRANSFORMATION QL ---
        p = diag(m)
        c = un
        s = zero
        do i = m-1, j, -1
            g = c*surdia(i)
            h = c*p
            if (abs(p) .ge. abs(surdia(i))) then
                c = surdia(i)/p
                r = sqrt(c*c+un)
                surdia(i+1) = s*p*r
                s = c/r
                c = un/r
            else
                c = p/surdia(i)
                r = sqrt(c*c+un)
                surdia(i+1) = s*surdia(i)*r
                s = un/r
                c = c*s
            end if
            p = c*diag(i)-s*g
            diag(i+1) = h+s*(c*g+s*diag(i))
            if (mxcmp .ge. neq) then
!              --- CALCUL DU VECTEUR PROPRE ---
                do k = 1, neq
                    h = vecpro(k, i+1)
                    vecpro(k, i+1) = s*vecpro(k, i)+c*h
                    vecpro(k, i) = c*vecpro(k, i)-s*h
                end do
            end if
        end do
        surdia(j) = s*p
        diag(j) = c*p
        if (abs(surdia(j)) .gt. b) goto 20
55      continue
        diag(j) = diag(j)+f
    end do
    goto 99999
!
999 continue
    ier = j
99999 continue
end subroutine
