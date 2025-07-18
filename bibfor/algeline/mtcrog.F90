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
subroutine mtcrog(a, b, nmax, n, nbsc, &
                  c, wks, ier)
    implicit none
!
#include "jeveux.h"
#include "asterfort/mtcro1.h"
#include "asterfort/mtcro2.h"
#include "asterfort/mtcro3.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nmax, n, nbsc, ier
    real(kind=8) :: a(nmax, n), b(nmax, nbsc), c(nmax, nbsc), wks(nmax)
!     RESOLUTION APPROCHEE D UN SYSTEME LINEAIRE AX = B PAR LA METHODE
!     DE CROUT, POUR UNE MATRICE A QUELCONQUE DE DIMENSION N*N
!     SI B EST DE DIMENSION N*1, IL S AGIT D UN SIMPLE SYSTEME
!     LINEAIRE. SI B EST DE DIMENSION N*N ET VAUT L IDENTITE, IL S AGIT
!     DE L INVERSION D UNE MATRICE
! ----------------------------------------------------------------------
!     OPERATEUR APPELANT: OP0144, FLUST3, MEFIST, MEFEIG, MEFREC, MEFCIR
! ----------------------------------------------------------------------
! IN  : NMAX   : DIMENSIONS DES TABLEAUX
! IN  : N      : DIMENSIONS EFFECTIVE DES MATRICES
! IN  : NBSC   : DIMENSIONS DU SYSTEM A RESOUDRE
! IN  : A      : MATRICE A DU SYSTEME LINEAIRE (MODIFIEE LORS DE LA
!                RESOLUTION)
! IN  : B      : MATRICE B DU SYSTEME LINEAIRE
! OUT : C      : SOLUTION DU SYSTEME LINEAIRE (SI IER = 0)
!                SI B EST L IDENTITE ET NBSC = N, C EST L INVERSE DE A
! --  : WKS    : TABLEAU DE TRAVAIL
! OUT : IER    : INDICE D ERREUR (0: RESOLUTION CORRECTE
!                                 1: PAS DE CONVERGENCE OU ERREUR)
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: id, i, j, k, l, iapro
    real(kind=8) :: det, prec, x, y
! ----------------------------------------------------------------------
!
!
! --- PREC EST LA PRECISION MACHINE.
! --- PREC   : LE PLUS PETIT REEL POSITIF, TEL QUE: 1.0 + PREC > 1.0
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    prec = .1110223024625156651d-15
!
    ier = 0
!
! --- TEST DE VERIFICATION DES DIMENSIONS
!
    if (nmax .lt. n) then
        call utmess('F', 'ALGELINE2_14')
    end if
!
!
! --- SAUVEGARDE DE LA MATRICE B
!
    do i = 1, n
        do j = 1, nbsc
            c(i, j) = b(i, j)
        end do
    end do
!
! --- A EST DECOMPOSEE EN A = LU, OU L EST UNE MATRICE TRIANGULAIRE
! --- INFERIEURE, ET U, UNE MATRICE TRIANGULAIRE SUPERIEURE UNITAIRE.
! --- L ET U SONT ECRITES A LA PLACE DE A
! --- WKS EST LA TABLE DES PERMUTATIONS EFFECTUEES
! --- SI LA MATRICE A EST SINGULIERE OU PRESQUE SINGULIERE, IER = 1.
! --- IER = 0 SINON, ET DET EST LE DETERMINANT
!
!
    do i = 1, n
        wks(i) = 0.0d0
    end do
!
    do j = 1, n
        do i = 1, n
            wks(i) = wks(i)+a(i, j)**2
        end do
    end do
!
    do i = 1, n
        if (wks(i) .le. 0.0d0) goto 240
        wks(i) = 1.0d0/sqrt(wks(i))
    end do
!
    det = 1.0d0
    id = 0
    do k = 1, n
        l = k
        x = 0.0d0
        do i = k, n
            y = abs(a(i, k)*wks(i))
            if (y .le. x) goto 100
            x = y
            l = i
100         continue
        end do
        if (l .ne. k) then
            det = -det
            do j = 1, n
                y = a(k, j)
                a(k, j) = a(l, j)
                a(l, j) = y
            end do
            wks(l) = wks(k)
        end if
        wks(k) = l
        det = det*a(k, k)
        if (x .lt. 8.0d0*prec) then
            ier = 1
            goto 240
        end if
160     continue
        if (abs(det) .lt. 1.0d0) goto 180
        det = det*0.0625d0
        id = id+4
        goto 160
180     continue
        if (abs(det) .ge. 0.0625d0) goto 200
        det = det*16.0d0
        id = id-4
        goto 180
200     continue
        if (k .lt. n) then
            call mtcro1(k, a, nmax, a(1, k+1))
            call mtcro3(n-k, k, a(k+1, 1), nmax, a(1, k+1), &
                        a(k+1, k+1))
        end if
    end do
240 continue
!
!
! --- RESOLUTION D UN SYSTEME LINEAIRE AX = B
! --- B UNE EST UNE MATRICE DE DIMENSION NMAX*NBSC
! --- A EST SOUS LA FORME A = LU, OU L EST UNE MATRICE TRIANGULAIRE
! --- INFERIEURE, ET U, UNE MATRICE TRIANGULAIRE SUPERIEURE UNITAIRE.
! --- L ET U SONT ECRITES A LA PLACE DE A
! --- P EST LA TABLE DES PERMUTATIONS EFFECTUEES AUPARAVANT.
! --- AX = B EST RESOLU EN TROIS ETAPES, PERMUTATION DES ELEMENTS DE B,
! --- CALCUL DE LY = B, ET CALCUL DE  UX = Y. LES MATRICES Y ET X SONT
! --- ECRITES A LA PLACE DE B
!
    if (ier .eq. 0) then
! ---    PERMUTATION DES ELEMENTS DE B
!
        do i = 1, n
            iapro = int(wks(i)+0.5d0)
            if (iapro .eq. i) goto 340
            do k = 1, nbsc
                x = b(i, k)
                b(i, k) = b(iapro, k)
                b(iapro, k) = x
            end do
340         continue
        end do
        do k = 1, nbsc
! ---       RECHERCHE DE LA SOLUTION DE LY = B
            call mtcro1(n, a, nmax, b(1, k))
! ---       RECHERCHE DE LA SOLUTION DE UX = Y
            call mtcro2(n, a, nmax, b(1, k))
        end do
!
! --- RESTAURATION DE LA MATRICE B, ET DE LA SOLUTION
!
        do i = 1, n
            do j = 1, nbsc
                x = b(i, j)
                b(i, j) = c(i, j)
                c(i, j) = x
            end do
        end do
!
    else
        ier = 1
    end if
!
!
end subroutine
