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

subroutine fgequi(tz, typz, ndim, equi)
    implicit none
#include "asterfort/jacobi.h"
#include "asterfort/lchydr.h"
#include "asterfort/lciv2e.h"
#include "asterfort/lciv2s.h"
#include "asterfort/lcqeqv.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim
    real(kind=8) :: tz(*), equi(*)
    character(len=*) :: typz
! person_in_charge: thomas.de-soza at edf.fr
!
!     BUT:
!       CALCULER LES GRANDEURS EQUIVALENTES SUIVANTES
!       . CONTRAINTES EQUIVALENTES               (= 17 VALEURS)
!          . VON MISES                               (= 1 VALEUR)
!          . TRESCA                                  (= 1 VALEUR)
!          . CONTRAINTES PRINCIPALES                 (= 3 VALEURS)
!          . VON-MISES * SIGNE (PRESSION)            (= 1 VALEUR)
!          . DIRECTIONS DES CONTRAINTES PRINCIPALES  (= 3*3 VALEURS)
!          . TRACE                                   (= 1 VALEUR)
!          . TAUX DE TRIAXIALITE                     (= 1 VALEUR)
!
!       . DEFORMATIONS EQUIVALENTES              (= 14 VALEURS)
!          . SECOND INVARIANT                        (= 1 VALEUR)
!          . DEFORMATIONS PRINCIPALES                (= 3 VALEURS)
!          . 2EME INV. * SIGNE (1ER.INV.)            (= 1 VALEUR)
!          . DIRECTIONS DES DEFORMATIONS PRINCIPALES (= 3*3 VALEURS)
!
! ----------------------------------------------------------------------
!     IN     TZ    TENSEUR CONTRAINTE OU DEFORMATION (XX YY ZZ XY XZ YZ)
!            TYPZ  TYPE DU TENSEUR 'SIGM' OU 'EPSI' SUIVI EVENTUELLEMENT
!            DES CARACTERES _DIR POUR CALCULER LES VECTEURS DIRECTEURS
!            NDIM  DIMENSION ESPACE 3 OU 2
!     OUT    EQUI  VECTEUR DES GRANDEURS EQUIVALENTES
! ----------------------------------------------------------------------
    real(kind=8) :: t(6), tn(6), tr(6), tu(6), vecp(3, 3), nul(6)
    real(kind=8) :: rac2, hyd, jacaux(3)
    real(kind=8) :: tol, toldyn
    integer(kind=8) :: nbvec, nperm
    integer(kind=8) :: type, iordre
    character(len=8) :: typ
    common/tdim/nt, nd
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, nd, nitjac, nt
!-----------------------------------------------------------------------
    data nul/6*0.d0/
    data nperm, tol, toldyn/12, 1.d-10, 1.d-2/
! ----------------------------------------------------------------------
    typ = typz
!
    if (ndim .eq. 3) then
        nt = 6
        nd = 3
    else if (ndim .eq. 2) then
        nt = 4
        nd = 3
    else if (ndim .eq. 1) then
        nt = 1
        nd = 1
    end if
!
    call r8inir(6, 0.d0, t, 1)
    do i = 1, nt
        t(i) = tz(i)
    end do
!
! --- TENSEUR TN = (XX YY ZZ RAC2.XY RAC2.XZ RAC2.YZ) (POUR LCIV2E)
!
    rac2 = sqrt(2.d0)
    do i = 1, nd
        tn(i) = t(i)
    end do
    do i = nd+1, nt
        tn(i) = rac2*t(i)
    end do
!
! --- MATRICE TR = (XX XY XZ YY YZ ZZ) (POUR JACOBI)
!
    tr(1) = t(1)
    tr(2) = t(4)
    tr(3) = t(5)
    tr(4) = t(2)
    tr(5) = t(6)
    tr(6) = t(3)
!
! --- MATRICE UNITE = (1 0 0 1 0 1) (POUR JACOBI)
!
    tu(1) = 1.d0
    tu(2) = 0.d0
    tu(3) = 0.d0
    tu(4) = 1.d0
    tu(5) = 0.d0
    tu(6) = 1.d0
!
! --- VALEURS PRINCIPALES
    nbvec = 3
!
    do i = 1, nbvec
        do j = 1, nbvec
            vecp(j, i) = 0.0d0
        end do
    end do
!
! --- DEFORMATIONS
!
    if (typ(1:4) .eq. 'EPSI') then
! ------ SECOND INVARIANT
        equi(1) = lciv2e(tn)
        if (lcqeqv(tr, nul) .eq. 'OUI') then
            equi(2) = 0.d0
            equi(3) = 0.d0
            equi(4) = 0.d0
        else
            type = 0
            iordre = 0
            call jacobi(nbvec, nperm, tol, toldyn, tr, &
                        tu, vecp, equi(2), jacaux, nitjac, &
                        type, iordre)
        end if
! ------ PREMIER INVARIANT
        call lchydr(tn, hyd)
! ------ EQUIVALENT FATIGUE = SECOND INVARIANT * SIGNE(PREMIER INV)
        if (hyd .ge. 0.d0) equi(5) = equi(1)
        if (hyd .lt. 0.d0) equi(5) = -equi(1)
!
! ------ DIRECTION DES DEFORMATIONS PRINCIPALES DANS EQUI
! -      DANS L ORDRE LES COORDONNEES DU VECTEUR PUIS PASSAGE
! -      A L AUTRE VECTEUR
        if (typ(5:8) .eq. '_DIR') then
            do i = 1, nbvec
                do j = 1, nbvec
                    equi(5+((i-1)*nbvec)+j) = vecp(j, i)
                end do
            end do
        end if
!
! --- CONTRAINTES
!
    else if (typ(1:4) .eq. 'SIGM') then
! ------ VON MISES = SECOND INVARIANT
        equi(1) = lciv2s(tn)
!
        if (lcqeqv(tr, nul) .eq. 'OUI') then
            equi(3) = 0.d0
            equi(4) = 0.d0
            equi(5) = 0.d0
        else
            type = 0
            iordre = 0
            call jacobi(nbvec, nperm, tol, toldyn, tr, &
                        tu, vecp, equi(3), jacaux, nitjac, &
                        type, iordre)
        end if
! ------ TRESCA = MAX DIFF VALEURS PRINCIPALES
        equi(2) = max(abs(equi(3)-equi(4)), abs(equi(3)-equi(5)), abs(equi(4)-equi(5)))
!
! ------ PREMIER INVARIANT
        call lchydr(tn, hyd)
! ------ EQUIVALENT FATIGUE = SECOND INVARIANT * SIGNE(PREMIER INV)
        if (hyd .ge. 0.d0) equi(6) = equi(1)
        if (hyd .lt. 0.d0) equi(6) = -equi(1)
!
! ------ DIRECTION DES CONTRAINTES PRINCIPALES DANS EQUI
! -      DANS L ORDRE LES COORDONNEES DU VECTEUR PUIS PASSAGE
! -      A L AUTRE VECTEUR
        if (typ(5:8) .eq. '_DIR') then
            do i = 1, nbvec
                do j = 1, nbvec
                    equi(6+((i-1)*nbvec)+j) = vecp(j, i)
                end do
            end do
!
! ------    TRACE DES CONTRAINTES : TRSIG
            equi(16) = 3.d0*hyd
!
! ------    TRIAXIALITE DES CONTRAINTES : TRIAX
!
            if (equi(1) .gt. tol*abs(equi(16))) then
                equi(17) = hyd/equi(1)
            else
                equi(17) = 0.d0
            end if
        end if
    end if
!
end subroutine
