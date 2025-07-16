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
subroutine mefasc(ndim, nbcyl, nbgrp, nbtron, numgrp, &
                  idir, igrp, som, rint, dcent, &
                  ficent, d, fi, a, b)
    implicit none
!
#include "asterfort/mefac1.h"
#include "asterfort/mefac2.h"
    integer(kind=8) :: ndim(14), nbcyl, nbgrp, nbtron, numgrp(*), idir, igrp
    real(kind=8) :: dcent(nbcyl), ficent(nbcyl), rint(*), som(9)
    real(kind=8) :: d(nbcyl, nbcyl), fi(nbcyl, nbcyl)
    real(kind=8) :: a(2*nbtron*(nbcyl+1), *), b(*)
!     ASSEMBLAGE POUR L ENCEINTE CIRCULAIRE
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST, MEFCIR
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : NBCYL  : NOMBRE DE CYLINDRES
! IN  : NBGRP  : NOMBRE DE GROUPES D EQUIVALENCE
! IN  : NBTRON : ORDRE DE TRONCATURE DES SERIES DE LAURENT DANS LA BASE
!                MODALE
! IN  : NUMGRP : INDICES DES GROUPES D EQUIVALENCE
! IN  : IDIR   : INDICES DE CYLINDRE
! IN  : IGRP   : INDICES DE GROUPE DE CYLINDRE
! IN  : SOM    : XEXT,YEXT,REXT
! IN  : RINT   : RAYONS DES CYLINDRES
! IN  : DCENT  : DISTANCE DU CENTRE DES CYLINDRES AU CENTRE DE
!                L ENCEINTE
! IN  : FICENT : ANGLE POLAIRE PAR RAPPORT AU CENTRE DE L ENCEINTE
! IN  : D      : DISTANCE RELATIVE ENTRE LES CENTRES DES CYLINDRES
! IN  : FI     : ANGLE POLAIRE RELATIF PAR RAPPORT AU CENTRE DE CHAQUE
!                CYLINDRE
! IN  : A      : TABLEAU DE TRAVAIL: SOUS MATRICE DU SYSTEME A.X = B
! IN  : B      : TABLEAU DE TRAVAIL: SECOND MEMBRE DU SYSTEME A.X = B
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, l, ni, nj, nk, nl
    real(kind=8) :: coef
! ----------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    real(kind=8) :: rext
!-----------------------------------------------------------------------
    nbcyl = ndim(3)
    nbgrp = ndim(4)
    nbtron = ndim(5)
!
!
    rext = som(3)
!
!
    do j = 1, nbtron
        nj = 2*j
        do k = 1, nbcyl
            nk = 2*nbtron*k
            do l = 1, j
                nl = nk+2*l
                if (dcent(k) .eq. 0.d0 .and. j .eq. l) then
                    coef = mefac1(j, l)*(rint(k)**(l+1))/(rext**(j+1))
                else if (dcent(k) .eq. 0.d0 .and. j .ne. l) then
                    coef = 0.d0
                else
                    coef = mefac1(j, l)*(dcent(k)**(j-l))*(rint(k)**(l+1))/(rext**(j+1))
                end if
                a(nj-1, nl-1) = -coef*cos((j-l)*ficent(k))
                a(nj, nl-1) = coef*sin((j-l)*ficent(k))
                a(nj-1, nl) = coef*sin((j-l)*ficent(k))
                a(nj, nl) = coef*cos((j-l)*ficent(k))
            end do
        end do
        a(nj-1, nj-1) = j
        a(nj, nj) = -j
!
    end do
!
    do i = 1, nbcyl
        ni = 2*nbtron*i
        do j = 1, nbtron
            nj = ni+2*j
            do k = 1, nbcyl
                nk = 2*nbtron*k
                if (k .ne. i) then
                    do l = 1, nbtron
                        nl = nk+2*l
                        coef = mefac2(l, j)*(rint(i)**(j-1))*(rint(k)**(l+1))/(d(i, k)**(l+j))
                        coef = coef*((-1)**l)
                        a(nj-1, nl-1) = coef*cos((j+l)*fi(i, k))
                        a(nj, nl-1) = coef*sin((j+l)*fi(i, k))
                        a(nj-1, nl) = coef*sin((j+l)*fi(i, k))
                        a(nj, nl) = -coef*cos((j+l)*fi(i, k))
                    end do
                else
                    nl = nk+2*j
                    a(nj-1, nl-1) = -j
                    a(nj, nl) = -j
                end if
            end do
!
            do l = j, nbtron
                nl = 2*l
                if (dcent(i) .eq. 0.d0 .and. j .eq. l) then
                    coef = mefac1(l, j)*(rint(i)**(j-1))/(rext**(l-1))
                else if (dcent(i) .eq. 0.d0 .and. j .ne. l) then
                    coef = 0.d0
                else
                    coef = mefac1(l, j)*(rint(i)**(j-1))*(dcent(i)**(l-j))/(rext**(l-1))
                end if
                a(nj-1, nl-1) = coef*cos((l-j)*ficent(i))
                a(nj, nl-1) = -coef*sin((l-j)*ficent(i))
                a(nj-1, nl) = coef*sin((l-j)*ficent(i))
                a(nj, nl) = coef*cos((l-j)*ficent(i))
            end do
!
        end do
        if (numgrp(i) .eq. igrp) then
            b(2*nbtron*i+idir) = 1.d0
        end if
    end do
!
end subroutine
