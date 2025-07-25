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
subroutine mefger(ndim, som, xint, yint, rint, &
                  sgn, orig, beta)
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/trigom.h"
    integer(kind=8) :: ndim(14), sgn(*), orig(*)
    real(kind=8) :: som(9), xint(*), yint(*), rint(*), beta(*)
!     MISE EN FORME DES DONNEES POUR LA PRISE EN COMPTE DES CONDITIONS
!     AUX LIMITES PAR UNE METHODE DERIVEE DE LA METHODE DES IMAGES
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : SOM    : COORDONNEES DES SOMMETS DE L'ENCEINTE RECTANGULAIRE
!                OU XEXT,YEXT,REXT
! IN  : XINT   : COORDONNEES 'X' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : YINT   : COORDONNEES 'Y' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : RINT   : RAYONS DES CYLINDRES
! OUT : SGN    : -1 OU +1, COEFFICIENT INTERVENANT DANS LA DECOMPOSITION
!                EN SERIE DE LAURENT, SELON LE NIVEAU D IMAGE
! OUT : ORIG   : NUMERO DU CYLINDRE D ORIGINE DES CYLINDRES REELS OU
!                IMAGES
! OUT : BETA   : ANGLE CUMULE INTERVENANT DANS LA DECOMPOSITION EN
!                SERIE DE LAURENT, POUR LES CYLINDRES IMAGES
! ----------------------------------------------------------------------
    real(kind=8) :: xsom(4), ysom(4)
    real(kind=8) :: xcent, ycent
    real(kind=8) :: x12, y12, x23, y23
    real(kind=8) :: alph12, alph23, long12, long23
    real(kind=8) :: x0, y0
    real(kind=8) :: pi
! ----------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iencei, j, k, nbcyl, nbtot, nima
    integer(kind=8) :: nima2, nj, np
!-----------------------------------------------------------------------
    nbcyl = ndim(3)
    iencei = ndim(6)
    nima = ndim(7)
    nima2 = ndim(8)
    nbtot = nbcyl*(2*nima+1)*(2*nima+1)
!
!
    pi = r8pi()
!
! --- INITIALISATIONS
!
    do i = 1, nbtot
        beta(i) = 0.d0
        sgn(i) = 0
        orig(i) = 0
    end do
!
! --- CONSTRUCTION DES IMAGES
!
    do i = 1, nbcyl
        orig(i) = i
        sgn(i) = 1
        beta(i) = 0.0d0
    end do
!
    if (iencei .eq. 2) then
!
        do i = 1, 4
            xsom(i) = som(2*i-1)
            ysom(i) = som(2*i)
        end do
!
! ---    DEFINITION DES DROITES DE SYMETRIES
!
        x12 = xsom(2)-xsom(1)
        y12 = ysom(2)-ysom(1)
        long12 = x12*x12+y12*y12
!
        x23 = xsom(3)-xsom(2)
        y23 = ysom(3)-ysom(2)
        long23 = x23*x23+y23*y23
!
        if (x12/sqrt(long12) .gt. 1.d0) then
            alph12 = 0.d0
        else if (x12/sqrt(long12) .lt. -1.d0) then
            alph12 = pi
        else
            alph12 = trigom('ACOS', x12/sqrt(long12))
        end if
        if (y12 .lt. 0.d0) then
            alph12 = pi-alph12
        end if
!
        if (x23/sqrt(long23) .gt. 1.d0) then
            alph23 = 0.d0
        else if (x23/sqrt(long23) .lt. -1.d0) then
            alph23 = pi
        else
            alph23 = trigom('ACOS', x23/sqrt(long23))
        end if
        if (y23 .lt. 0.d0) then
            alph23 = pi-alph23
        end if
!
        xcent = (xsom(3)+xsom(1))/2.d0
        ycent = (ysom(3)+ysom(1))/2.d0
!
    end if
!
    nj = nbcyl
    x0 = xsom(1)
    y0 = ysom(1)
!
    do i = 1, nima
        do j = 1, nbcyl
            nj = nj+1
            if (i .eq. 1) then
                np = nj-nbcyl
            else
                np = nj-8*(i-1)*nbcyl
            end if
            xint(nj) = xint(np)-2.d0*x23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
            yint(nj) = yint(np)-2.d0*y23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
            beta(nj) = -beta(np)+2.d0*alph12
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = (-1)**i
            rint(nj) = rint(np)
            orig(nj) = orig(np)
        end do
!
        x0 = x0-x23
        y0 = y0-y23
!
        do j = 1, i
            x0 = x0+x12
            y0 = y0+y12
            do k = 1, nbcyl
                nj = nj+1
                np = nj-nbcyl
                xint(nj) = xint(np)-2.d0*x12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                yint(nj) = yint(np)-2.d0*y12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                beta(nj) = -beta(np)+2.d0*alph23
                beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
                sgn(nj) = -sgn(np)
                rint(nj) = rint(np)
                orig(nj) = orig(np)
            end do
        end do
!
        x0 = x0+x12
        y0 = y0+y12
!
        do j = 1, 2*i
            x0 = x0+x23
            y0 = y0+y23
            do k = 1, nbcyl
                nj = nj+1
                np = nj-nbcyl
                xint(nj) = xint(np)-2.d0*x23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
                yint(nj) = yint(np)-2.d0*y23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
                beta(nj) = -beta(np)+2.d0*alph12
                beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
                sgn(nj) = -sgn(np)
                rint(nj) = rint(np)
                orig(nj) = orig(np)
            end do
        end do
!
        x0 = x0+x23
        y0 = y0+y23
!
        do j = 1, 2*i
            x0 = x0-x12
            y0 = y0-y12
            do k = 1, nbcyl
                nj = nj+1
                np = nj-nbcyl
                xint(nj) = xint(np)-2.d0*x12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                yint(nj) = yint(np)-2.d0*y12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                beta(nj) = -beta(np)+2.d0*alph23
                beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
                sgn(nj) = -sgn(np)
                rint(nj) = rint(np)
                orig(nj) = orig(np)
            end do
        end do
!
        x0 = x0-x12
        y0 = y0-y12
!
        do j = 1, 2*i
            x0 = x0-x23
            y0 = y0-y23
            do k = 1, nbcyl
                nj = nj+1
                np = nj-nbcyl
                xint(nj) = xint(np)-2.d0*x23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
                yint(nj) = yint(np)-2.d0*y23/long23*(x23*(xint(np)-x0)+y23*(yint(np)-y0))
                beta(nj) = -beta(np)+2.d0*alph12
                beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
                sgn(nj) = -sgn(np)
                rint(nj) = rint(np)
                orig(nj) = orig(np)
            end do
        end do
!
        x0 = x0-x23
        y0 = y0-y23
!
        do j = 1, i-1
            x0 = x0+x12
            y0 = y0+y12
            do k = 1, nbcyl
                nj = nj+1
                np = nj-nbcyl
                xint(nj) = xint(np)-2.d0*x12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                yint(nj) = yint(np)-2.d0*y12/long12*(x12*(xint(np)-x0)+y12*(yint(np)-y0))
                beta(nj) = -beta(np)+2.d0*alph23
                beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
                sgn(nj) = -sgn(np)
                rint(nj) = rint(np)
                orig(nj) = orig(np)
            end do
        end do
!
        x0 = x0+x12
        y0 = y0+y12
!
    end do
!
!
    nj = nbtot
!
    do i = 1, nima2
        nj = nj+1
        xint(nj) = xcent-(nima+i)*x23
        yint(nj) = ycent-(nima+i)*y23
        if (i .eq. 1) then
            beta(nj) = -beta(nj-8*nima*nbcyl)+2.d0*alph12
        else
            beta(nj) = -beta(nj-8*(nima+i-1))+2.d0*alph12
        end if
        beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
        sgn(nj) = (-1)**(nima+i)
!
        do j = 1, nima+i
!
            nj = nj+1
            xint(nj) = xint(nj-1)+x12
            yint(nj) = yint(nj-1)+y12
            beta(nj) = -beta(nj-1)+2.d0*alph23
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = -sgn(nj-1)
        end do
!
        do j = 1, 2*(nima+i)
            nj = nj+1
            xint(nj) = xint(nj-1)+x23
            yint(nj) = yint(nj-1)+y23
            beta(nj) = -beta(nj-1)+2.d0*alph12
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = -sgn(nj-1)
        end do
!
        do j = 1, 2*(nima+i)
            nj = nj+1
            xint(nj) = xint(nj-1)-x12
            yint(nj) = yint(nj-1)-y12
            beta(nj) = -beta(nj-1)+2.d0*alph23
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = -sgn(nj-1)
        end do
!
        do j = 1, 2*(nima+i)
            nj = nj+1
            xint(nj) = xint(nj-1)-x23
            yint(nj) = yint(nj-1)-y23
            beta(nj) = -beta(nj-1)+2.d0*alph12
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = -sgn(nj-1)
        end do
!
        do j = 1, nima+i-1
            nj = nj+1
            xint(nj) = xint(nj-1)+x12
            yint(nj) = yint(nj-1)+y12
            beta(nj) = -beta(nj-1)+2.d0*alph23
            beta(nj) = beta(nj)-int(beta(nj)/2.d0/pi)*2.d0*pi
            sgn(nj) = -sgn(nj-1)
        end do
!
    end do
!
!
end subroutine
