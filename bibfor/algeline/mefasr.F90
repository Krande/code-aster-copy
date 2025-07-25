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
subroutine mefasr(ndim, nbcyl, nbgrp, nbtron, numgrp, &
                  idir, igrp, xint, yint, rint, &
                  sgn, orig, beta, a, b)
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/mefac2.h"
#include "asterfort/trigom.h"
    integer(kind=8) :: ndim(14), nbcyl, nbgrp, nbtron, numgrp(*), idir, igrp
    integer(kind=8) :: sgn(*), orig(*)
    real(kind=8) :: xint(*), yint(*), rint(*), beta(*)
    real(kind=8) :: a(2*nbtron*nbcyl, *), b(*)
!     ASSEMBLAGE POUR L ENCEINTE RECTANGULAIRE
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST, MEFREC
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
! IN  : XINT   : COORDONNEES 'X' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : YINT   : COORDONNEES 'Y' DES CENTRES DES CYLINDRES DANS
!                LE REPERE AXIAL
! IN  : RINT   : RAYONS DES CYLINDRES
! IN  : SGN    : -1 OU +1, COEFFICIENT INTERVENANT DANS LA DECOMPOSITION
!                EN SERIE DE LAURENT, SELON LE NIVEAU D IMAGE
! IN  : ORIG   : NUMERO DU CYLINDRE D ORIGINE DES CYLINDRES REELS OU
!                IMAGES
! IN  : BETA   : ANGLE CUMULE INTERVENANT DANS LA DECOMPOSITION EN
!                SERIE DE LAURENT, POUR LES CYLINDRES IMAGES
! IN  : A      : TABLEAU DE TRAVAIL: SOUS MATRICE DU SYSTEME A.X = B
! IN  : B      : TABLEAU DE TRAVAIL: SECOND MEMBRE DU SYSTEME A.X = B
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, l, nj, nl, m, nm
    real(kind=8) :: fic, dc, dc1
    real(kind=8) :: coef, coef1, coef2
    real(kind=8) :: rayk, rayi, arg
    real(kind=8) :: pi
! ----------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    integer(kind=8) :: nbfin, nbtot, nima, nima2
    real(kind=8) :: epsit
!-----------------------------------------------------------------------
    nbcyl = ndim(3)
    nbgrp = ndim(4)
    nbtron = ndim(5)
    nima = ndim(7)
    nima2 = ndim(8)
    nbtot = nbcyl*(2*nima+1)*(2*nima+1)
    nbfin = nbtot+4*(nima2)*(nima2+2*nima+1)
!
!
    pi = r8pi()
    epsit = 1.d-5
!
!
    do i = 1, nbcyl
        rayi = 1.d0/rint(i)
!
        do j = 1, nbtron
            nj = 2*j+2*nbtron*(i-1)
            rayi = rayi*rint(i)
!
            do k = 1, nbtot
                dc = sqrt( &
                     (xint(k)-xint(i))*(xint(k)-xint(i))+(yint(k)-yint(i))*(yint(k)-yint(i)))
!
                if (dc .ne. 0.d0) then
                    if (abs(abs(dc)-abs(xint(k)-xint(i))) .gt. epsit) then
                        fic = trigom('ACOS', (xint(k)-xint(i))/dc)
                        if ((yint(k)-yint(i)) .lt. 0.d0) then
                            fic = 2.d0*pi-fic
                        end if
                    else
                        fic = pi/2.d0*( &
                              1.d0-((xint(k)-xint(i))/dc)/abs(((xint(k)-xint(i))/dc)))
                        if ((yint(k)-yint(i)) .lt. 0.d0) then
                            fic = 2.d0*pi-fic
                        end if
                    end if
                else
                    fic = 0.d0
                end if
!
                if (k .ne. i) then
                    rayk = rint(k)
                    dc1 = dc**j
!
                    do l = 1, nbtron
                        rayk = rint(k)*rayk
                        dc1 = dc*dc1
                        nl = 2*l+2*nbtron*(orig(k)-1)
                        coef = mefac2(l, j)*rayi*rayk/dc1
                        coef = coef*((-1)**l)
                        arg = (j+l)*fic-l*beta(k)
                        a(nj-1, nl-1) = coef*cos(arg)+a(nj-1, nl-1)
                        a(nj, nl-1) = coef*sin(arg)+a(nj, nl-1)
                        a(nj-1, nl) = sgn(k)*coef*sin(arg)+a(nj-1, nl)
                        a(nj, nl) = -sgn(k)*coef*cos(arg)+a(nj, nl)
                    end do
!
                else
                    nl = 2*j+2*nbtron*(orig(k)-1)
                    a(nj-1, nl-1) = a(nj-1, nl-1)-j
                    a(nj, nl) = a(nj, nl)-j
                end if
            end do
!
            do k = nbtot+1, nbfin
                dc = sqrt( &
                     (xint(k)-xint(i))*(xint(k)-xint(i))+(yint(k)-yint(i))*(yint(k)-yint(i)))
!
                if (dc .ne. 0.d0) then
                    if (abs(abs(dc)-abs(xint(k)-xint(i))) .gt. epsit) then
                        fic = trigom('ACOS', (xint(k)-xint(i))/dc)
                        if ((yint(k)-yint(i)) .lt. 0.d0) then
                            fic = 2.d0*pi-fic
                        end if
                    else
                        fic = pi/2.d0*( &
                              1.d0-((xint(k)-xint(i))/dc)/abs(((xint(k)-xint(i))/dc)))
                        if ((yint(k)-yint(i)) .lt. 0.d0) then
                            fic = 2.d0*pi-fic
                        end if
                    end if
                else
                    fic = 0.d0
                end if
                dc1 = dc**j
!
                do l = 1, nbtron
                    dc1 = dc1*dc
                    coef = mefac2(l, j)*(rayi)/dc1
                    coef = coef*((-1)**l)
                    arg = (l+j)*fic-l*beta(k)
                    coef1 = coef*cos(arg)
                    coef2 = coef*sin(arg)
!
                    do m = 1, nbcyl
                        nm = 2*nbtron*(m-1)+2*l
                        arg = rint(m)**(l+1)
                        a(nj-1, nm-1) = coef1*arg+a(nj-1, nm-1)
                        a(nj, nm-1) = coef2*arg+a(nj, nm-1)
                        a(nj-1, nm) = sgn(k)*coef2*arg+a(nj-1, nm)
                        a(nj, nm) = -sgn(k)*coef1*arg+a(nj, nm)
                    end do
                end do
            end do
        end do
!
        if (numgrp(i) .eq. igrp) then
            b(2*nbtron*(i-1)+idir) = 1.d0
        end if
    end do
!
end subroutine
