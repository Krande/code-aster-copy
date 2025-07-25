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
subroutine mefver(ndim, som, xint, yint, rint)
! aslint: disable=
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndim(14)
    real(kind=8) :: som(9), xint(*), yint(*), rint(*)
!     VERIFICATION DE L'ORDRE ET DE LA BONNE DISPOSITION DES SOMMETS DE
!     L ENCEINTE RECTANGULAIRE
!     VERIFICATION DE L INCLUSION DES FAISCEAUX DANS L ENCEINTE
!     CIRCULAIRE
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
! ----------------------------------------------------------------------
    integer(kind=8) :: ind(3)
    real(kind=8) :: xsom(4), ysom(4), ux(4), uy(4), norm, a1, a(4)
    real(kind=8) :: vect(4), long(4)
    character(len=3) :: note
!     ------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iencei, j, nbcyl
    real(kind=8) :: diff, epsit, pi, pis2, proj, rext
    real(kind=8) :: xext, yext
!-----------------------------------------------------------------------
    nbcyl = ndim(3)
    iencei = ndim(6)
!
!
    pi = r8pi()
    pis2 = pi/2.d0
    epsit = 1.d-5
!
    if (iencei .eq. 2) then
        do i = 1, 4
            xsom(i) = som(2*i-1)
            ysom(i) = som(2*i)
        end do
!
!
! ---    MISE EN ORDRE DES SOMMETS DE L ENCEINTE
!
        ux(1) = xsom(2)-xsom(1)
        uy(1) = ysom(2)-ysom(1)
        ux(2) = xsom(3)-xsom(1)
        uy(2) = ysom(3)-ysom(1)
        ux(3) = xsom(4)-xsom(1)
        uy(3) = ysom(4)-ysom(1)
!
        do i = 2, 3
            norm = (ux(i)*ux(i)+uy(i)*uy(i))*(ux(1)*ux(1)+uy(1)*uy(1))
            norm = sqrt(norm)
            if (norm .eq. 0.d0) then
                call utmess('F', 'ALGELINE_88')
            end if
            a(i-1) = acos((ux(i)*ux(1)+uy(i)*uy(1))/norm)
            a1 = asin((ux(1)*uy(i)-uy(1)*ux(i))/norm)
            if (a1 .lt. 0.d0) a(i-1) = 2*pi-a(i-1)
        end do
!
        if (a(1) .lt. a(2) .and. a(2) .lt. pi) then
            ind(1) = 2
            ind(2) = 3
            ind(3) = 4
        else if (a(1) .gt. a(2) .and. a(1) .lt. pi) then
            ind(1) = 2
            ind(2) = 4
            ind(3) = 3
        else if (a(1) .lt. pis2 .and. a(2) .gt. pi) then
            ind(1) = 4
            ind(2) = 2
            ind(3) = 3
        else if (a(2) .lt. pis2 .and. a(1) .gt. pi) then
            ind(1) = 3
            ind(2) = 2
            ind(3) = 4
        else if (a(1) .lt. a(2) .and. a(1) .gt. pi) then
            ind(1) = 3
            ind(2) = 4
            ind(3) = 2
        else if (a(1) .gt. a(2) .and. a(2) .gt. pi) then
            ind(1) = 4
            ind(2) = 3
            ind(3) = 2
        else
            call utmess('F', 'ALGELINE_89')
        end if
!
        do i = 1, 3
            som(2*(i+1)-1) = xsom(ind(i))
            som(2*(i+1)) = ysom(ind(i))
        end do
!
! ---    ON VERIFIE QUE LES QUATRES SOMMETS FORMENT BIEN UN RECTANGLE
        do i = 1, 4
            xsom(i) = som(2*i-1)
            ysom(i) = som(2*i)
        end do
!
        ux(1) = xsom(2)-xsom(1)
        uy(1) = ysom(2)-ysom(1)
        ux(2) = xsom(4)-xsom(1)
        uy(2) = ysom(4)-ysom(1)
        ux(3) = xsom(4)-xsom(3)
        uy(3) = ysom(4)-ysom(3)
        ux(4) = xsom(3)-xsom(2)
        uy(4) = ysom(3)-ysom(2)
!
        do i = 1, 2
            vect(i) = ux(i)*uy(i+2)-uy(i)*ux(i+2)
        end do
!
        norm = (ux(2)*ux(2)+uy(2)*uy(2))*(ux(1)*ux(1)+uy(1)*uy(1))
        norm = sqrt(norm)
        if (norm .eq. 0.d0) then
            call utmess('F', 'ALGELINE_88')
        end if
        a(1) = acos((ux(2)*ux(1)+uy(2)*uy(1))/norm)
        if ((abs(a(1)-pis2)+abs(vect(1))+abs(vect(2))) .gt. epsit) then
            call utmess('F', 'ALGELINE_89')
        end if
!
!
! ---    VERIFICATION DE L INCLUSION DES CYLINDRES DANS L ENCEINTE
! ---    RECTANGULAIRE
!
! ---    NORMALISATION DES VECTEURS U(1) ET U(2)
        do i = 1, 2
            long(i) = sqrt(ux(i)*ux(i)+uy(i)*uy(i))
            ux(i) = ux(i)/long(i)
            uy(i) = uy(i)/long(i)
        end do
!
        do i = 1, nbcyl
            do j = 1, 2
                proj = ux(j)*(xint(i)-xsom(1))+uy(j)*(yint(i)-ysom(1))
                if ((proj-rint(i)) .lt. 0.d0 .or. (proj+rint(i)) .gt. long(j)) then
                    write (note(1:3), '(I3.3)') i
                    call utmess('F', 'ALGELINE_90', sk=note)
!
                end if
            end do
        end do
!
! ---    VERIFICATION DE L INCLUSION DES CYLINDRES DANS L ENCEINTE
! ---    CIRCULAIRE
!
    else if (iencei .eq. 1) then
        xext = som(1)
        yext = som(2)
        rext = som(3)
        do i = 1, nbcyl
            diff = sqrt((xext-xint(i))**2+(yext-yint(i))**2)
            if ((diff+rint(i)) .gt. rext) then
                write (note(1:3), '(I3.3)') i
                call utmess('F', 'ALGELINE_81', sk=note)
            end if
        end do
    end if
!
!
end subroutine
