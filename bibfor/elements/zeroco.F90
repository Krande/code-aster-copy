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

subroutine zeroco(x, y)
!
    implicit none
#include "asterc/r8maem.h"
#include "asterfort/utmess.h"
#include "asterfort/zerodi.h"
    real(kind=8) :: x(4), y(4)
! ----------------------------------------------------------------------
!  RESOLUTION D'EQUATIONS SCALAIRES PAR UNE METHODE DE CORDES
! ----------------------------------------------------------------------
! VAR X(1) BORNE DE L'INTERVALLE DE RECHERCHE  TQ Y(1) < 0
! VAR X(2) BORNE DE L'INTERVALLE DE RECHERCHE  TQ Y(2) > 0
! VAR X(3) SOLUTION X(N-1) PUIS SOLUTION EN X(N)
! VAR X(4) SOLUTION X(N)   PUIS SOLUTION EN X(N+1)
! VAR Y(I) VALEUR DE LA FONCTION EN X(I)
! ----------------------------------------------------------------------
!
    real(kind=8) :: xp, p, s
    real(kind=8) :: infini
!
!    TEST DES PRE-CONDITIONS
    if (y(1) .gt. 0 .or. y(2) .lt. 0) then
        call utmess('F', 'ELEMENTS4_61')
    end if
!
!
!    TRAITEMENT DES SITUATIONS AVEC BORNES INFINIES -> DICHOTOMIE
    infini = r8maem()
    if (abs(y(3)) .eq. infini .or. abs(y(4)) .eq. infini) then
        call zerodi(x, y)
        goto 9999
    end if
!
!    REACTUALISATION DE L'INTERVALLE DE RECHERCHE
    if (y(4) .lt. 0.d0) then
        x(1) = x(4)
        y(1) = y(4)
    end if
!
    if (y(4) .gt. 0.d0) then
        x(2) = x(4)
        y(2) = y(4)
    end if
!
!    CONSTRUCTION D'UN NOUVEL ESTIME
    if (x(3) .eq. x(4)) then
        call utmess('A', 'ALGORITH9_84')
        goto 9999
    end if
    p = (y(4)-y(3))/(x(4)-x(3))
    if (p .ne. 0.d0) then
        xp = x(3)-y(3)/p
        s = (xp-x(1))/(x(2)-x(1))
    else
        s = -2.d0
    end if
!
!    CORRECTION DU NOUVEL ESTIME SI NECESSAIRE (EN DEHORS DES BORNES)
    if (s .le. 0.d0 .or. s .ge. 1.d0) then
        if (abs(y(1)) .eq. infini .or. abs(y(2)) .eq. infini) then
            xp = (x(1)+x(2))/2
        else
            p = (y(2)-y(1))/(x(2)-x(1))
            xp = x(1)-y(1)/p
        end if
    end if
!
!    DECALAGE DES ITERES
    x(3) = x(4)
    x(4) = xp
    y(3) = y(4)
!
9999 continue
end subroutine
