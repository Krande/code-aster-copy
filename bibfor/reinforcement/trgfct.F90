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

subroutine trgfct(nb, fcttab)
    !_____________________________________________________________________
    !
    !     TRGFCT
    !
    !      CALCUL DES PARAMETRES TRIGONOMÉTRIQUES POUR LES (36) FACETTES
    !
    !_____________________________________________________________________
    !
    !
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8rddg.h"

    !
    !     ! FACETTES POUR METHODE DE CAPRA ET MAURY
    !
    !       nb FACETTES
    integer(kind=8) :: nb
    !       NOMBRE DE DIVISIONS ENTRE -PI/2 ET +PI/2
    real(kind=8) :: fcttab(nb, 6)
    !
    !       ANGLE DE LA FACETTE (-PI/2 <= X < +PI/2)
    real(kind=8) :: angle
    real(kind=8) :: pas, pi
    integer(kind=8) :: i
    !
    pi = r8pi()
    fcttab(1, 1) = 0d0
    fcttab(1, 2) = 1d0
    fcttab(1, 3) = 0d0
    fcttab(1, 4) = 0d0
    fcttab(1, 5) = -1d0
    fcttab(1, 6) = -90.d0
    !
    !       -PI/2
    angle = -pi/2.d0
    !       2PI/N
    pas = -2.d0*angle/nb
    !
    !   POUR CHAQUE FACETTE, LES VALEURS SONT
    !     1 = COS^2
    !     2 = SIN^2
    !     3 = 2 SIN COS
    !     4 = SIN
    !     5 = COS
    do 20 i = 2, nb
        angle = angle+pas
        fcttab(i, 1) = cos(angle)*cos(angle)
        fcttab(i, 2) = sin(angle)*sin(angle)
        fcttab(i, 3) = -2.d0*sin(angle)*cos(angle)
        fcttab(i, 4) = cos(angle)
        fcttab(i, 5) = sin(angle)
        fcttab(i, 6) = angle*r8rddg()
20      continue
        end subroutine
