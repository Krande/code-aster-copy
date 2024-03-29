! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine xmcoor(jcesd, jcesv, jcesl, ifiss, ndim, &
                  npte, nummae, ifac, xp, yp, &
                  coord)
!
!
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer :: jcesd(10), jcesv(10), jcesl(10)
    integer :: ndim, nummae, ifac, npte, ifiss
    real(kind=8) :: xp, yp
    real(kind=8) :: coord(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (CONTACT - GRANDS GLISSEMENTS)
!
! CALCUL DES COORDONNEES LOCAUX DU PT D'INTEGRATION OU DE SON PROJETE
! DANS LES MAILLES MAITRES OU ESCLAVES RESPECTIVEMENT
!
! TRAVAIL EFFECTUE EN COLLABORATION AVEC L'I.F.P.
!
! ----------------------------------------------------------------------
!
!
!
!  JCES*(3)  : POINTEURS DE LA SD SIMPLE DES COOR DES PT D'INTER
!  JCES*(4)  : POINTEURS DE LA SD SIMPLE DE CONNECTIVITÉ DES FACETTES
! IN IFISS  : NUMÉRO DE FISSURE LOCALE DANS NUMMAE
! IN  NDIM  : DIMENSION DU PROBLEME
! IN  NUMMAE: POSITION DE LA MAILLE ESCLAVE OU MAITRE
! IN  IFAC  : NUMERO LOCAL DE LA FACETTE ESCLAVE OU MAITRE
! IN  XP    : COORDONNEE X DU POINT D'INTEGRATION DE CONTACT SUR
!             LA MAILLE ESCLAVE OU MAITRE
! IN  YP    : COORDONNEE Y DU POINT D'INTEGRATION DE CONTACT SUR
!             LA MAILLE ESCLAVE OU MAITRE
! OUT COORD : COORDONNEES DU POINT D'INTEGRATION DANS L'ELEMENT
!             PARENT
!
!
!
!
    real(kind=8) :: coor(npte)
    integer :: i, j, iad, numpi(npte)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- RECUPERATION DES NUM LOCAUX DES PTS D'INTER DE LA FACETTE
!
    coor(:) = 0.d0
    do i = 1, npte
        numpi(i) = 0
    end do
    coord(:) = 0.d0
!
    do i = 1, npte
        call cesexi('S', jcesd(4), jcesl(4), nummae, 1, &
                    ifiss, (ifac-1)*ndim+i, iad)
        ASSERT(iad .gt. 0)
        numpi(i) = zi(jcesv(4)-1+iad)
    end do
    do i = 1, ndim
! --- BOUCLE SUR LES DIMENSIONS
        do j = 1, npte
! --- BOUCLE SUR LES POINTS D'INTERSECTIONS
! --- RECUPERATION DE LA COMPOSANTE LOCALE I DE CHACUN DES POINTS
! --- D'INTERSECTIONS J DE LA FACETTE
            call cesexi('S', jcesd(3), jcesl(3), nummae, 1, &
                        ifiss, ndim*(numpi(j)-1)+i, iad)
            ASSERT(iad .gt. 0)
            coor(j) = zr(jcesv(3)-1+iad)
        end do
! --- CALCUL DE LA COMPOSANTE I POUR LE POINT DE CONTACT DANS LA
! --- MAILLE PARENTE
        if (ndim .eq. 2) then
            if (npte .le. 2) then
                coord(i) = coor(1)*(1-xp)/2+coor(2)*(1+xp)/2
            else if (npte .eq. 3) then
                coord(i) = -coor(1)*xp*(1-xp)/2+coor(2)*xp*(1+xp)/2+coor(3)*(1+xp)*(1-xp)
            end if
        else if (ndim .eq. 3) then
            coord(i) = coor(1)*(1-xp-yp)+coor(2)*xp+coor(3)*yp
        end if
    end do
!
    call jedema()
end subroutine
