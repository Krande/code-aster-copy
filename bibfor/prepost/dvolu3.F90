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
function dvolu3(coord, norm, coord1)
!
!**********************************************************
!              BUT DE CETTE ROUTINE :                     *
! CALCULER LE VOLUME DE L INTERSECTION CYLINDRE-TETRAEDRE *
!**********************************************************
!
! IN   COORD  : COORDONNEES DES NOEUDS DU TETRA
!               ET DES INTERSECTIONS AVEC LE CYLINDRE
! IN   NORM   : POSITION DES NOEUDS PAR RAPPORT AU CYLINDRE
! IN   COORD1 : COORDONNEES DE O1
! OUT  DVOLU3 : VOLUME DE L INTERSECTION
!
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
!
! DECLARATION GLOBALE
!
    integer(kind=8) :: norm(2, 4)
    real(kind=8) :: coord(3, 12), coord1(3), dvolu3
!
! DECLARATION LOCALE
!
    integer(kind=8) :: i, j, k, l, m, e
    real(kind=8) :: vol3, vol4, xio1, yio1, zio1, dio1
    aster_logical :: lnoeu
!
! 1 - RECHERCHE DES DEUX POINTS INTERNES
!     RQ : 2 POINTS DEDANS
!          3 INTERSECTIONS
!          0 OU 1 POINT HORS PLAN
!
    i = 0
    do k = 1, 4
        if (norm(1, k) .eq. 1 .and. i .gt. 0) j = k
        if (norm(1, k) .eq. 1 .and. i .eq. 0) i = k
    end do
!
! 2 - RECHERCHE DU POINT HORS PLAN
!
    do k = 1, 4
        if (norm(2, k) .ne. 0) e = k
    end do
!
! 3 - NOEU1 ET NOEU2 SONT CONFONDUS AVEC LES 2 SOMMETS I ET J
! RECHERCHE DE LA CORRESPONDANCE SI LNOEU ALORS I=NOEU1 SINON I=NOEU2
!
    xio1 = coord(1, i)-coord1(1)
    yio1 = coord(2, i)-coord1(2)
    zio1 = coord(3, i)-coord1(3)
    dio1 = (xio1**2+yio1**2+zio1**2)
    if (dio1 .le. 1.0d-6) then
        if (norm(2, e) .eq. 1) then
            lnoeu = .true.
        else
            lnoeu = .false.
        end if
    else
        if (norm(2, e) .eq. -1) then
            lnoeu = .true.
        else
            lnoeu = .false.
        end if
    end if
!
! 4 - DETERMINATION DANS LE TABLEAU DES POSITIONS DES INTERSECTIONS
!
    if (i .eq. 1 .and. j .eq. 2) then
        if (lnoeu .and. e .eq. 4) then
            k = 8
            l = 6
            m = 9
        else if (.not. lnoeu .and. e .eq. 3) then
            k = 6
            l = 9
            m = 7
        else if (lnoeu .and. e .eq. 3) then
            k = 8
            l = 7
            m = 9
        else if (.not. lnoeu .and. e .eq. 4) then
            k = 6
            l = 8
            m = 7
        end if
    else if (i .eq. 1 .and. j .eq. 3) then
        if (lnoeu .and. e .eq. 2) then
            k = 10
            l = 7
            m = 8
        else if (.not. lnoeu .and. e .eq. 4) then
            k = 7
            l = 8
            m = 5
        else if (lnoeu .and. e .eq. 4) then
            k = 10
            l = 5
            m = 8
        else if (.not. lnoeu .and. e .eq. 2) then
            k = 7
            l = 10
            m = 5
        end if
    else if (i .eq. 1 .and. j .eq. 4) then
        if (lnoeu .and. e .eq. 3) then
            k = 9
            l = 5
            m = 10
        else if (.not. lnoeu .and. e .eq. 2) then
            k = 5
            l = 10
            m = 6
        else if (lnoeu .and. e .eq. 2) then
            k = 9
            l = 6
            m = 10
        else if (.not. lnoeu .and. e .eq. 3) then
            k = 5
            l = 9
            m = 6
        end if
    else if (i .eq. 2 .and. j .eq. 3) then
        if (lnoeu .and. e .eq. 4) then
            k = 6
            l = 5
            m = 10
        else if (.not. lnoeu .and. e .eq. 1) then
            k = 5
            l = 10
            m = 9
        else if (lnoeu .and. e .eq. 1) then
            k = 6
            l = 9
            m = 10
        else if (.not. lnoeu .and. e .eq. 4) then
            k = 5
            l = 6
            m = 9
        end if
    else if (i .eq. 2 .and. j .eq. 4) then
        if (lnoeu .and. e .eq. 1) then
            k = 10
            l = 8
            m = 7
        else if (.not. lnoeu .and. e .eq. 3) then
            k = 5
            l = 7
            m = 8
        else if (lnoeu .and. e .eq. 3) then
            k = 10
            l = 5
            m = 7
        else if (.not. lnoeu .and. e .eq. 1) then
            k = 8
            l = 10
            m = 5
        end if
    else if (i .eq. 3 .and. j .eq. 4) then
        if (lnoeu .and. e .eq. 2) then
            k = 7
            l = 6
            m = 9
        else if (.not. lnoeu .and. e .eq. 1) then
            k = 6
            l = 9
            m = 8
        else if (lnoeu .and. e .eq. 1) then
            k = 7
            l = 8
            m = 9
        else if (.not. lnoeu .and. e .eq. 2) then
            k = 6
            l = 7
            m = 8
        end if
    end if
!
! 5 - CALCUL DU VOLUME
!
    vol3 = (coord(1, l)-coord(1, i))*&
     &     ((coord(2, j)-coord(2, i))*(coord(3, k)-coord(3, i))-&
     &      (coord(3, j)-coord(3, i))*(coord(2, k)-coord(2, i)))
    vol3 = vol3+( &
         coord(2, l)-coord(2, i))*((coord(1, k)-coord(1, i))*(coord(3, j)-coord(3, i))-(coord(3, k)&
           &-coord(3, i))*(coord(1, j)-coord(1, i)) &
           )
    vol3 = vol3+( &
         coord(3, l)-coord(3, i))*((coord(1, j)-coord(1, i))*(coord(2, k)-coord(2, i))-(coord(2, j)&
           &-coord(2, i))*(coord(1, k)-coord(1, i)) &
           )
!
    if (abs(vol3) .lt. 1.0d-10) vol3 = 0.0d0
    if (vol3 .lt. 0.d0) then
        vol3 = -vol3
    end if
    vol4 = (coord(1, m)-coord(1, i))*&
     &     ((coord(2, j)-coord(2, i))*(coord(3, l)-coord(3, i))-&
     &      (coord(3, j)-coord(3, i))*(coord(2, l)-coord(2, i)))
    vol4 = vol4+( &
         coord(2, m)-coord(2, i))*((coord(1, l)-coord(1, i))*(coord(3, j)-coord(3, i))-(coord(3, l)&
           &-coord(3, i))*(coord(1, j)-coord(1, i)) &
           )
    vol4 = vol4+( &
         coord(3, m)-coord(3, i))*((coord(1, j)-coord(1, i))*(coord(2, l)-coord(2, i))-(coord(2, j)&
           &-coord(2, i))*(coord(1, l)-coord(1, i)) &
           )
!
    if (abs(vol4) .lt. 1.0d-10) vol4 = 0.0d0
    if (vol4 .lt. 0.d0) then
        vol4 = -vol4
    end if
    dvolu3 = vol3+vol4
    dvolu3 = dvolu3/6.d0
    if (dvolu3 .gt. 1.d+6) then
        call utmess('A', 'PREPOST_28')
    end if
!
end function
