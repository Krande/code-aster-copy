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
subroutine conare(typma, ar, nbar)
    implicit none
#include "asterfort/utmess.h"
    integer(kind=8) :: ar(12, 3), nbar
    character(len=8) :: typma
!                       RENVOIE LA MATRICE DE CONNECTIVITÉ DES
!                       ARETES D'UNE MAILLE DE TYPE TYPMA
!
!    ENTREE :
!              TYPMA : TYPE DE LA MAILLE (TYPE_MAILLE)
!
!    SORTIE :
!              AR   : MATRICE DE CONNECTIVITÉ DES ARETES
!              NBAR : NOMBRE D'ARETES
!......................................................................
!
    integer(kind=8) :: i, j
!......................................................................
!
    do i = 1, 12
        do j = 1, 3
            ar(i, j) = 0
        end do
    end do
!
    if (typma .eq. 'HEXA8') then
        nbar = 12
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE HEXA8
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(4, 1) = 4
        ar(4, 2) = 1
        ar(5, 1) = 5
        ar(5, 2) = 6
        ar(6, 1) = 6
        ar(6, 2) = 7
        ar(7, 1) = 7
        ar(7, 2) = 8
        ar(8, 1) = 8
        ar(8, 2) = 5
        ar(9, 1) = 1
        ar(9, 2) = 5
        ar(10, 1) = 2
        ar(10, 2) = 6
        ar(11, 1) = 3
        ar(11, 2) = 7
        ar(12, 1) = 4
        ar(12, 2) = 8
    else if (typma .eq. 'HEXA20' .or. typma .eq. 'HEXA27') then
        nbar = 12
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE HEXA20 / HEXA27
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 9
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(2, 3) = 10
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(3, 3) = 11
        ar(4, 1) = 4
        ar(4, 2) = 1
        ar(4, 3) = 12
        ar(5, 1) = 5
        ar(5, 2) = 6
        ar(5, 3) = 17
        ar(6, 1) = 6
        ar(6, 2) = 7
        ar(6, 3) = 18
        ar(7, 1) = 7
        ar(7, 2) = 8
        ar(7, 3) = 19
        ar(8, 1) = 8
        ar(8, 2) = 5
        ar(8, 3) = 20
        ar(9, 1) = 1
        ar(9, 2) = 5
        ar(9, 3) = 13
        ar(10, 1) = 2
        ar(10, 2) = 6
        ar(10, 3) = 14
        ar(11, 1) = 3
        ar(11, 2) = 7
        ar(11, 3) = 15
        ar(12, 1) = 4
        ar(12, 2) = 8
        ar(12, 3) = 16
    else if (typma .eq. 'PENTA6') then
        nbar = 9
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE PENTA6 OU PENTA15
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(3, 1) = 3
        ar(3, 2) = 1
        ar(4, 1) = 4
        ar(4, 2) = 5
        ar(5, 1) = 5
        ar(5, 2) = 6
        ar(6, 1) = 6
        ar(6, 2) = 4
        ar(7, 1) = 1
        ar(7, 2) = 4
        ar(8, 1) = 2
        ar(8, 2) = 5
        ar(9, 1) = 3
        ar(9, 2) = 6
    else if (typma .eq. 'PENTA15' .or. typma .eq. 'PENTA18') then
        nbar = 9
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE PENTA15 OU PENTA18
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 7
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(2, 3) = 8
        ar(3, 1) = 3
        ar(3, 2) = 1
        ar(3, 3) = 9
        ar(4, 1) = 4
        ar(4, 2) = 5
        ar(4, 3) = 13
        ar(5, 1) = 5
        ar(5, 2) = 6
        ar(5, 3) = 14
        ar(6, 1) = 6
        ar(6, 2) = 4
        ar(6, 3) = 15
        ar(7, 1) = 1
        ar(7, 2) = 4
        ar(7, 3) = 10
        ar(8, 1) = 2
        ar(8, 2) = 5
        ar(8, 3) = 11
        ar(9, 1) = 3
        ar(9, 2) = 6
        ar(9, 3) = 12
    else if (typma .eq. 'PYRAM5') then
        nbar = 8
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE PYRAM5 OU PYRAM13
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(4, 1) = 1
        ar(4, 2) = 4
        ar(5, 1) = 1
        ar(5, 2) = 5
        ar(6, 1) = 2
        ar(6, 2) = 5
        ar(7, 1) = 3
        ar(7, 2) = 5
        ar(8, 1) = 4
        ar(8, 2) = 5
    else if (typma .eq. 'PYRAM13') then
        nbar = 8
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE PYRAM5 OU PYRAM13
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 6
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(2, 3) = 7
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(3, 3) = 8
        ar(4, 1) = 1
        ar(4, 2) = 4
        ar(4, 3) = 9
        ar(5, 1) = 1
        ar(5, 2) = 5
        ar(5, 3) = 10
        ar(6, 1) = 2
        ar(6, 2) = 5
        ar(6, 3) = 11
        ar(7, 1) = 3
        ar(7, 2) = 5
        ar(7, 3) = 12
        ar(8, 1) = 4
        ar(8, 2) = 5
        ar(8, 3) = 13
    else if (typma .eq. 'TETRA4') then
        nbar = 6
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE TETRA4 OU TETRA10
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 1
        ar(2, 2) = 3
        ar(3, 1) = 1
        ar(3, 2) = 4
        ar(4, 1) = 2
        ar(4, 2) = 3
        ar(5, 1) = 2
        ar(5, 2) = 4
        ar(6, 1) = 3
        ar(6, 2) = 4
    else if (typma .eq. 'TETRA10') then
        nbar = 6
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE TETRA4 OU TETRA10
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 5
        ar(2, 1) = 1
        ar(2, 2) = 3
        ar(2, 3) = 7
        ar(3, 1) = 1
        ar(3, 2) = 4
        ar(3, 3) = 8
        ar(4, 1) = 2
        ar(4, 2) = 3
        ar(4, 3) = 6
        ar(5, 1) = 2
        ar(5, 2) = 4
        ar(5, 3) = 9
        ar(6, 1) = 3
        ar(6, 2) = 4
        ar(6, 3) = 10
    else if (typma .eq. 'QUAD4') then
        nbar = 4
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE QUAD4
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(4, 1) = 4
        ar(4, 2) = 1
    else if (typma .eq. 'QUAD8' .or. typma .eq. 'QUAD9') then
        nbar = 4
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE QUAD8
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 5
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(2, 3) = 6
        ar(3, 1) = 3
        ar(3, 2) = 4
        ar(3, 3) = 7
        ar(4, 1) = 4
        ar(4, 2) = 1
        ar(4, 3) = 8
    else if (typma .eq. 'TRIA3') then
        nbar = 3
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE TRIA3
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(3, 1) = 3
        ar(3, 2) = 1
    else if (typma .eq. 'TRIA6' .or. typma .eq. 'TRIA7') then
        nbar = 3
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE TRIA6
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 4
        ar(2, 1) = 2
        ar(2, 2) = 3
        ar(2, 3) = 5
        ar(3, 1) = 3
        ar(3, 2) = 1
        ar(3, 3) = 6
    else if (typma .eq. 'SEG2') then
        nbar = 1
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE SEG2
        ar(1, 1) = 1
        ar(1, 2) = 2
    else if (typma .eq. 'SEG3') then
        nbar = 1
!       CONNECTIVITÉ DES ARETES POUR UNE MAILLE SEG3
        ar(1, 1) = 1
        ar(1, 2) = 2
        ar(1, 3) = 3
    else if (typma .eq. 'POI1') then
        nbar = 0
    else
        call utmess('F', 'ALGORITH2_23', sk=typma)
    end if
!
end subroutine
