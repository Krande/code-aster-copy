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

subroutine xdecoupe(elp, cnset, nse, nnose)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"

    integer(kind=8)          :: cnset(:), nse, nnose
    character(len=8) :: elp

! person_in_charge: samuel.geniaut at edf.fr
!
!                      CONNECTIVITÉ DES ÉLÉMENTS TETRAS À PARTIR
!                               D'UN ÉLÉMENT PARENT X-FEM
!                          (VOIR BOOK III 19/04/04)
!
!     ENTREE
!       ELP     : TYPE DE MAILLE
!
!     SORTIE
!       CNSET   : CONNECTIVITÉ DES NOEUDS DE LA MAILLE
!       NSE     : NOMBRE DE SOUS-TÉTRAS (SOUS TRIA)
!       NNOSE   : NOMBRE DE NOEUDS DU SOUS TETRA (SOUS TRIA)
!     ------------------------------------------------------------------
!
    integer(kind=8) :: i, ino, ise, connec(144, 4), connec1(12, 6)
! ----------------------------------------------------------------------

    call jemarq()
    if (elp .eq. 'HE8') then
!       Nombre de pentahèdre dans un hexahèdre
!       1
        connec1(1, 1) = 4
        connec1(1, 2) = 1
        connec1(1, 3) = 2
        connec1(1, 4) = 3
        connec1(1, 5) = 5
        connec1(1, 6) = 6
!       2
        connec1(2, 1) = 5
        connec1(2, 2) = 8
        connec1(2, 3) = 7
        connec1(2, 4) = 6
        connec1(2, 5) = 4
        connec1(2, 6) = 3
!       3
        connec1(3, 1) = 7
        connec1(3, 2) = 6
        connec1(3, 3) = 5
        connec1(3, 4) = 8
        connec1(3, 5) = 2
        connec1(3, 6) = 1
!       4
        connec1(4, 1) = 2
        connec1(4, 2) = 3
        connec1(4, 3) = 4
        connec1(4, 4) = 1
        connec1(4, 5) = 7
        connec1(4, 6) = 8
!       5
        connec1(5, 1) = 2
        connec1(5, 2) = 1
        connec1(5, 3) = 5
        connec1(5, 4) = 6
        connec1(5, 5) = 4
        connec1(5, 6) = 8
!       6
        connec1(6, 1) = 4
        connec1(6, 2) = 3
        connec1(6, 3) = 7
        connec1(6, 4) = 8
        connec1(6, 5) = 2
        connec1(6, 6) = 6
!       7
        connec1(7, 1) = 7
        connec1(7, 2) = 8
        connec1(7, 3) = 4
        connec1(7, 4) = 3
        connec1(7, 5) = 5
        connec1(7, 6) = 1
!       8
        connec1(8, 1) = 5
        connec1(8, 2) = 6
        connec1(8, 3) = 2
        connec1(8, 4) = 1
        connec1(8, 5) = 7
        connec1(8, 6) = 3
!       9
        connec1(9, 1) = 3
        connec1(9, 2) = 4
        connec1(9, 3) = 1
        connec1(9, 4) = 2
        connec1(9, 5) = 8
        connec1(9, 6) = 5
!       10
        connec1(10, 1) = 8
        connec1(10, 2) = 7
        connec1(10, 3) = 6
        connec1(10, 4) = 5
        connec1(10, 5) = 3
        connec1(10, 6) = 2
!       11
        connec1(11, 1) = 6
        connec1(11, 2) = 5
        connec1(11, 3) = 8
        connec1(11, 4) = 7
        connec1(11, 5) = 1
        connec1(11, 6) = 4
!       12
        connec1(12, 1) = 1
        connec1(12, 2) = 2
        connec1(12, 3) = 3
        connec1(12, 4) = 4
        connec1(12, 5) = 6
        connec1(12, 6) = 7

        nse = 12
        nnose = 6

!       12 tétraèdres par pentahèdre
        do i = 1, nse
            connec(((i-1)*12)+1, 1) = connec1(i, 1)
            connec(((i-1)*12)+1, 2) = connec1(i, 2)
            connec(((i-1)*12)+1, 3) = connec1(i, 4)
            connec(((i-1)*12)+1, 4) = connec1(i, 5)

            connec(((i-1)*12)+2, 1) = connec1(i, 1)
            connec(((i-1)*12)+2, 2) = connec1(i, 3)
            connec(((i-1)*12)+2, 3) = connec1(i, 4)
            connec(((i-1)*12)+2, 4) = connec1(i, 5)

            connec(((i-1)*12)+3, 1) = connec1(i, 3)
            connec(((i-1)*12)+3, 2) = connec1(i, 4)
            connec(((i-1)*12)+3, 3) = connec1(i, 5)
            connec(((i-1)*12)+3, 4) = connec1(i, 6)

            connec(((i-1)*12)+4, 1) = connec1(i, 2)
            connec(((i-1)*12)+4, 2) = connec1(i, 4)
            connec(((i-1)*12)+4, 3) = connec1(i, 5)
            connec(((i-1)*12)+4, 4) = connec1(i, 6)

            connec(((i-1)*12)+5, 1) = connec1(i, 2)
            connec(((i-1)*12)+5, 2) = connec1(i, 3)
            connec(((i-1)*12)+5, 3) = connec1(i, 4)
            connec(((i-1)*12)+5, 4) = connec1(i, 6)

            connec(((i-1)*12)+6, 1) = connec1(i, 1)
            connec(((i-1)*12)+6, 2) = connec1(i, 2)
            connec(((i-1)*12)+6, 3) = connec1(i, 3)
            connec(((i-1)*12)+6, 4) = connec1(i, 5)

            connec(((i-1)*12)+7, 1) = connec1(i, 1)
            connec(((i-1)*12)+7, 2) = connec1(i, 2)
            connec(((i-1)*12)+7, 3) = connec1(i, 5)
            connec(((i-1)*12)+7, 4) = connec1(i, 6)

            connec(((i-1)*12)+8, 1) = connec1(i, 1)
            connec(((i-1)*12)+8, 2) = connec1(i, 3)
            connec(((i-1)*12)+8, 3) = connec1(i, 5)
            connec(((i-1)*12)+8, 4) = connec1(i, 6)

            connec(((i-1)*12)+9, 1) = connec1(i, 1)
            connec(((i-1)*12)+9, 2) = connec1(i, 3)
            connec(((i-1)*12)+9, 3) = connec1(i, 4)
            connec(((i-1)*12)+9, 4) = connec1(i, 6)

            connec(((i-1)*12)+10, 1) = connec1(i, 1)
            connec(((i-1)*12)+10, 2) = connec1(i, 2)
            connec(((i-1)*12)+10, 3) = connec1(i, 4)
            connec(((i-1)*12)+10, 4) = connec1(i, 6)

            connec(((i-1)*12)+11, 1) = connec1(i, 1)
            connec(((i-1)*12)+11, 2) = connec1(i, 2)
            connec(((i-1)*12)+11, 3) = connec1(i, 3)
            connec(((i-1)*12)+11, 4) = connec1(i, 6)

            connec(((i-1)*12)+12, 1) = connec1(i, 2)
            connec(((i-1)*12)+12, 2) = connec1(i, 3)
            connec(((i-1)*12)+12, 3) = connec1(i, 4)
            connec(((i-1)*12)+12, 4) = connec1(i, 5)
        end do

        nse = 144
        nnose = 4

    else if (elp .eq. 'PE6') then

!           12 tétraèdres dans un pentahèdre

        connec(1, 1) = 1
        connec(1, 2) = 2
        connec(1, 3) = 4
        connec(1, 4) = 5

        connec(2, 1) = 1
        connec(2, 2) = 3
        connec(2, 3) = 4
        connec(2, 4) = 5

        connec(3, 1) = 3
        connec(3, 2) = 4
        connec(3, 3) = 5
        connec(3, 4) = 6

        connec(4, 1) = 2
        connec(4, 2) = 4
        connec(4, 3) = 5
        connec(4, 4) = 6

        connec(5, 1) = 2
        connec(5, 2) = 3
        connec(5, 3) = 4
        connec(5, 4) = 6

        connec(6, 1) = 1
        connec(6, 2) = 2
        connec(6, 3) = 3
        connec(6, 4) = 5

        connec(7, 1) = 1
        connec(7, 2) = 2
        connec(7, 3) = 5
        connec(7, 4) = 6

        connec(8, 1) = 1
        connec(8, 2) = 3
        connec(8, 3) = 5
        connec(8, 4) = 6

        connec(9, 1) = 1
        connec(9, 2) = 3
        connec(9, 3) = 4
        connec(9, 4) = 6

        connec(10, 1) = 1
        connec(10, 2) = 2
        connec(10, 3) = 4
        connec(10, 4) = 6

        connec(11, 1) = 1
        connec(11, 2) = 2
        connec(11, 3) = 3
        connec(11, 4) = 6

        connec(12, 1) = 2
        connec(12, 2) = 3
        connec(12, 3) = 4
        connec(12, 4) = 5

        nse = 12
        nnose = 4
    else if (elp .eq. 'PY5') then

!       4 tétraèdres dasn une pyramide PY5

        connec(1, 1) = 1
        connec(1, 2) = 2
        connec(1, 3) = 3
        connec(1, 4) = 5

        connec(2, 1) = 1
        connec(2, 2) = 3
        connec(2, 3) = 4
        connec(2, 4) = 5

        connec(3, 1) = 1
        connec(3, 2) = 2
        connec(3, 3) = 4
        connec(3, 4) = 5

        connec(4, 1) = 2
        connec(4, 2) = 3
        connec(4, 3) = 4
        connec(4, 4) = 5

        nse = 4
        nnose = 4
    else if (elp .eq. 'TE4') then
        connec(1, 1) = 1
        connec(1, 2) = 2
        connec(1, 3) = 3
        connec(1, 4) = 4

        nse = 1
        nnose = 4
    else if (elp .eq. 'QU4') then
!       4 triangles
        connec(1, 1) = 1
        connec(1, 2) = 2
        connec(1, 3) = 4

        connec(2, 1) = 2
        connec(2, 2) = 3
        connec(2, 3) = 4

        connec(3, 1) = 1
        connec(3, 2) = 2
        connec(3, 3) = 3

        connec(4, 1) = 3
        connec(4, 2) = 1
        connec(4, 3) = 4

        nse = 4
        nnose = 3
    else if (elp .eq. 'TR3') then
        connec(1, 1) = 1
        connec(1, 2) = 2
        connec(1, 3) = 3

        nse = 1
        nnose = 3

!       My added quadratic triangle !!!!!!!!!!!!!!!!!
    else if (elp .eq. 'TR6') then

        connec(1, 1) = 1
        connec(1, 2) = 4
        connec(1, 3) = 6

        connec(2, 1) = 4
        connec(2, 2) = 2
        connec(2, 3) = 5

        connec(3, 1) = 6
        connec(3, 2) = 5
        connec(3, 3) = 3

        connec(4, 1) = 4
        connec(4, 2) = 5
        connec(4, 3) = 6

        nse = 4
        nnose = 3

!        My added quadratic quad !!!!!!!!!!!!!!!!!!!!!
    else if (elp .eq. 'QU8') then

        connec(1, 1) = 1
        connec(1, 2) = 5
        connec(1, 3) = 8

        connec(2, 1) = 5
        connec(2, 2) = 2
        connec(2, 3) = 6

        connec(3, 1) = 6
        connec(3, 2) = 3
        connec(3, 3) = 7

        connec(4, 1) = 7
        connec(4, 2) = 4
        connec(4, 3) = 8

        connec(5, 1) = 5
        connec(5, 2) = 6
        connec(5, 3) = 8

        connec(6, 1) = 6
        connec(6, 2) = 7
        connec(6, 3) = 5

        connec(7, 1) = 6
        connec(7, 2) = 7
        connec(7, 3) = 8

        connec(8, 1) = 7
        connec(8, 2) = 8
        connec(8, 3) = 5

        nse = 8
        nnose = 3
!        My added quadratic tetrahedra !!!!!!!!!!!!!!!!!!!!!
    else if (elp .eq. 'TE10') then

        connec(1, 1) = 1
        connec(1, 2) = 5
        connec(1, 3) = 7
        connec(1, 4) = 8

        connec(2, 1) = 5
        connec(2, 2) = 2
        connec(2, 3) = 6
        connec(2, 4) = 9

        connec(3, 1) = 8
        connec(3, 2) = 9
        connec(3, 3) = 10
        connec(3, 4) = 4

        connec(4, 1) = 7
        connec(4, 2) = 5
        connec(4, 3) = 3
        connec(4, 4) = 8

        connec(5, 1) = 5
        connec(5, 2) = 6
        connec(5, 3) = 3
        connec(5, 4) = 9

        connec(6, 1) = 8
        connec(6, 2) = 9
        connec(6, 3) = 10
        connec(6, 4) = 3

        connec(7, 1) = 5
        connec(7, 2) = 9
        connec(7, 3) = 8
        connec(7, 4) = 3

        nse = 7
        nnose = 4

    else
!       TYPE D'ELEMENT FINI PAS TRAITE
        ASSERT(.false.)

    end if
!
    do ise = 1, nse
        do ino = 1, nnose
            cnset(nnose*(ise-1)+ino) = connec(ise, ino)
        end do
    end do
!
    call jedema()
end subroutine
