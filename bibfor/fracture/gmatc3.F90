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

subroutine gmatc3(nnoff, milieu, connex, &
                  abscur, matr)

    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gmate3.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"

    integer(kind=8)           :: nnoff
    character(len=24) :: abscur
    character(len=24) :: matr
    aster_logical     :: milieu, connex

!       CALCUL DE LA MATRICE DU SYSTEME LINEAIRE [A] {GS} = {GTHI}
!       POUR LA METHODE THETA-LAGRANGE ET G-LAGRANGE
!
! ENTREE
!
!   NNOFF    --> NOMBRE DE NOEUDS DU FOND DE FISSURE
!   ABSCUR   --> ABSCISSES CURVILIGNES S
!   MILIEU   --> .TRUE.  : ELEMENT QUADRATIQUE
!                .FALSE. : ELEMENT LINEAIRE
!   CONNEX   --> .TRUE.  : FOND DE FISSURE FERME
!                .FALSE. : FOND DE FISSURE OUVERT

! SORTIE
!
!   MATR     --> MATRICE DU SYTEME A RESOUDRE
!
! ......................................................................

    integer(kind=8)          :: nseg, iseg, imatr
    integer(kind=8)          :: i, j, ij, nno
    integer(kind=8)          :: conn(3), conn2(3)
    real(kind=8)     :: mele(3, 3)
    character(len=8) :: elrefe

! ......................................................................
!
    call jemarq()

    conn(1:3) = 0

!   NOMBRE DE SEGMENT DU FOND DE FISSURE
    if (milieu) then
        nseg = (nnoff-1)/2
        elrefe = 'SE3'
    else
        nseg = nnoff-1
        elrefe = 'SE2'
    end if

!   CREA OBJET TEMP POUR LA VAL DE GLEGEN A ABSC CURV S
    call wkvect(matr, 'V V R8', nnoff*nnoff, imatr)
!
!  BOUCLE SUR LES SEGMENTS
    do iseg = 1, nseg

        if (milieu) then
            conn(1) = 2*iseg-1
            conn(2) = 2*iseg+1
            conn(3) = 2*iseg
        else
            conn(1) = iseg
            conn(2) = iseg+1
        end if
!
!       CALCUL DE LA MATRICE ELEMENTAIRE POUR L'ELEMENT COURANT
        call gmate3(abscur, elrefe, conn, nno, mele)

!       AJOUT DE LA CONTRIBUTION DE DE LA MATRICE DE MASSE ELEMENTAIRE
!       A LA MATRICE DE MASSE ASSEMBLEE
        do i = 1, nno
            do j = 1, nno
                ij = (conn(i)-1)*nnoff+conn(j)
                zr(imatr+ij-1) = zr(imatr+ij-1)+mele(i, j)
            end do
        end do

!       CAS CONNEX >> FISSURE FERME
        if (connex) then

!           AJOUT D'UNE CONTRIBUTION SUR LA LIGNE NDIMTE
            if (conn(1) .eq. 1) then
                conn2 = conn
                conn2(1) = nnoff

                i = 1
                do j = 1, nno
                    ij = (conn2(i)-1)*nnoff+conn2(j)
                    zr(imatr+ij-1) = zr(imatr+ij-1)+mele(i, j)
                end do
            end if

!           AJOUT D'UNE CONTRIBUTION SUR LA LIGNE 1
            if (conn(2) .eq. nnoff) then
                conn2 = conn
                conn2(2) = 1

                i = 2
                do j = 1, nno
                    ij = (conn2(i)-1)*nnoff+conn2(j)
                    zr(imatr+ij-1) = zr(imatr+ij-1)+mele(i, j)
                end do
            end if
        end if

    end do

    call jedema()

end subroutine
