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
subroutine zreord(zmat, nbddg, nbmod, nbmob, nbddr, &
                  axok, liax, nbliax, zvec)
    implicit none
#include "asterf_types.h"
!
!  BUT:  < REMISE EN ORDRE DES COORDONNEES COMPLEXES >
!
!   CETTE ROUTINE MET EN ORDRE LES COORDONNEES COMPLEXES DES MODE
!  DANS LE CAS D'ASSEMBLAGE PARTIEL DU A  UNE INTEFACE AXE EN
!  CYCLIQUE CRAIG-BAMPTON
!-----------------------------------------------------------------------
!
! ZMAT      /M/: MATRICE DES COORDONNEES GENERALISEES DES MODES PROPRES
! NBDDG     /I/: NOMBRE DE DEGRE DE LIBETE GENERALISE TOTAL
! NBMOD     /I/: NOMBRE DE MODE CYCLIQUES CALCULES
! NBMOB     /I/: NOMBRE DE MODE PROPRES DE LA BASE MODALE
! NBDDR     /I/: NOMBRE DDL GENERALISE DE DROITE (=GAUCHE)
! AXOK      /I/: INDICATEUR ASSEMBLAGE AXE
! LIAX      /I/: LISTE NUMERO DDL AXE A ASSEMBLES
! NBLIAX    /I/: NOMBRE DE DDL AXE A ASSEMBLES
! ZVEC      /M/: VECTEUR DE TRAVAIL COMPLEXE DIMENSIONNE A NBDDG
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, nbddg, nbddr, nbliax, nbmob, nbmod
    integer(kind=8) :: liax(nbliax)
    complex(kind=8) :: zmat(nbddg, nbmod), zvec(nbddg)
    aster_logical :: axok
!-----------------------------------------------------------------------
!
!  CAS DE PRESENCE DDL AXE
!
!-----------------------------------------------------------------------
    if (axok) then
        do j = 1, nbmod
!
            do i = 1, nbliax
                zvec(i) = zmat(i+nbmob+nbddr, j)
            end do
!
            if ((nbmob+nbddr) .lt. nbddg) then
                do i = nbmob+nbddr+1, nbddg
                    zmat(i, j) = dcmplx(0.d0, 0.d0)
                end do
            end if
!
            do i = 1, nbliax
                zmat(nbmob+nbddr+liax(i), j) = zvec(i)
            end do
!
        end do
!
!  AUTRE CAS
!
    else
        if ((nbmob+nbddr) .lt. nbddg) then
            do j = 1, nbmod
                do i = nbmob+nbddr+1, nbddg
                    zmat(i, j) = dcmplx(0.d0, 0.d0)
                end do
            end do
        end if
!
    end if
end subroutine
