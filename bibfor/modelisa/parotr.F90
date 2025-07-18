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
subroutine parotr(nomma, iageom, ima, nbno, o, &
                  mrot, t, coor)
    implicit none
#include "jeveux.h"
#include "asterfort/pacoor.h"
    character(len=8) :: nomma
    integer(kind=8) :: ima, nbno
    real(kind=8) :: o(3), mrot(3, 3), t(3), coor(*)
!     BUT: DONNER LA LISTE DES COORDONNEES DES NBNO 1ERS NOEUDS DE LA
!          MAILLE IMA DU MAILLAGE NOMMA APRES UNE ROTATION (O,MROT)
!          ET UNE TRANSLATION T OU DU NOEUD IMA SI NBNO=0
!     VERIFICTION : NBNO < OU = NBRE DE NOUDS DE LA MAILLE
! ARGUMENTS D'ENTREE:
! IN   NOMMA  K8  : NOM DU MAILLAGE
! IN   IAGEOM I   : ADRESSE DE L'OBJET MAILLAGE//'.COORDO   .VALE'
! IN   IMA    I   : NUMERO DE LA MAILLE OU DU NOEUD SI NBNO=0
! IN   NBNO   I   : NOMBRE DE NOEUDS DE LA MAILLE A EXTRAIRE OU 0
! IN   O      R(3): CENTRE DE ROTATION
! IN   MROT   R(9): MATRICE DE ROTATION
! IN   T      R(3): VECTEUR DE TRANSLATION
! OUT  COOR   R(*): COORDONNEES DES NBNO 1ERS NOEUDS DE LA MAILLE
!                   APRES TRANSFORMATION OU DU NOEUD
!                   POUR INO = 1,NBNO OU INO = IMA SI NBNO=0
!                   COOR(3*(INO-1)+1)= X1'(INO)
!                   COOR(3*(INO-1)+2)= X2'(INO)
!                   COOR(3*(INO-1)+3)= X3'(INO) ( EN 2D 0)
    real(kind=8) :: x(3)
! --- DEBUT
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iageom, icoor, ino, j, nbnot
!-----------------------------------------------------------------------
    if (nbno .eq. 0) then
        nbnot = 1
        coor(1) = zr(iageom-1+3*(ima-1)+1)
        coor(2) = zr(iageom-1+3*(ima-1)+2)
        coor(3) = zr(iageom-1+3*(ima-1)+3)
    else
        nbnot = nbno
        call pacoor(nomma, ima, nbno, coor)
    end if
!
    do ino = 1, nbnot
        icoor = 3*(ino-1)
        do i = 1, 3
            x(i) = 0.d0
            do j = 1, 3
                x(i) = x(i)+mrot(j, i)*(coor(icoor+j)-o(j))
            end do
        end do
        do i = 1, 3
            coor(icoor+i) = x(i)+o(i)+t(i)
        end do
    end do
end subroutine
