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
subroutine xtyhea(nfiss, ifiss, ima, nno, jconx1, &
                  jconx2, jstnl, jstnv, nbheav)
! person_in_charge: patrick.massin at edf.fr
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: nfiss, ifiss, ima, nno, nbheav
    integer(kind=8) :: jconx1, jconx2, jstnl(nfiss), jstnv(nfiss)
! ----------------------------------------------------------------------
!
! --- ROUTINE XFEM
!
! --- TYPE D'ELEMENT D'INTERSECTION XFEM
!
! ----------------------------------------------------------------------
!
! OUT NBHEAV   : NOMBRE D'ENRICHISSEMENTS HEAVISIDES DANS L'ÉLÉMENT
!
!
!
!
    integer(kind=8) :: ino, stno(nno), ifis
    integer(kind=8) :: nngl
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATION
    nbheav = 0
    do ino = 1, nno
        stno(ino) = 0
    end do
!
! --- BOUCLE SUR LES NOEUDS DE LA MAILLE
!
    do ino = 1, nno
        nngl = zi(jconx1-1+zi(jconx2+ima-1)+ino-1)
!
! --- BOUCLE SUR LES FISSURES PRECEDENTES
!
        do ifis = 1, ifiss
! --- RECUPERATION DES STATUTS DES NOEUDS
            if (zl(jstnl(ifis)-1+nngl) .and. zi(jstnv(ifis)-1+nngl) .gt. 0) then
                stno(ino) = stno(ino)+1
            end if
            nbheav = max(nbheav, stno(ino))
        end do
    end do
    call jedema()
end subroutine
