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
subroutine pacoor(nomma, ima, nbno, coor)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: nomma
    integer(kind=8) :: ima, nbno
    real(kind=8) :: coor(*)
!     BUT: DONNER LA LISTE DES COORDONNEES DES NBNO 1ERS NOEUDS DE LA
!          MAILLE IMA DU MAILLAGE NOMMA OU D'UN NOEUD SI NBNO = 0
!     VERIFICTION : NBNO < OU = NBRE DE NOUDS DE LA MAILLE
! ATTENTION IL FAUT QUE DIM DE COOR >= 3 MEME POUR UN NOEUD EN 2D
!---------------------------------------------------------------------
! ARGUMENTS D'ENTREE:
! IN   NOMMA  K8  : NOM DU MAILLAGE
! IN   IMA    I   : NUMERO DE LA MAILLE OU DU NOEUD SI NBNO = 0
! IN   NBNO   I   : NOMBRE DE NOEUDS DE LA MAILLE A EXTRAIRE, OU 0 POUR
!                   UN NOEUD
! OUT  COOR   R(*): COORDONNEES DES NBNO 1ERS NOEUDS DE LA MAILLE
!                   OU COORDONNEES DU NOEUD IMA
!                   POUR INO = 1,NBNO  OU INO = IMA SI NBNO = 0
!                   COOR(3*(INO-1)+1)= X1(INO)
!                   COOR(3*(INO-1)+2)= X2(INO)
!                   COOR(3*(INO-1)+3)= X3(INO) ( EN 2D 0)
    character(len=24) :: desc, vale, connex
    real(kind=8) :: x(3)
! --- DEBUT
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icmp, icoor, idconn, iddesc, idino, idvale
    integer(kind=8) :: ino, inoma, nbcmp, nbnomx
!-----------------------------------------------------------------------
    call jemarq()
    desc = nomma(1:8)//'.COORDO    .DESC'
    vale = nomma(1:8)//'.COORDO    .VALE'
    connex = nomma(1:8)//'.CONNEX'
    call jeveuo(desc, 'L', iddesc)
    nbcmp = -zi(iddesc+1)
    if (nbcmp .eq. 2) then
        if (nbno .gt. 0) x(3) = 0.d0
        if (nbno .eq. 0) coor(3) = 0.d0
    end if
    call jeveuo(vale, 'L', idvale)
    if (nbno .gt. 0) then
        call jeveuo(jexnum(connex, ima), 'L', idconn)
        call jelira(jexnum(connex, ima), 'LONMAX', nbnomx)
        if (nbno .gt. nbnomx) then
            call utmess('F', 'MODELISA6_5')
        end if
        do inoma = 1, nbno
            ino = zi(idconn-1+inoma)
            idino = idvale+nbcmp*(ino-1)-1
            do icmp = 1, nbcmp
                x(icmp) = zr(idino+icmp)
            end do
            icoor = 3*(inoma-1)
            do i = 1, 3
                coor(icoor+i) = x(i)
            end do
        end do
    else if (nbno .eq. 0) then
        idino = idvale+nbcmp*(ima-1)-1
        do icmp = 1, nbcmp
            coor(icmp) = zr(idino+icmp)
        end do
    else
        call utmess('F', 'MODELISA6_6')
    end if
    call jedema()
end subroutine
