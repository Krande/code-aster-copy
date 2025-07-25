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
subroutine caexno(lvavz, nomaz, motfac, mcgrno, mcno, &
                  iocc)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/palino.h"
!
    character(len=24) :: lvav, noma
    character(len=*) :: motfac, mcgrno, mcno, lvavz, nomaz
    integer(kind=8) :: iocc
!                 DATE 04/01/93   C.MENGONI
!
!     BUT : RETIRER DE LA LISTE DE VIS A VIS LVAV TOUS LES COUPLES
!           DONT UN DES TERMES EST APPARTIENT A LA LISTE DE NOEUDS
!           GENEREE PAR LES MOTS CLES MCGRNO ET MCNO POUR L'OCCURENCE
!           IOCC DU MOT CLE FACTEUR MOTFAC.
! IN  LVAVZ  K*(*): NOM UTILISATEUR DE LA LISTE DES VIS A VIS
!            OJB V V I DIM = 2 * NBCOUPLE + 1
!            LVAV(ILAD) = NBCOUPLE (BOUCLE SUR I <= NBCOUPLE)
!            LVAV(ILAD+2*(I-1)+1) = NUM1 DU NOEUD 1
!            LVAV(ILAD+2*(I-1)+2) = NUM2 DU NOEUD 2
!                   OU PREFIXE DES OJB .CONI ET .CONR
! IN  NOMAZ  K*(*): NOM DU MAILLAGE
! IN  MOTFAC K16  : MOT CLE FACTEUR A TRAITER
! IN  MCGRNO K16  : MOT CLE REGROUPANT LES GROUP_NO
! IN  MCNO   K16  : MOT CLE REGROUPANT LES NOEUDS
! IN  IOCC   I    : SI >0 ON TRAITE L'OCCURENCE IOCC DE MOTFAC
!                   SI <0 OU =0 ERREUR FATALE
!
!
!
    character(len=16) :: motf, mcgr, mcn
    character(len=24) :: listex
    integer(kind=8) :: nbcpl, nbex
! --- DEBUT
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idlex, idlvav, j, l, nlino
!-----------------------------------------------------------------------
    call jemarq()
    lvav = lvavz
    noma = nomaz
    motf = motfac
    mcgr = mcgrno
    mcn = mcno
    if (motf .ne. 'LIAISON_GROUP') then
        ASSERT(.false.)
    end if
    if (iocc .le. 0) then
        ASSERT(.false.)
    end if
    call getfac(motf, nlino)
    if ((nlino .eq. 0) .or. (iocc .gt. nlino)) goto 999
!
! --- LECTURE DE LA LISTE DES NOEUDS EXCLUS
!
    listex = '&&CAEXNO.LISTENOEUD'
    call palino(noma, motfac, mcgr, mcn, iocc, &
                listex)
!
! --- ELIMINATION
!
    call jeveuo(listex, 'L', idlex)
    nbex = zi(idlex)
    if (nbex .eq. 0) goto 998
    call jeveuo(jexnum(lvav, iocc), 'E', idlvav)
    nbcpl = zi(idlvav)
    l = 0
    do i = 1, nbcpl
        l = l+1
        do j = 1, nbex
            if ((zi(idlvav+2*(i-1)+1) .eq. zi(idlex+j)) .or. &
                (zi(idlvav+2*(i-1)+2) .eq. zi(idlex+j))) then
                l = l-1
                goto 2
            end if
        end do
        zi(idlvav+2*(l-1)+1) = zi(idlvav+2*(i-1)+1)
        zi(idlvav+2*(l-1)+2) = zi(idlvav+2*(i-1)+2)
2       continue
    end do
    zi(idlvav) = l
998 continue
    call jedetr('&&CAEXNO.LISTENOEUD')
999 continue
! FIN -----------------------------------------------------------------
    call jedema()
end subroutine
