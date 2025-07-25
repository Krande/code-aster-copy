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
subroutine resmod(bmodal, nbmode, neq, numgen, mdgene, &
                  noecho, modsst)
    implicit none
!  C. VARE     DATE 16/10/95
!-----------------------------------------------------------------------
!  BUT : < CALCUL DEPLACEMENT MODAL >
!  CALCULER LES DEPLACEMENTS MODAUX D'UN NOEUD D'UNE SOUS-STRUCTURE
!-----------------------------------------------------------------------
!
! BMODAL /I/ : BASE MODALE DE LA STRUCTURE COMPLETE
! NBMODE /I/ : NOMBRE DE MODES DE LA STRUCTURE COMPLETE
! NEQ    /I/ : NOMBRE D'EQUATIONS DU MODELE GENERALISE
! NUMGEN /I/ : NUMEROTATION DU PROBLEME GENERALISE
! MDGENE /I/ : MODELE GENERALISE
! NOECHO /I/ : NOEUD A RESTITUER : NOECHO(1) = NOEUD_1
!                                  NOECHO(2) = SOUS_STRUC_1
!                                  NOECHO(3) = NUME_1
! MODSST /O/ : DEPL PHYSIQUES DES MODES AUX NOEUDS DE CHOC
!
!
!
#include "jeveux.h"
#include "asterfort/dcapno.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/orient.h"
#include "asterfort/posddl.h"
#include "asterfort/wkvect.h"
!
!
!
    integer(kind=8) :: nbmode, ddl(6), neq
    character(len=8) :: basmod, nomsst, soutr, numddl, noeud, noecho(3)
    character(len=16) :: depl
    character(len=24) :: chamba, mdgene, numgen
    real(kind=8) :: bmodal(neq, *), coord(3), modsst(nbmode, 6)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, ieq, j, jcoord, lchab, llors
    integer(kind=8) :: llprs, neqgen, nsst, nunoe, nusst, nutars
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data soutr/'&SOUSSTR'/
!-----------------------------------------------------------------------
!
    call jemarq()
    noeud = noecho(1)
    nomsst = noecho(2)
    numddl = noecho(3)
    numgen = numgen(1:14)//'.NUME'
!
! --- NUMERO DE SOUS-STRUCTURE ET DU NOEUD TARDIF CORRESPONDANT
!
    call jenonu(jexnom(mdgene(1:14)//'.MODG.SSNO', nomsst), nusst)
    call jenonu(jexnom(numgen(1:19)//'.LILI', soutr), ibid)
    call jeveuo(jexnum(numgen(1:19)//'.ORIG', ibid), 'L', llors)
    call jelira(jexnum(numgen(1:19)//'.ORIG', ibid), 'LONMAX', nsst)
    do i = 1, nsst
        if (zi(llors+i-1) .eq. nusst) nutars = i
    end do
!
    call jenonu(jexnom(numgen(1:19)//'.LILI', soutr), ibid)
    call jeveuo(jexnum(numgen(1:19)//'.PRNO', ibid), 'L', llprs)
    neqgen = zi(llprs+(nutars-1)*2+1)
    ieq = zi(llprs+(nutars-1)*2)
!
! --- BASE MODALE ET DU NBRE D'EQUATIONS DE LA SOUS-STRUCTURE
!
    call mgutdm(mdgene, nomsst, ibid, 'NOM_BASE_MODALE', ibid, &
                basmod)
!
! --- RESTITUTION PROPREMENT DITE
!
    call posddl('NUME_DDL', numddl, noeud, 'DX', nunoe, &
                ddl(1))
    call posddl('NUME_DDL', numddl, noeud, 'DY', nunoe, &
                ddl(2))
    call posddl('NUME_DDL', numddl, noeud, 'DZ', nunoe, &
                ddl(3))
    call posddl('NUME_DDL', numddl, noeud, 'DRX', nunoe, &
                ddl(4))
    call posddl('NUME_DDL', numddl, noeud, 'DRY', nunoe, &
                ddl(5))
    call posddl('NUME_DDL', numddl, noeud, 'DRZ', nunoe, &
                ddl(6))
!
    call wkvect('&&RESMOD.COORDO', 'V V R', 3, jcoord)
    do i = 1, nbmode
        modsst(i, 1) = 0.d0
        modsst(i, 2) = 0.d0
        modsst(i, 3) = 0.d0
        modsst(i, 4) = 0.d0
        modsst(i, 5) = 0.d0
        modsst(i, 6) = 0.d0
        zr(jcoord) = 0.d0
        zr(jcoord+1) = 0.d0
        zr(jcoord+2) = 0.d0
        do j = 1, neqgen
            call dcapno(basmod, depl, j, chamba)
            call jeveuo(chamba, 'E', lchab)
            modsst(i, 1) = modsst(i, 1)+bmodal(ieq+j-1, i)*zr(lchab+ddl(1)- &
                                                              1)
            modsst(i, 2) = modsst(i, 2)+bmodal(ieq+j-1, i)*zr(lchab+ddl(2)- &
                                                              1)
            modsst(i, 3) = modsst(i, 3)+bmodal(ieq+j-1, i)*zr(lchab+ddl(3)- &
                                                              1)
            if (ddl(4) .ne. 0 .and. ddl(5) .ne. 0 .and. ddl(6) .ne. 0) then
                modsst(i, 4) = modsst(i, 4)+bmodal(ieq+j-1, i)*zr(lchab+ &
                                                                  ddl(4)-1)
                modsst(i, 5) = modsst(i, 5)+bmodal(ieq+j-1, i)*zr(lchab+ &
                                                                  ddl(5)-1)
                modsst(i, 6) = modsst(i, 6)+bmodal(ieq+j-1, i)*zr(lchab+ &
                                                                  ddl(6)-1)
            end if
        end do
        zr(jcoord) = modsst(i, 1)
        zr(jcoord+1) = modsst(i, 2)
        zr(jcoord+2) = modsst(i, 3)
        call orient(mdgene, nomsst, jcoord, 1, coord, &
                    0)
        modsst(i, 1) = coord(1)
        modsst(i, 2) = coord(2)
        modsst(i, 3) = coord(3)
        zr(jcoord) = modsst(i, 4)
        zr(jcoord+1) = modsst(i, 5)
        zr(jcoord+2) = modsst(i, 6)
        call orient(mdgene, nomsst, jcoord, 1, coord, &
                    0)
        modsst(i, 4) = coord(1)
        modsst(i, 5) = coord(2)
        modsst(i, 6) = coord(3)
    end do
    call jedetr('&&RESMOD.COORDO')
!
    call jedema()
end subroutine
