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

subroutine mtcopy(matin, matout, ier)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdscr.h"
#include "asterfort/utmess.h"
#include "asterfort/vrrefe.h"
    character(len=*) :: matin, matout
    integer(kind=8) :: ier
!     RECOPIE LES VALEURS DE LA MATRICE MATIN  DANS LA MATRICE MATOUT
!     ------------------------------------------------------------------
!     PRECAUTION D'EMPLOI :
!        1) LA MATRICE "MATOUT" DOIT EXISTER ET AVOIR LA MEME STRUCTURE
!     QUE "MATIN"
!        2) ON RECOPIE LE .CCID DE MATIN DANS
!     MATOUT, SI MATOUT POSSEDAIT DEJA CE CHAMP ON LE DETRUIT.
!     ------------------------------------------------------------------
!     RAPPEL :   UNE MATRICE  "MAT" EXISTE
!          S'IL EXISTE UN OBJET SIMPLE  MAT//"REFE"
!          ET UNE COLLECTION NUMEROTEE  MAT//"VALE"
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: lmatou, lmatin, nimpou
    character(len=8) :: nomddl
    character(len=19) :: mati19, mato19
    character(len=24) :: nmatou, nmatin
    character(len=24) :: valk(2)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    data nomddl/'        '/
!     ------------------------------------------------------------------
!
!     --- CONTROLE DES REFERENCES ---
    call jemarq()
    call vrrefe(matin, matout, ier)
    mati19 = matin
    mato19 = matout
    if (ier .ne. 0) then
        valk(1) = mati19
        valk(2) = mato19
        call utmess('F', 'ALGELINE2_11', nk=2, valk=valk)
!
    else
!        --- TYPE DES VALEURS, NOMBRE DE BLOCS, LONGUEUR D'UN BLOC ---
        call mtdscr(matin)
        nmatin = matin(1:19)//'.&INT'
        call jeveuo(matin(1:19)//'.&INT', 'E', lmatin)
        call mtdscr(matout)
        nmatou = matout(1:19)//'.&INT'
        call jeveuo(matout(1:19)//'.&INT', 'E', lmatou)
!
! --- GESTION DES .CCID .CCLL .CCVA
!
        nimpou = zi(lmatou+7)
        if (nimpou .ne. 0) then
            call jedetr(mato19//'.CCID')
            call jedetr(mato19//'.CCLL')
            call jedetr(mato19//'.CCII')
            call jedetr(mato19//'.CCVA')
            zi(lmatou+7) = 0
            zi(lmatou+15) = 0
            zi(lmatou+16) = 0
        end if
!
!
! --- RECOPIE DU .VALE ET DE .CCID, .CCLL, .CCVA
        call mtcmbl(1, 'R', [1.d0], nmatin, nmatou, &
                    nomddl, ' ', 'ELIM=')
    end if
!
    call jedema()
end subroutine
