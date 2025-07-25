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
subroutine crpcvg(ma1, ma2, gma1, gma2, tran, &
                  prec, lima1, lima2, linoeu)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    real(kind=8) :: tran(3), prec
    integer(kind=8) :: linoeu(*)
    character(len=8) :: ma1, ma2
    character(len=24) :: gma1, gma2
    character(len=*) :: lima1, lima2
!
!     COMMANDE:  CREA_RESU
!     TRAITEMENT DU MOT CLE FACTEUR "PERM_CHAMP"
!
! ----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: nbma1, nbma2, jtyma1, jtyma2, jgma1, jgma2, ima, ima1, ima2, ino
    integer(kind=8) :: ino1, ino2, nutyp1, iamac1, ilmac1, nutyp2
    integer(kind=8) :: iamac2, ilmac2, jcoor1, jcoor2, jnum1, jnum2, jma
    integer(kind=8) :: ibid
    real(kind=8) :: x1, y1, z1, x2, y2, z2, v1, v2, v3
    real(kind=8) :: valr(3)
    aster_logical :: erreur
    character(len=8) :: noma1
    character(len=24) :: valk(5), nom_gr1, nom_gr2
    character(len=24) :: grpma1, grpma2, coova1, coova2, typma1, typma2, conne1
    character(len=24) :: conne2
!
!     NBNOMA(IMA) = NOMBRE DE NOEUDS DE LA MAILLE IMA
#define nbnoma(ima) zi(ilmac1-1+ima+1) - zi(ilmac1-1+ima)
!
!     NUMGLM(IMA,INO) = NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
#define numgl1(ima,ino) zi(iamac1-1+zi(ilmac1+ima-1)+ino-1)
#define numgl2(ima,ino) zi(iamac2-1+zi(ilmac2+ima-1)+ino-1)
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
    grpma1 = ma1//'.GROUPEMA       '
    grpma2 = ma2//'.GROUPEMA       '
    coova1 = ma1//'.COORDO    .VALE'
    coova2 = ma2//'.COORDO    .VALE'
    typma1 = ma1//'.TYPMAIL        '
    typma2 = ma2//'.TYPMAIL        '
    conne1 = ma1//'.CONNEX         '
    conne2 = ma2//'.CONNEX         '
!
    call jeveuo(coova1, 'L', jcoor1)
    call jeveuo(coova2, 'L', jcoor2)
!
    call jeveuo(typma1, 'L', jtyma1)
    call jeveuo(typma2, 'L', jtyma2)
!
    call jeveuo(conne1, 'L', iamac1)
    call jeveuo(jexatr(conne1, 'LONCUM'), 'L', ilmac1)
    call jeveuo(conne2, 'L', iamac2)
    call jeveuo(jexatr(conne2, 'LONCUM'), 'L', ilmac2)
!
    call jelira(jexnom(grpma1, gma1), 'LONUTI', nbma1)
    call jelira(jexnom(grpma2, gma2), 'LONUTI', nbma2)
    if (nbma1 .ne. nbma2) then
        valk(1) = gma1
        valk(2) = gma2
        call utmess('F', 'CALCULEL5_67', nk=2, valk=valk)
    end if
!
    call jeveuo(jexnom(grpma1, gma1), 'L', jgma1)
    call jeveuo(jexnom(grpma2, gma2), 'L', jgma2)
!
    call wkvect(lima1, 'V V I', nbma1, jnum1)
    call wkvect(lima2, 'V V I', nbma1, jnum2)
!
!   boucle sur les mailles de GROUP_MA_INIT
    do ima = 1, nbma1
!
!       numero et type de la maille courante de GROUP_MA_INIT
        ima1 = zi(jgma1+ima-1)
        nutyp1 = zi(jtyma1-1+ima1)
!
!       boucle sur les mailles de GROUP_MA_FINAL
        do jma = 1, nbma2
!
            ima2 = zi(jgma2+jma-1)
            nutyp2 = zi(jtyma2-1+ima2)
!
            do ino = 1, nbnoma(ima1)
                ino1 = numgl1(ima1, ino)
                ino2 = numgl2(ima2, ino)
                x1 = zr(jcoor1-1+3*(ino1-1)+1)
                y1 = zr(jcoor1-1+3*(ino1-1)+2)
                z1 = zr(jcoor1-1+3*(ino1-1)+3)
                x2 = zr(jcoor2-1+3*(ino2-1)+1)
                y2 = zr(jcoor2-1+3*(ino2-1)+2)
                z2 = zr(jcoor2-1+3*(ino2-1)+3)
                v1 = abs(x2-x1-tran(1))
                v2 = abs(y2-y1-tran(2))
                v3 = abs(z2-z1-tran(3))
                erreur = .false.
                if (v1 .gt. prec) erreur = .true.
                if (v2 .gt. prec) erreur = .true.
                if (v3 .gt. prec) erreur = .true.
!
!               en cas d'erreur, on passe a la maille suivante de
!               GROUP_MA_FINAL (iteration suivante de la boucle 100)
                if (erreur) then
                    goto 100
                end if
!
                linoeu(ino2) = ino1
!
            end do
!
!           a ce niveau, on a trouve la maille ima2 en
!           correspondance avec la maille ima1
!           on verifie leur type
            if (nutyp1 .ne. nutyp2) then
                valk(1) = gma1
                valk(2) = gma2
                call utmess('F', 'CALCULEL5_68', nk=2, valk=valk)
            end if
!
!           on stocke ima1 et ima2
!           et on passe a la maille suivante dans GROUP_MA_INIT
            zi(jnum1+ima-1) = ima1
            zi(jnum2+ima-1) = ima2
            goto 10
!
100         continue
        end do
!
!       si on passe ici, c'est que l'on n'a trouve aucune maille de
!       GROUP_MA_FINAL en correspondance avec la maille ima1 de GROUP_MA_INIT
        call getvtx('PERM_CHAM', 'GROUP_MA_INIT', iocc=1, nbval=1, vect=nom_gr1, &
                    nbret=ibid)
        call getvtx('PERM_CHAM', 'GROUP_MA_FINAL', iocc=1, nbval=1, vect=nom_gr2, &
                    nbret=ibid)
        noma1 = int_to_char8(ima1)
        valk(1) = nom_gr1
        valk(2) = nom_gr2
        valk(3) = noma1
        valk(4) = ma1
        valk(5) = ma2
        valr(1) = tran(1)
        valr(2) = tran(2)
        valr(3) = tran(3)
        call utmess('F', 'CALCULEL5_69', nk=5, valk=valk, nr=3, &
                    valr=valr)
10      continue
    end do
!
    call jedema()
end subroutine
