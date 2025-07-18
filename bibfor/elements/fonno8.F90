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

subroutine fonno8(resu, noma, tablev, vect)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/panbno.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: noma, resu
    integer(kind=8) :: tablev(2)
    real(kind=8) :: vect(3)
!
!      DETERMINATION D'UN VECTEUR SE DIRIGEANT VERS LA LEVRE SUPERIEURE
!       ----------------------------------------------------------------
!    ENTREES
!       RESU   : NOM DU CONCEPT RESULTAT DE L'OPERATEUR
!       NOMA   : NOM DU MAILLAGE
!       TABLEV : VECTEUR CONTENANT LES NUMEROS DES DEUX MAILLES
!                CONNECTEES AU NOEUD SOMMET COURANT DU PREMIER SEGMENT
!                ET AUX LEVRES
!    SORTIE
!       VECT   : VECTEUR SE DIRIGEANT VERS LA LEVRE SUPERIEURE
!                IL SERT A REORIENTER LE VECTEUR NORMAL POUR QU'IL AILLE
!                DE LA LEVRE INF VERS LA LEVRE SUP
!
!
    integer(kind=8) :: comp
    integer(kind=8) :: iamase, iatyma, ifon, ilev, inn, inn2, inp, iret, ityp, itypma
    integer(kind=8) ::  jconx, jfon
    integer(kind=8) :: nblev, nn, nn2, nbnott(3)
    real(kind=8) :: xg, yg, zg
    character(len=8) :: type
    character(len=8), pointer :: mail(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
!     RECUPERATION DES DONNES SUR LE MAILLAGE
    call jeveuo(noma//'.TYPMAIL', 'L', iatyma)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
    call jeveuo(resu//'.LEVRESUP.MAIL', 'L', vk8=mail)
    call jelira(resu//'.LEVRESUP.MAIL', 'LONUTI', nblev)
!
!
!     DETERMINATION D'UN VECTEUR ALLANT DU PREMIER NOEUD DU FOND
!     FISSURE A UN NOEUD DE LA LEVRE_SUP
!
    do inp = 1, 2
!
        call jeveuo(jexnum(noma//'.CONNEX', tablev(inp)), 'L', jconx)
        itypma = iatyma-1+tablev(inp)
!
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(itypma)), type)
        call dismoi('NBNO_TYPMAIL', type, 'TYPE_MAILLE', repi=nn2)
!
        do ilev = 1, nblev
            comp = 0
            iret = char8_to_int(mail(ilev))
            call jeveuo(jexnum(noma//'.CONNEX', iret), 'L', iamase)
            ityp = iatyma-1+iret
!
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
            call dismoi('NBNO_TYPMAIL', type, 'TYPE_MAILLE', repi=nn)
!
            do inn = 1, nn
                do inn2 = 1, nn2
                    if (zi(jconx-1+inn2) .eq. zi(iamase-1+inn)) then
                        comp = comp+1
                        if (comp .eq. nn) goto 300
                    end if
                end do
            end do
        end do
    end do
!
300 continue
!
!     CALCUL DES COORDONNEES DU CENTRE DE GRAVITE
    call panbno(zi(itypma), nbnott)
!
    xg = 0
    yg = 0
    zg = 0
    do inn2 = 1, nbnott(1)
        xg = xg+vale((zi(jconx-1+inn2)-1)*3+1)
        yg = yg+vale((zi(jconx-1+inn2)-1)*3+2)
        zg = zg+vale((zi(jconx-1+inn2)-1)*3+3)
    end do
!
    call jeexin(resu//'.FOND.NOEU', ifon)
    if (ifon .ne. 0) then
        call jeveuo(resu//'.FOND.NOEU', 'L', jfon)
    else
        ASSERT(.FALSE.)
    end if
    iret = char8_to_int(zk8(jfon))
!
    vect(1) = xg/nbnott(1)-vale((iret-1)*3+1)
    vect(2) = yg/nbnott(1)-vale((iret-1)*3+2)
    vect(3) = zg/nbnott(1)-vale((iret-1)*3+3)
!
    call jedema()
end subroutine
