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
subroutine dismzc(questi, nomobz, repi, repkz, ierd)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi, nomobz, repkz
!     --     DISMOI( 'Z_CST', MODELE, ... )
!    IN:
!       QUESTI : 'Z_CST'
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
! ----------------------------------------------------------------------
!
    character(len=8) :: ma, typma
    character(len=19) :: nolig
    character(len=24) :: nema
    character(len=32) :: repk
    integer(kind=8) :: idnema, ier
    integer(kind=8) :: ii, ilmaco, ima, ino, iocc, itypm
    integer(kind=8) :: jima, jnbno, nbma, nbnoma, nbnot, nbpt, numail
    integer(kind=8) :: nunoel, nunota, nutioc
    real(kind=8) :: z1
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    character(len=8), pointer :: typema(:) => null()
    integer(kind=8), pointer :: nbno(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    real(kind=8), pointer :: vale(:) => null()
! -----  FONCTIONS FORMULES
!     NUMGLM(IMA,INO)=NUMERO GLOBAL DU NOEUD INO DE LA MAILLE IMA
!                     IMA ETANT UNE MAILLE DU MAILLAGE.
#define numglm(numail,ino) connex(zi(ilmaco+numail-1)+ino-1)
! --------------------------------------------------------------------
    call jemarq()
    ASSERT(questi .eq. 'Z_CST')
!
    nolig = nomobz
    repk = ' '
    repi = 0
    ierd = 0
!
! --- LE MODELE
!
    call jeveuo(nolig//'.LGRF', 'L', vk8=lgrf)
    call jelira(nolig//'.LIEL', 'NUTIOC', nutioc)
    nema = nolig//'.NEMA'
    call jeexin(nema, ier)
!
! --- LE MAILLAGE
!
    ma = lgrf(1)
    call jeveuo(ma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilmaco)
    call jeveuo(ma//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo(ma//'.DIME', 'L', vi=dime)
    nbnoma = dime(1)
!
! --- ON CREE UN TABLEAU DONT LA COMPOSANTE I VAUDRA 1 SI LE NOEUD I
!     APPARTIENT AU MODELE
!
    call wkvect('&&DISMZC.TRAV.NOEUDS', 'V V I', nbnoma, jnbno)
!
    call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typema)
    call jeveuo('&CATA.TM.NBNO', 'L', vi=nbno)
!
    do iocc = 1, nutioc
!
        call jelira(jexnum(nolig//'.LIEL', iocc), 'LONMAX', nbma)
        call jeveuo(jexnum(nolig//'.LIEL', iocc), 'L', jima)
        typma = typema(zi(jima+nbma-1))
        call jenonu(jexnom('&CATA.TM.NOMTM', typma), itypm)
        nbpt = nbno(itypm)
        nbma = nbma-1
!
        do ii = 1, nbma
!
            numail = zi(jima+ii-1)
!
            if (numail .lt. 0) then
! --------- MAILLE TARDIVE: ON RECUPERE NEMA
                if (ier .eq. 0) then
                    call utmess('F', 'UTILITAI_71')
                end if
                ima = -numail
                call jeveuo(jexnum(nema, ima), 'L', idnema)
                call jelira(jexnum(nema, ima), 'LONMAX', nbnot)
! --------- NEMA NOUS DONNE DIRECTEMENT LE NUMERO DU NOEUD
                nbnot = nbnot-1
                do ino = 1, nbnot
                    nunota = zi(idnema+ino-1)
                    if (nunota .lt. 0) then
                        call utmess('F', 'UTILITAI_72')
                    end if
                    zi(jnbno+nunota-1) = 1
                end do
            else
!
! --------- RECUPERATION DU NOMBRE DE NOEUDS ET DE LA LISTE
!           DES NOEUDS DE LA MAILLE NUMAIL
!
                do ino = 1, nbpt
                    nunoel = numglm(numail, ino)
                    zi(jnbno+nunoel-1) = 1
                end do
            end if
        end do
    end do
!
! --- ON RECUPERE LA COORDONNEE Z DU PREMIER NOEUD DU MAILLAGE
!     CONTENU DANS LE MODELE POUR TESTER LES Z SUIVANTS
!
    do ino = 1, nbnoma
        if (zi(jnbno+ino-1) .ne. 0) then
            z1 = vale(3*(ino-1)+3)
            goto 26
        end if
    end do
26  continue
!
    do ino = 1, nbnoma
        if (zi(jnbno+ino-1) .ne. 0) then
            if (vale(3*(ino-1)+3) .ne. z1) goto 30
        end if
    end do
    repk = 'OUI'
    goto 999
30  continue
    repk = 'NON'
!
999 continue
    repkz = repk
    call jedetr('&&DISMZC.TRAV.NOEUDS')
    call jedema()
end subroutine
