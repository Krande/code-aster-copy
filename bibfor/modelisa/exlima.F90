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

subroutine exlima(motfaz, iocc, base, modelz, ligrel)
    implicit none
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnoms2.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
    character(len=*) :: motfaz, base, modelz, ligrel
    integer(kind=8) :: iocc
! but  :  scruter les mots cle tout/group_ma/maille pour creer
!         un ligrel "reduit" a partir du ligrel du modele modelz
!
! in  : modelz : nom du modele
!
! out/jxout   : ligrel  : ligrel reduit
!     attention :
!          - le nom de ligrel est toujours "out"
!          - parfois on rend ligrel=ligrel(modele) :
!             - alors on ne tient donc pas compte de 'base'
!             - il ne faut pas le detruire !
!          - parfois on en cree un nouveau sur la base 'base'
!             - le nom du ligrel est obtenu par gnomsd
!  -----------------------------------------------------------------
!
    integer(kind=8) :: n1, jma, nbma
    character(len=8) :: modele, noma
    character(len=16) :: motfac, motcle(2), typmcl(2), oper, k16b
    character(len=19) :: ligrmo
    character(len=24) :: lismai, noojb
!  -----------------------------------------------------------------
!
    motfac = motfaz
    modele = modelz
    if (modele .eq. ' ') then
        call utmess('F', 'UTILITAI8_10')
    end if
!
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    lismai = '&&EXLIMA.LISTE_MAILLES'
!
!
!     --  SI ON DOIT TOUT PRENDRE , LIGREL = LIGRMO
!     ------------------------------------------------------
    if (motfac .ne. ' ') then
        if (getexm(motfac, 'TOUT') .eq. 1) then
            call getvtx(motfac, 'TOUT', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) goto 9998
        end if
    else
        if (getexm(' ', 'TOUT') .eq. 1) then
            call getvtx(' ', 'TOUT', nbval=0, nbret=n1)
            if (n1 .ne. 0) goto 9998
        end if
    end if
!
!
!
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
! --- CREATION ET AFFECTATION DU VECTEUR DE K8 DE NOM LISMAI
!     CONTENANT LES NOMS DES MAILLES FORMANT LE LIGREL A CREER
!     --------------------------------------------------------
    call reliem(modele, noma, 'NU_MAILLE', motfac, iocc, &
                2, motcle(1), typmcl(1), lismai, nbma)
!
!     -- SI LES MOTS CLES GROUP_MA ET MAILLE N'ONT PAS ETE UTILISES:
    if (nbma .eq. 0) goto 9998
!
!
!
! --- CREATION DU LIGREL
!     ---------------------------------
    call getres(k16b, k16b, oper)
    if (oper .ne. 'IMPR_RESU') then
        noojb = '12345678.LIGR000000.LIEL'
        call gnomsd(' ', noojb, 14, 19)
    else
!     -- DANS LE CAS IMPR_RESU, GNOMSD NE PEUT PAS SERVIR CAR
!        LA COMMANDE NE CREE PAS DE CONCEPT
        ASSERT(base .eq. 'V')
        noojb = '&&EXLIMA.LIGR000000.LIEL'
        call gnoms2(noojb, 14, 19)
    end if
    ligrel = noojb(1:19)
    ASSERT(ligrel(1:8) .ne. ' ')
    call jeveuo(lismai, 'L', jma)
    call exlim1(zi(jma), nbma, modele, base, ligrel)
    call jedetr(lismai)
    goto 999
!
!
9998 continue
    ligrel = ligrmo
!
999 continue
!
!
end subroutine
