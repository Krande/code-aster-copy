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

subroutine caraff(noma, gran, base, cartz)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8vide.h"
#include "asterfort/alcart.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
!
!
    character(len=1) :: base
    character(len=8) :: noma, gran
    character(len=*) :: cartz
! BUT :
!  - TRAITER L'OPTION 'AFFE' DE LA COMMANDE CREA_CHAMP
!    POUR LES CARTES ET LES CHAM_ELEM (SAUF POUR VARI_R)
!  - CREER LA CARTE  (CARTZ)
!-----------------------------------------------------------------------
    integer(kind=8) :: gd, nocc, ncmpmx, nbtou, n1, vali(2)
    integer(kind=8) :: iad, jncmp, jvalv, jmail, nbcmp, k, iocc, nbmail, nbvar
    real(kind=8) :: rvid
    character(len=8) :: k8b, tsca, typmcl(2)
    character(len=16) :: motclf, motcls(2)
    character(len=19) :: carte
    character(len=24) :: mesmai
!     ------------------------------------------------------------------
    call jemarq()
!
    if (noma .eq. ' ') then
        call utmess('F', 'UTILITAI_10')
    end if
!
    if (gran .eq. 'VARI_R') then
        call utmess('F', 'UTILITAI_11')
    end if
!
    call dismoi('TYPE_SCA', gran, 'GRANDEUR', repk=tsca)
!
    motclf = 'AFFE'
    call getfac(motclf, nocc)
!
    mesmai = '&&CARAFF.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
!     1- ALLOCATION DE LA CARTE
!     --------------------------------------------
    carte = cartz
    call alcart(base, carte, noma, gran)
    call jeveuo(carte//'.NCMP', 'E', jncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
    call jenonu(jexnom('&CATA.GD.NOMGD', gran), gd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad)
!
    if (gran .eq. 'VAR2_R') then
        rvid = r8vide()
        nbcmp = ncmpmx
        do k = 1, ncmpmx
            zk8(jncmp-1+k) = zk8(iad-1+k)
            zr(jvalv-1+k) = rvid
        end do
        call nocart(carte, 1, nbcmp)
    end if
!
!     2- BOUCLE SUR LES OCCURENCES DU MOT CLE AFFE
!     --------------------------------------------
    do iocc = 1, nocc
!
        call getvtx(motclf, 'NOEUD', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call utmess('F', 'UTILITAI_12')
        end if
!
        call getvtx(motclf, 'GROUP_NO', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call utmess('F', 'UTILITAI_13')
        end if
!
        call getvtx(motclf, 'NOM_CMP', iocc=iocc, nbval=0, nbret=nbcmp)
!
        if (tsca .eq. 'R') then
            call getvr8(motclf, 'VALE', iocc=iocc, nbval=0, nbret=nbvar)
        else if (tsca .eq. 'I') then
            call getvis(motclf, 'VALE_I', iocc=iocc, nbval=0, nbret=nbvar)
        else if (tsca .eq. 'C') then
            call getvc8(motclf, 'VALE_C', iocc=iocc, nbval=0, nbret=nbvar)
        else if (tsca .eq. 'K8') then
            call getvid(motclf, 'VALE_F', iocc=iocc, nbval=0, nbret=nbvar)
        else
            call utmess('F', 'UTILITAI_14', sk=tsca)
        end if
!
!       TEST SUR LES DONNEES INTRODUITES
        if (nbvar .ne. nbcmp) then
            call utmess('F', 'UTILITAI_15')
        else if (-nbvar .gt. ncmpmx) then
            vali(1) = -nbvar
            vali(2) = ncmpmx
            call utmess('F', 'UTILITAI_8', ni=2, vali=vali)
        else
            nbcmp = -nbcmp
            nbvar = -nbvar
            call getvtx(motclf, 'NOM_CMP', iocc=iocc, nbval=nbcmp, vect=zk8(jncmp))
            if (tsca .eq. 'R') then
                call getvr8(motclf, 'VALE', iocc=iocc, nbval=nbvar, vect=zr(jvalv))
            else if (tsca .eq. 'I') then
                call getvis(motclf, 'VALE_I', iocc=iocc, nbval=nbvar, vect=zi(jvalv))
            else if (tsca .eq. 'C') then
                call getvc8(motclf, 'VALE_C', iocc=iocc, nbval=nbvar, vect=zc(jvalv))
            else if (tsca .eq. 'K8') then
                call getvid(motclf, 'VALE_F', iocc=iocc, nbval=nbvar, vect=zk8(jvalv))
            end if
        end if
!
        call getvtx(motclf, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(carte, 1, nbcmp)
!
        else
            call reliem(' ', noma, 'NU_MAILLE', motclf, iocc, &
                        2, motcls, typmcl, mesmai, nbmail)
            if (nbmail .eq. 0) goto 30
            call jeveuo(mesmai, 'L', jmail)
            call nocart(carte, 3, nbcmp, mode='NUM', nma=nbmail, &
                        limanu=zi(jmail))
            call jedetr(mesmai)
        end if
30      continue
    end do
!
    call tecart(carte)
    call jedetr(carte//'.NCMP')
    call jedetr(carte//'.VALV')
!
    call jedema()
end subroutine
