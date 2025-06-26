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

subroutine pjelga(nomo2, cham1, ligre1, prolong, corres, &
                  leres1, ligre2, iret)
!
! --------------------------------------------------------------------------------------------------
!
!       COMMANDE:  PROJ_CHAMP
!
!               ROUTINE "CHAPEAU" CONCERNANT LA PROJECTION DE CHAM_ELEM (ELGA)
!
! --------------------------------------------------------------------------------------------------
!
!  ELLE EST CONSTITUEE DE TROIS TEMPS
!    APPEL A ECLPGC (ECLATEMENT DU CHAMP)
!    APPEL A PJXXCH (USUEL POUR TOUS LES CHAMPS, CF. OP0166 OU PJXXPR)
!    APPEL A PJCORR (RETOUR AUX POINTS DE GAUSS)
!
! --------------------------------------------------------------------------------------------------
!
    use proj_champ_module
    implicit none
!
#include "asterfort/celfpg.h"
#include "asterfort/cescel.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/eclpgc.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pjcorr.h"
#include "asterfort/pjxxch.h"
#include "asterfort/titre.h"
#include "asterfort/xcesrd.h"
#include "asterfort/xnpgxx.h"
!
    character(len=8)    :: nomo2
    type(prolongation)  :: prolong
    character(len=16)   :: corres
    character(len=19)   :: cham1, ligre1
    character(len=19)   :: leres1, ligre2
    integer(kind=8)             :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)             :: nncp, jcnsv
    character(len=3)    :: exixfm
    character(len=4)    :: tycha2
    character(len=8)    :: ma1p, prol0, nompar
    character(len=16)   :: option
    character(len=19)   :: cham1e, chauxs, chbid, prfchn, cns1, ch2s, chsnpg
    character(len=24)   :: nomfpg
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call titre()
!
    chbid = cham1
    nomfpg = '&&OP0166.NOMFPG'
    call celfpg(cham1, nomfpg, iret)
!
!   NOM DE L'OPTION ET NOM DU PARAMETRE POUR LE CHAM_ELEM VIVANT SUR LE MAILLAGE 1
    call dismoi('NOM_OPTION', cham1, 'CHAM_ELEM', repk=option)
    call dismoi('NOM_PARAM', cham1, 'CHAM_ELEM', repk=nompar)
!
!   ECLATEMENT DU CHAMP VIVANT SUR LE MAILLAGE 1
!   EN UN CHAMP VIVANT SUR LE MAILLAGE 1 PRIME (MAILLAGE ECLATE)
    ma1p = '&&PJELC1'
    cham1e = '&&PJELGA'//'.CHAM1E'
!
    prfchn = '&&PJELGA.PRFCHN'
    call eclpgc(cham1, cham1e, ligre1, ma1p, prfchn, nomfpg)
!
    chauxs = '&&PJELGA'//'.CHAS'
    tycha2 = ' '
    call pjxxch(corres, cham1e, chauxs, tycha2, ' ', prolong, ' ', 'G', iret)
!
!   ON TRANSFORME LE CHAM_NO PROJETE EN UN CHAM_NO_S
    cns1 = '&&PJELGA'//'.CH1S'
    call cnocns(chauxs, 'G', cns1)
    call jeveuo(cns1//'.CNSV', 'L', jcnsv)
!
!   IL FAUT MAINTENANT REVENIR AUX CHAM_ELEM (ELGA)
    ch2s = '&&OP0166'//'.CH2S'
    call pjcorr(nomo2, chbid, cns1, ch2s, ligre2, corres, option, nompar, iret)
!
!   construction du CHAM_ELEM_S contenant le nombre de points
!   de Gauss réellement utilisés par chaque élément, dans le
!   cas d'un champs ELGA, base sur la famille "XFEM"
    chsnpg = '&&CHPCHD.CHSNPG'
    call xnpgxx(nomo2, ligre2, option, nompar, chsnpg, exixfm)
!
    if (exixfm .eq. 'OUI') then
!        si le champ ELGA s'appuie sur la famille "XFEM", on
!        annule toutes les composantes associées aux points de Gauss inutilisés
        call xcesrd(ch2s, chsnpg)
    end if
!
    call detrsd('CHAM_ELEM_S', chsnpg)
!
!
!   'prolong' est donné alors 'prol0' ne sert pas
    prol0 = '????'
    call cescel(ch2s, ligre2, option, nompar, prol0, nncp, 'G', leres1, 'A', iret, prolong)
!
    call detrsd('MAILLAGE', ma1p)
    call detrsd('CHAM_NO_S', cns1)
!
    call dismoi('NUME_EQUA', cham1e, 'CHAM_NO', repk=prfchn)
    call detrsd('NUME_EQUA', prfchn)
    call detrsd('CHAM_NO', cham1e)
!
    call dismoi('NUME_EQUA', chauxs, 'CHAM_NO', repk=prfchn)
    call detrsd('NUME_EQUA', prfchn)
    call detrsd('CHAM_NO', chauxs)
!
    call jedema()
end subroutine
