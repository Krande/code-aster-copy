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

subroutine cescns(cesz, celfpz, base, cnsz, comp, &
                  cret)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterc/r8vide.h"
!
    character(len=*) :: cnsz, cesz, base, celfpz
    character(len=1) :: comp
    integer(kind=8) :: cret
! ------------------------------------------------------------------
! BUT: TRANSFORMER UN CHAM_ELEM_S EN CHAM_NO_S
! ------------------------------------------------------------------
!     ARGUMENTS:
! CESZ   IN/JXIN  K19 : SD CHAM_ELEM_S A TRANSFORMER
!
! CELFPZ IN/JXIN  K24 :
!    NOM DE L'OBJET DECRIVANT LES FAMILLES DE P.G. DE CESZ (OU ' ')
!    CET OBJET N'EST UTILISE QUE SI CESZ EST 'ELGA'
!    CET OBJET EST OBTENU PAR LA ROUTINE CELFPG.F
!
! CNSZ   IN/JXOUT K19 : SD CHAM_NO_S RESULTAT
! BASE   IN       K1  : BASE DE CREATION POUR CNSZ : G/V/L
! COMP   IN           : COMPORTEMENT EN PRESENCE DE SOUS-POINTS
!    'F' : EMISSION D'UNE ERREUR <F>
!    'A' : EMISSION D'UNE ALARME POUR PREVENIR L'UTILISATEUR
!    ' ' : SILENCE => CODE RETOUR
! CRET   OUT      I   : CODE RETOUR
!    0   : SI TOUT EST OK
!    100 : EN PRESENCE D'UN CHAMP ELEM A SOUS-POINTS
!-----------------------------------------------------------------------
!
!  PRINCIPES RETENUS POUR LA CONVERSION :
!
!  1) ON NE TRAITE QUE LES CHAM_ELEM_S REELS OU COMPLEXES
!
!  2) ON SE RAMENE TOUJOURS A UN CHAMP ELNO
!     PUIS ON FAIT LA MOYENNE ARITHMETIQUE DES MAILLES
!     QUI CONCOURRENT EN 1 MEME NOEUD.
!
!  3) S'IL Y A DES SOUS POINTS, ON S'ARRETE EN ERREUR <F>
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ima, ncmp, icmp, jcnsl, jcnsv, ispt
    integer(kind=8) :: jcesd, jcesv, jcesl, nbma, iret, nbno
    integer(kind=8) :: ino, nuno, nbpt, iad1, ilcnx1
    integer(kind=8) :: jcesk, nbnot, jnbno, ieq, nbsp
    character(len=3) :: tsca
    character(len=8) :: ma, nomgd
    character(len=19) :: ces, cns, ces1
    character(len=8), pointer :: cesc(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    ces = cesz
    cns = cnsz
!
    call jeveuo(ces//'.CESK', 'L', jcesk)
!
!     1- ON TRANSFORME CES EN CHAM_ELEM_S/ELNO:
!     --------------------------------------------
    ces1 = '&&CESCNS.CES1'
    call cesces(ces, 'ELNO', ' ', ' ', celfpz, &
                'V', ces1)
!
!
!     2. RECUPERATION DE :
!        MA     : NOM DU MAILLAGE
!        NOMGD  : NOM DE LA GRANDEUR
!        NCMP   : NOMBRE DE CMPS DE CES1
!        TSCA   : TYPE SCALAIRE DE LA GRANDEUR : R/C
!        NBMA   : NOMBRE DE MAILLES DU MAILLAGE
!        ILCNX1,IACNX1   : ADRESSES DE LA CONNECTIVITE DU MAILLAGE
!     --------------------------------------------------------------
    call exisd('CHAM_ELEM_S', ces1, iret)
    ASSERT(iret .gt. 0)
    call jeveuo(ces1//'.CESK', 'L', jcesk)
    call jeveuo(ces1//'.CESC', 'L', vk8=cesc)
    call jeveuo(ces1//'.CESD', 'L', jcesd)
    call jeveuo(ces1//'.CESV', 'L', jcesv)
    call jeveuo(ces1//'.CESL', 'L', jcesl)
    ma = zk8(jcesk-1+1)
    nomgd = zk8(jcesk-1+2)
!     TEST SI CHAMP ELNO
    ASSERT(zk8(jcesk-1+3) .eq. 'ELNO')
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbnot)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call jeveuo(ma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jelira(ces1//'.CESC', 'LONMAX', ncmp)
!
    cret = 0
!     EMISSION D'UN MESSAGE POUR SIGNIFIER QU'ON VA FILTRER
!     LES MAILLES CONTENANT DES SOUS-POINTS
    if (zi(jcesd-1+4) .gt. 1) then
        if (comp .ne. ' ') then
            call utmess(comp, 'UTILITAI_3')
        else
            cret = 100
        end if
    end if
!
!     ON ATTEND SEULEMENT DES REELS OU DES COMPLEXES
    ASSERT((tsca .eq. 'R') .or. (tsca .eq. 'C'))
!
!
!     3- ALLOCATION DE CNS :
!     -------------------------------------------
    if (nomgd .eq. 'VARI_R') nomgd = 'VAR2_R'
    call cnscre(ma, nomgd, ncmp, cesc, base, &
                cns)
!
!
!     4- REMPLISSAGE DE CNS.CNSL ET CNS.CNSV :
!     -------------------------------------------
    call jeveuo(cns//'.CNSL', 'E', jcnsl)
    call jeveuo(cns//'.CNSV', 'E', jcnsv)
!
    do icmp = 1, ncmp
        call jedetr('&&CESCNS.NBNO')
        call wkvect('&&CESCNS.NBNO', 'V V I', nbnot, jnbno)
!
        do ima = 1, nbma
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            nbsp = zi(jcesd-1+5+4*(ima-1)+2)
            nbno = zi(ilcnx1+ima)-zi(ilcnx1-1+ima)
!
            ASSERT(nbno .eq. nbpt)
            if (nbsp .eq. 1) then
                do ino = 1, nbno
                    call cesexi('C', jcesd, jcesl, ima, ino, &
                                1, icmp, iad1)
                    if (iad1 .le. 0) goto 10
!
                    nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                    ieq = (nuno-1)*ncmp+icmp
                    zl(jcnsl-1+ieq) = .true.
                    if (tsca .eq. 'R') then
                        if (zi(jnbno-1+nuno) .eq. 0) zr(jcnsv-1+ieq) = 0.d0
                        if (zr(jcesv-1+iad1) .eq. r8vide()) then
                            zr(jcnsv-1+ieq) = r8vide()
                            goto 10
                        end if
                        zr(jcnsv-1+ieq) = zr(jcnsv-1+ieq)+zr(jcesv-1+iad1)
                    else if (tsca .eq. 'C') then
                        if (zi(jnbno-1+nuno) .eq. 0) zc(jcnsv-1+ieq) = (0.d0, 0.d0)
                        zc(jcnsv-1+ieq) = zc(jcnsv-1+ieq)+zc(jcesv-1+iad1)
                    end if
                    zi(jnbno-1+nuno) = zi(jnbno-1+nuno)+1
!
10                  continue
                end do
            else
                do ino = 1, nbno
                    do ispt = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ino, &
                                    ispt, icmp, iad1)
                        if (iad1 .le. 0) goto 60
                        nuno = connex(1+zi(ilcnx1-1+ima)-2+ino)
                        ieq = (nuno-1)*ncmp+icmp
                        zl(jcnsl-1+ieq) = .false.
60                      continue
                    end do
                end do
            end if
        end do
!
        do nuno = 1, nbnot
            ieq = (nuno-1)*ncmp+icmp
            if (zl(jcnsl-1+ieq)) then
                if (tsca .eq. 'R') then
                    if (zr(jcnsv-1+ieq) .ne. r8vide()) then
                        zr(jcnsv-1+ieq) = zr(jcnsv-1+ieq)/zi(jnbno-1+nuno)
                    end if
                else if (tsca .eq. 'C') then
                    zc(jcnsv-1+ieq) = zc(jcnsv-1+ieq)/zi(jnbno-1+nuno)
                end if
            end if
        end do
!
    end do
!
!
!     7- MENAGE :
!     -----------
    call detrsd('CHAM_ELEM_S', ces1)
    call jedetr('&&CESCNS.NBNO')
!
    call jedema()
end subroutine
