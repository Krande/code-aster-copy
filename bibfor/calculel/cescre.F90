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

subroutine cescre(basez, cesz, typcez, maz, nomgdz, &
                  ncmpg, licmp, npg, nspt, ncmp, undf0_)
! person_in_charge: jacques.pellet at edf.fr
! A_UTIL
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeundf.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/verigd.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: maz, nomgdz, cesz, basez, typcez
    integer(kind=8) :: npg(*), nspt(*), ncmp(*)
    character(len=*) :: licmp(*)
    aster_logical, optional, intent(in) :: undf0_
! ------------------------------------------------------------------
! BUT : CREER UN CHAM_ELEM_S VIERGE (CESZ)
! ------------------------------------------------------------------
!     ARGUMENTS:
! BASEZ   IN       K1  : BASE DE CREATION POUR CESZ : G/V/L
! CESZ    IN/JXOUT K19 : SD CHAM_ELEM_S A CREER
! TYPCEZ  IN       K4  : TYPE DU CHAM_ELEM_S :
!                        / 'ELNO'
!                        / 'ELGA'
!                        / 'ELEM'
! MAZ     IN/JXIN  K8  : SD MAILLAGE ASSOCIEE A CESZ
! NOMGDZ  IN       K8  : NOM DE LA GRANDEUR DE CESZ
! NCMPG   IN       I   : DIMENSION DE LICMP (CI-DESSOUS)
!            SI NCMPG <= 0  :  ON NE SE SERT PAS DE LICMP
!              SI NOMGDZ /= 'VARI_*' :
!                ON PREND TOUTES LES CMPS DU CATALOGUE
!              SI NOMGDZ = 'VARI_*' :
!                ON PREND LES CMPS V1,V2,...,'V'//CHAR(-NCMPG)
!
!
! LICMP   IN       L_K8: NOMS DES CMPS VOULUES DANS CESZ
!                        SI NOMGD='VARI_*' :
!                        LES CMPS DOIVENT AVOIR LA FORME : 'V1','V2',...
!
! SI TYPCEZ = 'ELGA' (SINON NPG EST INUTILISE):
! NPG     IN       V(I) : NOMBRES DE POINTS DE GAUSS POUR LES MAILLES.
!    / NPG(1)<0 : LE TABLEAU NPG EST ALORS DE DIMENSION 1
!                 ET -NPG EST LE NOMBRE DE POINTS DE GAUSS POUR TOUTES
!                 LES MAILLES DU MAILLAGE.
!    / NPG(1)>=0 : LE TABLEAU NPG EST DE DIMENSION NB_MAILLES(MAZ)
!                 NPG(IMA) EST LE NOMBRE DE POINTS VOULUS POUR LA
!                 MAILLE IMA
!
! NSPT    IN       V(I) : NOMBRES DE SOUS-POINTS POUR LES MAILLES.
!    / NSPT(1)<0 : LE TABLEAU NSPT EST ALORS DE DIMENSION 1
!                 ET -NSPT EST LE NOMBRE DE SOUS-POINTS POUR TOUTES
!                 LES MAILLES DU MAILLAGE.
!    / NSPT(1)>=0 : LE TABLEAU NSPT EST DE DIMENSION NB_MAILLES(MAZ)
!                 NSPT(IMA) EST LE NOMBRE DE SOUS-POINTS VOULUS POUR LA
!                 MAILLE IMA
!
! NCMP    IN       V(I) : NOMBRES DE CMPS POUR LES MAILLES.
!    / NCMP(1)<0 : LE TABLEAU NCMP EST ALORS DE DIMENSION 1
!                 ET -NCMP EST LE NOMBRE DE CMPS POUR TOUTES
!                 LES MAILLES DU MAILLAGE.
!    / NCMP(1)>=0 : LE TABLEAU NCMP EST DE DIMENSION NB_MAILLES(MAZ)
!                 NCMP(IMA) EST LE NOMBRE DE CMPS VOULUES POUR LA
!                 MAILLE IMA
!
! UNDFO  IN       L : Initialisation à zéro du champ
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=1) :: base
    character(len=3) :: tsca
    character(len=4) :: typces
    character(len=8) :: ma, nomgd, nomcmp
    character(len=19) :: ces
    character(len=24) :: valk(2)
    integer(kind=8) :: gd, ncmpmx, nbma, jcmpgd, icmp, jcmp, jcesk, jcesd
    integer(kind=8) :: jcesc, k, jcesl, jcesv, ncmpg, ima, jlconx, decal
    integer(kind=8) :: nptma, nsptma, ncmpma, ncmp2, jlicmp, iret
    aster_logical :: undf0
!
!     FONCTION FORMULE:
!     NBNOMA(IMA)=NOMBRE DE NOEUDS DE LA MAILLE IMA
#define nbnoma(ima) zi(jlconx-1+ima+1) - zi(jlconx-1+ima)
!     ------------------------------------------------------------------
!
    call jemarq()
    ces = cesz
    base = basez
    nomgd = nomgdz
    ma = maz

    undf0 = ASTER_FALSE
    if (present(undf0_)) then
        undf0 = undf0_
    end if
!
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
!     -- SI CES EXISTE DEJA, ON LE DETRUIT :
    call detrsd('CHAM_ELEM_S', ces)
!
!
!------------------------------------------------------------------
!     1- QUELQUES VERIFS (+ RECUPERATION DE JLCONX):
!     ----------------------------------------------
!
    typces = typcez
    if (typces .eq. 'ELEM') then
    else if (typces .eq. 'ELGA') then
    else if (typces .eq. 'ELNO') then
        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', jlconx)
    else
        ASSERT(.false.)
    end if
!
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), gd)
    if (gd .eq. 0) then
        call utmess('F', 'CALCULEL_67', sk=nomgd)
    end if
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', jcmpgd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
!
!
!
!     -- ON CALCULE ET ON VERIFIE :  '&&CESCRE.LICMP' :
!     --------------------------------------------------
    if (ncmpg .eq. 0) then
        ASSERT(nomgd(1:5) .ne. 'VARI_')
        ncmp2 = ncmpmx
        call wkvect('&&CESCRE.LICMP', 'V V K8', ncmp2, jlicmp)
        do k = 1, ncmp2
            zk8(jlicmp-1+k) = zk8(jcmpgd-1+k)
        end do
!
    else if (ncmpg .gt. 0) then
        call verigd(nomgd, licmp, ncmpg, iret)
        ASSERT(iret .le. 0)
!
        ncmp2 = ncmpg
        call wkvect('&&CESCRE.LICMP', 'V V K8', ncmp2, jlicmp)
        do k = 1, ncmp2
            zk8(jlicmp-1+k) = licmp(k)
        end do
!
    else if (ncmpg .lt. 0) then
        ASSERT(nomgd(1:5) .eq. 'VARI_')
        ncmp2 = -ncmpg
        call wkvect('&&CESCRE.LICMP', 'V V K8', ncmp2, jlicmp)
        nomcmp(1:1) = 'V'
        do k = 1, ncmp2
            call codent(k, 'G', nomcmp(2:8))
            zk8(jlicmp-1+k) = nomcmp
        end do
    end if
!
    if (ncmp2 .eq. 1 .and. zk8(jlicmp) .eq. ' ') then
        continue
    else
        do icmp = 1, ncmp2
            if (nomgd(1:5) .ne. 'VARI_') then
                jcmp = indik8(zk8(jcmpgd), zk8(jlicmp-1+icmp), 1, ncmpmx)
            else
                if (zk8(jlicmp-1+icmp) (1:1) .ne. 'V') then
                    jcmp = 0
                else
                    jcmp = 1
                end if
            end if
            if (jcmp .eq. 0) then
                valk(1) = zk8(jlicmp-1+icmp)
                valk(2) = nomgd
                call utmess('F', 'CALCULEL_52', nk=2, valk=valk)
            end if
        end do
    end if
!
!
!------------------------------------------------------------------
!     2- CREATION DE CES.CESC :
!     -------------------------------------------
    call wkvect(ces//'.CESC', base//' V K8', ncmp2, jcesc)
    do k = 1, ncmp2
        zk8(jcesc-1+k) = zk8(jlicmp-1+k)
    end do
!
!------------------------------------------------------------------
!     3- CREATION DE CES.CESK:
!     ------------------------
    call wkvect(ces//'.CESK', base//' V K8', 3, jcesk)
    zk8(jcesk-1+1) = ma
    zk8(jcesk-1+2) = nomgd
    zk8(jcesk-1+3) = typces
!
!------------------------------------------------------------------
!     4- CREATION DE CES.CESD:
!     ------------------------
    call wkvect(ces//'.CESD', base//' V I', 5+4*nbma, jcesd)
    zi(jcesd-1+1) = nbma
    zi(jcesd-1+2) = ncmp2
    if (ncmp2 .eq. 1 .and. zk8(jlicmp) .eq. ' ') zi(jcesd-1+2) = 0
    zi(jcesd-1+3) = 0
    zi(jcesd-1+4) = 0
    zi(jcesd-1+5) = 0
    decal = 0
    do ima = 1, nbma
!
!       -- CALCUL DE NPT(IMA):
        if (typces .eq. 'ELEM') then
            nptma = 1
        else if (typces .eq. 'ELNO') then
            nptma = nbnoma(ima)
        else if (typces .eq. 'ELGA') then
            if (npg(1) .lt. 0) then
                nptma = -npg(1)
            else
                nptma = npg(ima)
            end if
        end if
!
!       -- CALCUL DE NSPT(IMA):
        if (nspt(1) .lt. 0) then
            nsptma = -nspt(1)
        else
            nsptma = nspt(ima)
        end if
!
!       -- CALCUL DE NCMP(IMA):
        if (ncmp(1) .lt. 0) then
            ncmpma = -ncmp(1)
        else
            ncmpma = ncmp(ima)
        end if
!
        zi(jcesd-1+5+4*(ima-1)+1) = nptma
        zi(jcesd-1+5+4*(ima-1)+2) = nsptma
        zi(jcesd-1+5+4*(ima-1)+3) = ncmpma
        zi(jcesd-1+5+4*(ima-1)+4) = decal
!
        decal = decal+nptma*nsptma*ncmpma
!
        zi(jcesd-1+3) = max(nptma, zi(jcesd-1+3))
        zi(jcesd-1+4) = max(nsptma, zi(jcesd-1+4))
        zi(jcesd-1+5) = max(ncmpma, zi(jcesd-1+5))
    end do
!
!     -- POUR POUVOIR CONTINUER SI DECAL=0 (CES VIDE):
    decal = max(decal, 1)
!
!------------------------------------------------------------------
!     5- CREATION DE CES.CESL:
!     ------------------------
    call wkvect(ces//'.CESL', base//' V L', decal, jcesl)
    call jeundf(ces//'.CESL')
!
!------------------------------------------------------------------
!     6- CREATION DE CES.CESV:
!     ------------------------
    call wkvect(ces//'.CESV', base//' V '//tsca, decal, jcesv)
    call jeundf(ces//'.CESV', undf0)
!
!
!------------------------------------------------------------------
!     7- MENAGE :
!     ------------------------
    call jedetr('&&CESCRE.LICMP')
!
    call jedema()
end subroutine
