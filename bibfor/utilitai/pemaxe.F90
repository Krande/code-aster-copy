! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine pemaxe(resu, nomcha, lieu, nomlie, modele,&
                  chpost, nbcmp, nomcmp, nuord, inst,&
                  nbmail, numemail)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
#include "asterfort/wkvect.h"
!
    integer :: nbcmp, nuord, nbmail, numemail(*)
    character(len=8) :: nomcmp(nbcmp), modele, lieu
    character(len=19) :: chpost, resu
    character(len=24) :: nomcha, nomlie
!
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "MINMAX"
!     ROUTINE D'APPEL : PEMIMA
!
!     BUT : EXTRAIRE LE MIN ET LE MAX D'UNE CMP D'UN CHAMELEM
!           ET LES STOCKER DANS LA TABLE
!
!     IN  RESU   : NOM DE LA TABLE
!     IN  NOMCHA : NOM SYMBOLIQUE DU CHAMP DU POST-TRAITEMENT
!     IN  LIEU   : LIEU DU POST-TRAITEMENT
!         (LIEU='TOUT'/'GROUP_MA'/'MAILLE')
!     IN  NOMLIE : NOM DU LIEU
!     IN  MODELE : NOM DU MODELE
!     IN  CHPOST  : NOM DU CHAMP DU POST-TRAITEMENT
!     IN  NBCMP   : NOMBRE DE COMPOSANTES
!     IN  NOMCMP  : NOM DES COMPOSANTES
!     IN  NUORD   : NUMERO D'ORDRE
!     IN  INST    : INSTANT
!     IN  IOCC    : NUMERO DE L'OCCURENCE
!     IN  NBAMIL  : NOMBRE DE MAILLE DE LA LISTE A TRAITER
!     IN  NUMEMAIL : NUMERO DES MAILLES A TRAITER
!     ------------------------------------------------------------------
!
    integer :: nbma, i, jcesl, jcesd
    integer :: nucmp, jcmpgd, ncmpm, iad
    integer :: ipt, nbsp, nbpt, icmp, ima, nbpara
    integer :: nmin, nmax, npara, pmax, pmin
    real(kind=8) :: vmin, vmax, inst
    complex(kind=8) :: cbid
    character(len=8) :: noma, k8b, nomgd, nomva, knmin, knmax
    character(len=19) :: cesout
    character(len=24) :: nommai
    aster_logical :: exist
! Tableaux automatiques F90
    real(kind=8) :: mima(2*nbcmp+2)
    character(len=24) :: nompar(6*nbcmp+5), mamax(2*nbcmp+3)
    integer :: ptmax(1+2*nbcmp)
    character(len=8), pointer :: cesk(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    integer, pointer :: list_ma(:) => null()
!
    call jemarq()
!
    cbid=(0.d0,0.d0)
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
!
    nommai = noma//'.NOMMAI'
!
! --- CREATION D'UN TABLEAU D'INDICES POUR REPERER
!     LES MAILLES DU POST TRAITEMENT
    call wkvect('&&PEMAXC_IND.MAILLE', 'V V I', nbma, vi=list_ma)
    if (lieu == 'GROUP_MA') then
        list_ma(:) = 0
        do i = 1, nbmail
            list_ma(numemail(i)) = 1
        end do
    elseif (lieu == 'TOUT') then
        list_ma(:) = 1
    else
        ASSERT(ASTER_FALSE)
    endif
!
    nompar(1)='CHAMP_GD'
    nompar(2)='NUME_ORDRE'
    nompar(3)='INST'
    nompar(4)=lieu
    mima(1)=inst
    mamax(1)=nomcha
    mamax(2)=nomlie
!
    call tbexip(resu, lieu, exist, k8b)
    if (.not.exist) then
        call tbajpa(resu, 1, nompar(4), 'K24')
    endif
!
! --- CALCULS DES CHAMPS SIMPLES:
    cesout='&&PEMAXC_CESOUT'
!
    call celces(chpost, 'V', cesout)
    call jeveuo(cesout//'.CESV', 'L', vr=cesv)
    call jeveuo(cesout//'.CESL', 'L', jcesl)
    call jeveuo(cesout//'.CESD', 'L', jcesd)
    call jeveuo(cesout//'.CESK', 'L', vk8=cesk)
!
!
! --- RECUPERATION DE LA LISTE DES CMPS DU CATALOGUE :
!     (POUR LA GRANDEUR VARI_* , IL FAUT CONSTITUER :(V1,V2,...,VN))
    nomgd = cesk(2)
    call jelira(cesout//'.CESC', 'LONMAX', ncmpm)
    if (nomgd(1:5) .ne. 'VARI_') then
        call jeveuo(cesout//'.CESC', 'L', jcmpgd)
    else
        call wkvect('&&PEMAXC.LIST_CMP', 'V V K8', ncmpm, jcmpgd)
        do i = 1, ncmpm
            nomva = 'V'
            call codent(i, 'G', nomva(2:8))
            zk8(jcmpgd-1+i) = nomva
        end do
    endif
!
    do icmp = 1, nbcmp
        nucmp=indik8(zk8(jcmpgd),nomcmp(icmp),1,ncmpm)
        vmin=r8maem()
        vmax=-r8maem()
!
        do ima = 1, nbma
            if (list_ma(ima) == 1) then
                nbpt=zi(jcesd-1+5+4*(ima-1)+1)
                nbsp=zi(jcesd-1+5+4*(ima-1)+2)
                ASSERT(nbsp.eq.1)
                do ipt = 1, nbpt
                    call cesexi('C', jcesd, jcesl, ima, ipt,&
                                1, nucmp, iad)
                    if (iad .gt. 0) then
                        if (vmax .lt. cesv(iad)) then
                            vmax=cesv(iad)
                            nmax=ima
                            pmax=ipt
                        endif
                        if (vmin .gt. cesv(iad)) then
                            vmin=cesv(iad)
                            nmin=ima
                            pmin=ipt
                        endif
                    endif
                end do
            endif
        end do
!
        nompar(4+6*(icmp-1)+1)='MAX_'//nomcmp(icmp)
        nompar(4+6*(icmp-1)+2)='MA_MAX_'//nomcmp(icmp)
        nompar(4+6*(icmp-1)+3)='PT_MAX_'//nomcmp(icmp)
        nompar(4+6*(icmp-1)+4)='MIN_'//nomcmp(icmp)
        nompar(4+6*(icmp-1)+5)='MA_MIN_'//nomcmp(icmp)
        nompar(4+6*(icmp-1)+6)='PT_MIN_'//nomcmp(icmp)
!
        mima(1+2*(icmp-1)+1)=vmax
        mima(1+2*(icmp-1)+2)=vmin
!
        call jenuno(jexnum(nommai, nmin), knmin)
        call jenuno(jexnum(nommai, nmax), knmax)
!
        mamax(2+2*(icmp-1)+1)=knmax
        mamax(2+2*(icmp-1)+2)=knmin
        ptmax(1+2*(icmp-1)+1)=pmax
        ptmax(1+2*(icmp-1)+2)=pmin
!
        call tbexip(resu, nompar(4+6*(icmp-1)+1), exist, k8b)
        if (.not.exist) then
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+1), 'R')
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+2), 'K16')
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+3), 'I')
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+4), 'R')
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+5), 'K16')
            call tbajpa(resu, 1, nompar(4+6*(icmp-1)+6), 'I')
        endif
    end do
!
    npara=6*nbcmp
    ptmax(1)=nuord
!
! --- ON REMPLIT LA TABLE
    nbpara=4+npara
    call tbajli(resu, nbpara, nompar, ptmax, mima,&
                [cbid], mamax, 0)
!
    call jedetr('&&PEMAXC_IND.MAILLE')
    call jedetr('&&PEMAXC_CESOUT')
    call jedetr('&&PEMAXC.LIST_CMP')
!
    call jedema()
!
end subroutine
