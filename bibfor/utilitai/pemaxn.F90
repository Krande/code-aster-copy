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

subroutine pemaxn(resu, nomcha, lieu, nomlie, modele,&
                  chpost, nbcmp, nomcmp, nuord, inst, nbmail, numemail)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
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
!     BUT : EXTRAIRE LE MIN ET LE MAX D'UNE CMP D'UN CHAMNO
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
!     ------------------------------------------------------------------
!
    integer :: i, jcesl, jcmpgd, ncmpm, nbnoma
    integer :: icmp, nbpara, nbno, numno, iacnex
    integer :: ino, nmin, nmax, npara, nbcmpm
    real(kind=8) :: vmin, vmax, inst
    complex(kind=8) :: cbid
    character(len=8) :: noma, k8b, nomgd, nomva, knmin, knmax
    character(len=19) :: cesout
    character(len=24) :: nomnoe
    aster_logical :: exist
! Tableaux automatiques F90
    real(kind=8) :: mima(2*nbcmp+2)
    character(len=24) :: nompar(4*nbcmp+5), nomax(2*nbcmp+3)
    integer, pointer :: list_no(:) => null()
    integer, pointer :: cnsd(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: cesc(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
!
    call jemarq()
    cbid=(0.d0,0.d0)
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
!
! --- CREATION D'UN TABLEAU D'INDICES POUR REPERER
!     LES MAILLES DU POST TRAITEMENT
    call wkvect('&&PEMAXC_IND.NOEUD', 'V V I', nbno, vi=list_no)
    if (lieu == 'GROUP_MA') then
        list_no(:) = 0
        do i = 1, nbmail
            call jeveuo(jexnum(noma//'.CONNEX', numemail(i)), 'L', iacnex)
            call jelira(jexnum(noma//'.CONNEX', numemail(i)), 'LONMAX', nbnoma)
            do ino = 1, nbnoma
                numno = zi(iacnex-1+ino)
                list_no(numno) = 1
            end do
        end do
    elseif (lieu == 'TOUT') then
        list_no(:) = 1
    else
        ASSERT(ASTER_FALSE)
    endif
!
    nomnoe = noma//'.NOMNOE         '
    nompar(1)='CHAMP_GD'
    nompar(2)='NUME_ORDRE'
    nompar(3)='INST'
    nompar(4)=lieu
    mima(1)=inst
    nomax(1)=nomcha
    nomax(2)=nomlie
!
    call tbexip(resu, lieu, exist, k8b)
    if (.not.exist) then
        call tbajpa(resu, 1, nompar(4), 'K24')
    endif
!
! --- CALCULS DES CHAMPS SIMPLES:
    cesout='&&PEMAXC_CESOUT'
    call cnocns(chpost, 'V', cesout)
    call jeveuo(cesout//'.CNSV', 'L', vr=cnsv)
    call jeveuo(cesout//'.CNSL', 'L', jcesl)
    call jeveuo(cesout//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cesout//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cesout//'.CNSC', 'L', vk8=cesc)
!
! --- RECUPERATION DE LA LISTE DES CMPS DU CATALOGUE :
!     (POUR LA GRANDEUR VARI_* , IL FAUT CONSTITUER :(V1,V2,...,VN))
    nomgd = cnsk(2)
    call jelira(cesout//'.CNSC', 'LONMAX', ncmpm)
    if (nomgd(1:5) .ne. 'VARI_') then
        call jeveuo(cesout//'.CNSC', 'L', jcmpgd)
    else
        call wkvect('&&PEMAXC.LIST_CMP', 'V V K8', ncmpm, jcmpgd)
        do i = 1, ncmpm
            nomva = 'V'
            call codent(i, 'G', nomva(2:8))
            zk8(jcmpgd-1+i) = nomva
        end do
    endif
!
    npara=4*nbcmp
    nbcmpm=cnsd(2)
!
    do i = 1, nbcmp
        vmin=r8maem()
        vmax=-r8maem()
        icmp=indik8(cesc,nomcmp(i),1,nbcmpm)
        ASSERT(icmp.gt.0)
        do ino = 1, nbno
            if (list_no(ino) == 1 .and. zl(jcesl+(ino-1)*nbcmpm+icmp-1)) then
                if (vmax .lt. cnsv(1+(ino-1)*nbcmpm+icmp-1)) then
                    vmax=cnsv(1+(ino-1)*nbcmpm+icmp-1)
                    nmax=ino
                endif
                if (vmin .gt. cnsv(1+(ino-1)*nbcmpm+icmp-1)) then
                    vmin=cnsv(1+(ino-1)*nbcmpm+icmp-1)
                    nmin=ino
                endif
            endif
        end do
        mima(1+2*(i-1)+1)=vmax
        mima(1+2*(i-1)+2)=vmin
        call jenuno(jexnum(nomnoe, nmin), knmin)
        call jenuno(jexnum(nomnoe, nmax), knmax)
        nomax(2+2*(i-1)+1)=knmax
        nomax(2+2*(i-1)+2)=knmin
!
        nompar(4+4*(i-1)+1)='MAX_'//nomcmp(i)
        nompar(4+4*(i-1)+2)='NO_MAX_'//nomcmp(i)
        nompar(4+4*(i-1)+3)='MIN_'//nomcmp(i)
        nompar(4+4*(i-1)+4)='NO_MIN_'//nomcmp(i)
!
! ---    ON AJOUTE LES PARAMETRES MANQUANTS DANS LA TABLE:
        call tbexip(resu, nompar(4+4*(i-1)+1), exist, k8b)
        if (.not.exist) then
            call tbajpa(resu, 1, nompar(4+4*(i-1)+1), 'R')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+2), 'K16')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+3), 'R')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+4), 'K16')
        endif
!
!
    end do
!
! --- ON REMPLIT LA TABLE
    nbpara=4+npara
    call tbajli(resu, nbpara, nompar, [nuord], mima,&
                [cbid], nomax, 0)
!
    call jedetr('&&PEMAXC_CESOUT')
    call jedetr('&&PEMAXC.LIST_CMP')
    call jedetr('&&PEMAXC_IND.NOEUD')
!
    call jedema()
!
end subroutine
