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

subroutine pemaxn(resu, nomcha, lieu, nomlie, modele, &
                  chpost, nbcmp, nomcmp, nuord, inst, nbmail, numemail)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/asmpi_allgather_char8.h"
#include "asterc/asmpi_allgather_r.h"
#include "asterc/asmpi_comm.h"
#include "asterc/indik8.h"
#include "asterc/r8maem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbcmp, nuord, nbmail, numemail(*)
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
    integer(kind=8) :: i, jcesl, jcmpgd, ncmpm, nbnoma
    integer(kind=8) :: icmp, nbpara, nbno, numno, iacnex, iproc
    integer(kind=8) :: ino, nmin, nmax, npara, nbcmpm, rank, nbproc
    real(kind=8) :: vmin, vmax, inst
    complex(kind=8) :: cbid
    character(len=7) :: chnuno
    character(len=8) :: noma, k8b, nomgd, nomva, knmin, knmax
    character(len=19) :: cesout
    aster_logical :: exist, l_pmesh
! Tableaux automatiques F90
    real(kind=8) :: mima(2*nbcmp+2)
    character(len=24) :: nompar(4*nbcmp+5), nomax(2*nbcmp+3)
    integer(kind=8), pointer :: list_no(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: cesc(:) => null()
    character(len=8), pointer :: v_name(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    mpi_int :: mrank, msize, count_send, count_recv, mpicom
    real(kind=8), pointer :: v_value(:) => null()
!
    call jemarq()
    cbid = (0.d0, 0.d0)
!
    call asmpi_comm('GET', mpicom)
    call asmpi_info(rank=mrank, size=msize)
    rank = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
    l_pmesh = isParallelMesh(noma)
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
    end if
!
    nompar(1) = 'CHAMP_GD'
    nompar(2) = 'NUME_ORDRE'
    nompar(3) = 'INST'
    nompar(4) = lieu
    mima(1) = inst
    nomax(1) = nomcha
    nomax(2) = nomlie
!
    call tbexip(resu, lieu, exist, k8b)
    if (.not. exist) then
        call tbajpa(resu, 1, nompar(4), 'K24')
    end if
!
! --- CALCULS DES CHAMPS SIMPLES:
    cesout = '&&PEMAXC_CESOUT'
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
    end if
!
    npara = 4*nbcmp
    nbcmpm = cnsd(2)
    if (l_pmesh) then
        call jeveuo(noma//".NUNOLG", "L", vi=nulg)
        AS_ALLOCATE(vr=v_value, size=nbproc)
        AS_ALLOCATE(vk8=v_name, size=nbproc)
    end if
!
    do i = 1, nbcmp
        vmin = r8maem()
        vmax = -r8maem()
        icmp = indik8(cesc, nomcmp(i), 1, nbcmpm)
        ASSERT(icmp .gt. 0)
        nmax = -1
        nmin = -1
        do ino = 1, nbno
            if (list_no(ino) == 1 .and. zl(jcesl+(ino-1)*nbcmpm+icmp-1)) then
                if (vmax .lt. cnsv(1+(ino-1)*nbcmpm+icmp-1)) then
                    vmax = cnsv(1+(ino-1)*nbcmpm+icmp-1)
                    nmax = ino
                end if
                if (vmin .gt. cnsv(1+(ino-1)*nbcmpm+icmp-1)) then
                    vmin = cnsv(1+(ino-1)*nbcmpm+icmp-1)
                    nmin = ino
                end if
            end if
        end do
        if (nmax .ne. -1) then
            if (l_pmesh) then
                ! Passage en numerotation globale
                call codent(nulg(nmax), "G", chnuno)
                knmax = 'N'//chnuno
            else
                knmax = int_to_char8(nmax)
            end if
        else
            knmax = ' '
        end if
        if (nmin .ne. -1) then
            if (l_pmesh) then
                ! Passage en numerotation globale
                call codent(nulg(nmin), "G", chnuno)
                knmin = 'N'//chnuno
            else
                knmin = int_to_char8(nmin)
            end if
        else
            knmin = ' '
        end if
        nomax(2+2*(i-1)+1) = knmax
        nomax(2+2*(i-1)+2) = knmin
        if (l_pmesh) then
            count_send = to_mpi_int(1)
            count_recv = count_send
            call asmpi_allgather_r([vmax], count_send, v_value, count_recv, mpicom)
            call asmpi_allgather_char8([knmax], count_send, v_name, count_recv, mpicom)
            vmax = -r8maem()
            do iproc = 1, nbproc
                if (v_value(iproc) .gt. vmax) then
                    vmax = v_value(iproc)
                    knmax = v_name(iproc)
                end if
            end do
            call asmpi_allgather_r([vmin], count_send, v_value, count_recv, mpicom)
            call asmpi_allgather_char8([knmin], count_send, v_name, count_recv, mpicom)
            vmin = r8maem()
            do iproc = 1, nbproc
                if (v_value(iproc) .lt. vmin) then
                    vmin = v_value(iproc)
                    knmin = v_name(iproc)
                end if
            end do
        end if
        mima(1+2*(i-1)+1) = vmax
        mima(1+2*(i-1)+2) = vmin
!
        nompar(4+4*(i-1)+1) = 'MAX_'//nomcmp(i)
        nompar(4+4*(i-1)+2) = 'NO_MAX_'//nomcmp(i)
        nompar(4+4*(i-1)+3) = 'MIN_'//nomcmp(i)
        nompar(4+4*(i-1)+4) = 'NO_MIN_'//nomcmp(i)
!
! ---    ON AJOUTE LES PARAMETRES MANQUANTS DANS LA TABLE:
        call tbexip(resu, nompar(4+4*(i-1)+1), exist, k8b)
        if (.not. exist) then
            call tbajpa(resu, 1, nompar(4+4*(i-1)+1), 'R')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+2), 'K16')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+3), 'R')
            call tbajpa(resu, 1, nompar(4+4*(i-1)+4), 'K16')
        end if
!
!
    end do
    if (l_pmesh) then
        AS_DEALLOCATE(vr=v_value)
        AS_DEALLOCATE(vk8=v_name)
    end if
!
! --- ON REMPLIT LA TABLE
    nbpara = 4+npara
    call tbajli(resu, nbpara, nompar, [nuord], mima, &
                [cbid], nomax, 0)
!
    call jedetr('&&PEMAXC_CESOUT')
    call jedetr('&&PEMAXC.LIST_CMP')
    call jedetr('&&PEMAXC_IND.NOEUD')
!
    call jedema()
!
end subroutine
