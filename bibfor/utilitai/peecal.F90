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

subroutine peecal(tych, resu, nomcha, lieu, nomlie, list_ma, nbma, &
                  modele, lFromResult, chpost, nbcmp, nomcmp, &
                  nomcp2, nuord, inst, iocc, ligrel, cespoi)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/chpond.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/teattr.h"
#include "asterfort/jexnum.h"
#include "asterfort/jenuno.h"
#include "asterfort/panbno.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/asmpi_comm_vect.h"
!
    integer(kind=8) :: nbcmp, nuord, iocc, nbma, list_ma(*)
    character(len=8) :: nomcmp(nbcmp), nomcp2(nbcmp), modele, lieu
    character(len=19) :: chpost, resu, cespoi, ligrel
    character(len=24) :: nomcha
    character(len=*) :: nomlie
    character(len=4) :: tych
    aster_logical, intent(in) :: lFromResult
!
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "INTEGRALE"
!     ROUTINE D'APPEL : PEEINT
!
!     BUT : REALISER LES CALCULS DE MOYENNE ET LES
!           STOCKER DANS LA TABLE
!
!     IN  RESU   : NOM DE LA TABLE
!     IN  NOMCHA : NOM SYMBOLIQUE DU CHAMP DU POST-TRAITEMENT
!                  OU NOM DU CHAM_GD
!     IN  TYCH   : TYPE DU CHAMP 'ELGA/ELEM/ELNO)
!     IN  LIEU   : LIEU DU POST-TRAITEMENT
!         (LIEU='TOUT'/'GROUP_MA'/'MAILLE')
!     IN  NOMLIE : NOM DU LIEU
!     IN  MODELE : NOM DU MODELE
!     IN  CHPOST  : NOM DU CHAMP DU POST-TRAITEMENT
!     IN  NBCMP   : NOMBRE DE COMPOSANTES
!     IN  NOMCMP  : NOM DES COMPOSANTES
!     IN  NOMCP2  : NOM DES COMPOSANTES A AFFICHER
!                   (DIFFERENT DE NOMCMP DANS CERTAINS CAS)
!     IN  NUORD   : NUMERO D'ORDRE
!     IN  INST    : INSTANT
!     IN  IOCC    : NUMERO DE L'OCCURENCE DE INTEGRALE
!     ------------------------------------------------------------------
!
    integer(kind=8) :: iret, i, jcesl, jcesd, jpoil, jpoid
    integer(kind=8) :: nucmp, jcmpgd, ncmpm, iad, jintr, jintk
    integer(kind=8) :: ipt, nbsp, nbpt, icmp, ima, nbpara, nume_ma
    integer(kind=8) :: ico, ind1, ind2, ifm, niv, ier, type_cell, nbnott(3)
    real(kind=8) :: vol, val, inst, volpt, rang, nbproc
    complex(kind=8) :: cbid
    character(len=8) :: noma, k8b, nomgd, nomva, type_inte
    character(len=4) :: dejain
    character(len=19) :: cesout
    character(len=24) :: valk(3), nomte
    aster_logical :: exist, l_red, l_pmesh
    real(kind=8), pointer :: pdsm(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: poiv(:) => null()
    integer(kind=8), pointer :: repe(:) => null()
    integer(kind=8), pointer :: v_model_elem(:) => null()
    integer(kind=8), pointer :: v_type_cell(:) => null()
    integer(kind=8), pointer :: v_maex(:) => null()
    mpi_int :: mrank, mnbproc
! -------------------------------------------------------------------------
    call jemarq()
    cbid = (0.d0, 0.d0)
    call infniv(ifm, niv)
!
    k8b = '        '
!
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    l_pmesh = isParallelMesh(noma)
    if (l_pmesh) then
        call jeveuo(noma//'.MAEX', 'L', vi=v_maex)
        call asmpi_info(rank=mrank, size=mnbproc)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(mnbproc)
    end if

    call jeveuo(ligrel//'.REPE', 'L', vi=repe)
!
! --- TABLEAUX DE TRAVAIL:
!     - TABLEAU DES PARAMETRES INTE_XXXX : ZK16(JINTK)
!     - TABLEAU DES VALEURS DES MOYENNES : ZR(JINTR)
    if (lFromResult) then
        call wkvect('&&PEECAL.INTE_R', 'V V R', 2*nbcmp+2, jintr)
        call wkvect('&&PEECAL.INTE_K', 'V V K16', 2*nbcmp+5, jintk)
        zk16(jintk) = 'NOM_CHAM'
        zk16(jintk+1) = 'NUME_ORDRE'
        zk16(jintk+2) = 'INST'
        zk16(jintk+3) = 'VOL'
        zk16(jintk+4) = lieu
        zr(jintr) = inst
        valk(1) = nomcha
        valk(2) = nomlie
        ind1 = 5
        ind2 = 1
    else
        call wkvect('&&PEECAL.INTE_R', 'V V R', 2*nbcmp, jintr)
        call wkvect('&&PEECAL.INTE_K', 'V V K16', 2*nbcmp+3, jintk)
        zk16(jintk) = 'CHAM_GD'
        zk16(jintk+1) = 'VOL'
        zk16(jintk+2) = lieu
        valk(1) = nomcha
        valk(2) = nomlie
        ind1 = 3
        ind2 = 0
    end if

!
! --- POUR LES CHAM_ELEM / ELEM : MOT CLE DEJA_INTEGRE:
    if (tych .eq. 'ELEM') then
        call getvtx('INTEGRALE', 'DEJA_INTEGRE', iocc=iocc, scal=dejain, nbret=iret)
        if (iret .eq. 0) then
            call utmess('F', 'UTILITAI7_13', sk=valk(1))
        end if
    end if
!
!
! --- CALCULS DES CHAMPS SIMPLES:
!      CESOUT: CHAMP ELXX CORRESPONDANT AU CHAMP CHPOST (SIMPLE) PONDERE
!              PAR LE POIDS*JACOBIEN.
!      CESPOI: CHAMP ELXX CORRESPONDANT AU POIDS*JACOBIEN
    cesout = '&&PEECAL.CESOUT'
    call chpond(tych, dejain, chpost, cesout, cespoi, ligrel, k8b)
    call jeveuo(cesout//'.CESV', 'L', vr=cesv)
    call jeveuo(cesout//'.CESL', 'L', jcesl)
    call jeveuo(cesout//'.CESD', 'L', jcesd)
    call jeveuo(cesout//'.CESK', 'L', vk8=cesk)
    call jeveuo(cespoi//'.CESV', 'L', vr=poiv)
    call jeveuo(cespoi//'.CESL', 'L', jpoil)
    call jeveuo(cespoi//'.CESD', 'L', jpoid)
    if (tych .ne. 'ELGA') call jeveuo(cespoi//'.PDSM', 'L', vr=pdsm)
!
!
! --- RECUPERATION DE LA LISTE DES CMPS DU CATALOGUE :
!     (POUR LA GRANDEUR VARI_* , IL FAUT CONSTITUER :(V1,V2,...,VN))
    nomgd = cesk(2)
    call jelira(cesout//'.CESC', 'LONMAX', ncmpm)
    if (nomgd(1:5) .ne. 'VARI_') then
        call jeveuo(cesout//'.CESC', 'L', jcmpgd)
    else
        call wkvect('&&PEECAL.LIST_CMP', 'V V K8', ncmpm, jcmpgd)
        do i = 1, ncmpm
            nomva = 'V'
            call codent(i, 'G', nomva(2:8))
            zk8(jcmpgd-1+i) = nomva
        end do
    end if
!
!     - INFOS
    if (niv .gt. 1) then
        write (6, *) '<PEECAL> NOMBRE DE MAILLES A TRAITER : ', nbma
        write (6, *) '<PEECAL> NOMBRE DE COMPOSANTES : ', ncmpm
    end if
!
!
    call jeveuo(modele//'.MAILLE', 'L', vi=v_model_elem)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=v_type_cell)
! --- CALCUL DE L'INTEGRALE ET DE LA MOYENNE(=INTEGRALE/VOLUME):
    do icmp = 1, nbcmp
        nucmp = indik8(zk8(jcmpgd), nomcmp(icmp), 1, ncmpm)
        val = 0.d0
        vol = 0.d0
        ico = 0
        do ima = 1, nbma
            nume_ma = list_ma(ima)
            if (l_pmesh) then
                if (v_maex(nume_ma) .ne. rang) cycle
            end if
            if (repe(2*(nume_ma-1)+1) .eq. 0) cycle
            nbpt = zi(jcesd-1+5+4*(nume_ma-1)+1)
            nbsp = zi(jcesd-1+5+4*(nume_ma-1)+2)
            l_red = ASTER_FALSE
            if (v_model_elem(nume_ma) .ne. 0) then
                call jenuno(jexnum('&CATA.TE.NOMTE', v_model_elem(nume_ma)), nomte)
                call teattr('C', 'INTTHM', type_inte, ier, typel=nomte)
                if (ier .eq. 0) then
                    l_red = type_inte .eq. 'RED'
                    if (l_red) then
                        type_cell = v_type_cell(nume_ma)
                        call panbno(type_cell, nbnott)
                        ! nbpt = nbpt - nbnott(1)
                    end if
                end if
            end if
            if (nbsp .gt. 1) then
                call utmess('F', 'UTILITAI8_60')
            end if
            do ipt = 1, nbpt
                call cesexi('C', jcesd, jcesl, nume_ma, ipt, &
                            1, nucmp, iad)
                ASSERT(iad .ge. 0)
                if (iad .eq. 0) cycle
!
                val = val+cesv(iad)
!
                if (tych .eq. 'ELGA') then
                    call cesexi('C', jpoid, jpoil, nume_ma, ipt, &
                                1, 1, iad)
                    ASSERT(iad .gt. 0)
                    volpt = poiv(iad)
                else if (tych .eq. 'ELEM') then
                    ASSERT(nbpt .eq. 1)
                    volpt = pdsm(nume_ma)
                else if (tych .eq. 'ELNO') then
                    ASSERT(nbpt .ge. 1)
                    volpt = pdsm(nume_ma)/nbpt
                end if

                if (.NOT. l_red) then
                    vol = vol+volpt
                else if ((l_red) .AND. (ipt < nbpt-nbnott(1)+1)) then
                    vol = vol+volpt
                end if

                ico = ico+1
            end do

        end do
        if (ico .eq. 0 .and. nbma .ne. 0) then
            valk(3) = nomcmp(icmp)
            call utmess('F', 'UTILITAI7_12', nk=3, valk=valk)
        end if
!
! --- Sum Value in HPC
        if (l_pmesh) then
            call asmpi_comm_vect("MPI_SUM", 'R', scr=vol)
            call asmpi_comm_vect("MPI_SUM", 'R', scr=val)
        end if
!
        if (icmp .eq. 1) zr(jintr+icmp+ind2-1) = vol
        zr(jintr+icmp+ind2) = val
        zk16(jintk+ind1+icmp-1) = 'INTE_'//nomcp2(icmp)
        zr(jintr+nbcmp+icmp+ind2) = val/vol
        zk16(jintk+ind1+nbcmp+icmp-1) = 'MOYE_'//nomcp2(icmp)
    end do
!
!
! --- ON AJOUTE LES PARAMETRES MANQUANTS DANS LA TABLE:
    call tbexip(resu, lieu, exist, k8b)
    if (.not. exist) then
        call tbajpa(resu, 1, zk16(jintk+ind1-1), 'K16')
    end if
    do icmp = 1, nbcmp*2
        call tbexip(resu, zk16(jintk+ind1+icmp-1), exist, k8b)
        if (.not. exist) then
            call tbajpa(resu, 1, zk16(jintk+ind1+icmp-1), 'R')
        end if
    end do
!
! --- ON REMPLIT LA TABLE
    nbpara = ind1+nbcmp*2
    call tbajli(resu, nbpara, zk16(jintk), [nuord], zr(jintr), &
                [cbid], valk, 0)
    call detrsd('CHAM_ELEM_S', '&&PEECAL.CESOUT')
    call jedetr('&&PEECAL.INTE_R')
    call jedetr('&&PEECAL.INTE_K')
    call jedetr('&&PEECAL.IND.MAILLE')
    call jedetr('&&PEECAL.LIST_CMP')
    call jedetr('&&PEECAL.MAILLES_FILTRE')
    call jedetr('&&PEECAL.MES_MAILLES')
!
    call jedema()
!
end subroutine
