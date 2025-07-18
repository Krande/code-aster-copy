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

subroutine rvcpnc(mcf, iocc, nch19, gd, typegd, &
                  nbcpc, nlscpc, nomojb, repere, option, &
                  quant, codir, dir, iret)
! aslint: disable=W1501
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/i2trgi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/num2k8.h"
#include "asterfort/numek8.h"
#include "asterfort/rvopti.h"
#include "asterfort/utmess.h"
#include "asterfort/utncmp.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=*) :: mcf
    character(len=24) :: nlscpc, nomojb, quant
    character(len=19) :: nch19
    character(len=16) :: option
    character(len=8) :: repere
    character(len=4) :: typegd
    integer(kind=8) :: iocc, iret, nbcpc, gd, codir
    real(kind=8) :: dir(*)
!
!     SAISIE DES CMP NECESSAIRES AU POST-TRAITEMENT
!     ------------------------------------------------------------------
! IN  IOCC   : I : NUMERO DE L' OCCURENCE TRAITEE
! IN  GD     : I : NUMERO DE LA GRANDEUR
! IN  NCH19  : K : NOM DU CHAMP A POST-TRAITER
! IN  TYPEGD : K : VAUT 'CHNO' OU 'CHLM'
! IN  NBCPC  : I : NOMBRE DE CMP A POST-TRAITER
! IN  NLSCPC : K : NOM OJB V K8 CONTENANT LE NOM DES CMP CANDIDATES
! OUT NOMOJB : K : NOM OJB V K8 CONTENANT LE NOM DES CMP NECESSAIRES
! OUT IRET   : I : CODE RETOUR : 1 RAS, 0 ERREUR (EMISSION MESSAGE)
! OUT REPER  : K : TYPE DE REPERE UTILISE
! OUT QUANT  : K : NOM DE LA QUANTITE A TRAITER
! OUT CODIR  : I : CODE DES DIR ACT CAS D' UNE TRACE DIRECTION
!     ------------------------------------------------------------------
!
!
!
    character(len=24) :: nomaux, nomnew, nomob1
    character(len=24) :: valk(23)
    character(len=8) :: nomgd, k8b, mailla
    character(len=4) :: docu
    integer(kind=8) :: acpgd, ntc, nin, ntn1, ntn2, nep, nnc, avk8, i, nbcpgd, nso
    integer(kind=8) :: alsi, n1, n2, n3, ntd1, pt, avicp, avinew, alscpc, ptnc
    integer(kind=8) :: ancpu, nbc, iadt1, iadt2, ibid, ntd2, iexi, nn
    aster_logical :: dirx, diry, dirz
    character(len=8), pointer :: tmp(:) => null()
!
!======================================================================
!
    call jemarq()
    codir = 0
    n1 = 0
    n2 = 0
    n3 = 0
    i = 0
    acpgd = 0
    avk8 = 0
    nbcpgd = 0
    ntc = 0
    nin = 0
    ntn1 = 0
    ntn2 = 0
    ntd1 = 0
    ntd2 = 0
    nep = 0
    nnc = 0
    nso = 0
    call jeexin(nomojb, iexi)
    if (iexi .ne. 0) call jedetr(nomojb)
    nomnew = '&&RVCPNC.NUM.CMP.COURANT'
    nomaux = '&&RVCPNC.CMP.NC.TEMPORAI'
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', acpgd)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nbcpgd)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nomgd)
    call jeveuo(nlscpc, 'L', alscpc)
    if (mcf(1:6) .eq. 'ACTION') then
        call getvtx(mcf, 'TOUT_CMP', iocc=iocc, nbval=0, nbret=ntc)
        call getvtx(mcf, 'NOM_CMP', iocc=iocc, nbval=0, nbret=nnc)
        call getvtx(mcf, 'INVARIANT', iocc=iocc, nbval=0, nbret=nin)
        call getvtx(mcf, 'ELEM_PRINCIPAUX', iocc=iocc, nbval=0, nbret=nep)
        call getvtx(mcf, 'TRAC_NOR', iocc=iocc, nbval=0, nbret=ntn2)
        call getvtx(mcf, 'TRAC_DIR', iocc=iocc, nbval=0, nbret=ntd2)
        call getvtx(mcf, 'RESULTANTE', iocc=iocc, nbval=0, nbret=nso)
        call getvtx(mcf, 'REPERE', iocc=iocc, scal=repere, nbret=i)
        if (repere .eq. 'UTILISAT') then
!           --- REPERE TRAITE DANS "EXTCHE" OU "EXTCHN" ---
            repere = 'GLOBAL  '
        else if (repere .eq. 'CYLINDRI') then
!           --- REPERE TRAITE DANS "EXTCHE" OU "EXTCHN" ---
            repere = 'GLOBAL  '
        end if
    else
        nin = 1
    end if
    if (nin .ne. 0) then
        quant = 'INVARIANTS_TENSORIELS'
        repere = 'GLOBAL'
    else if (nep .ne. 0) then
        quant = 'ELEMENTS_PRINCIPAUX'
        repere = 'GLOBAL'
    else if ((ntn1 .ne. 0) .or. (ntn2 .ne. 0)) then
        call dismoi('NOM_MAILLA', nch19, 'CHAMP', repk=mailla)
        call dismoi('Z_CST', mailla, 'MAILLAGE', repk=k8b)
        if (k8b(1:3) .eq. 'NON') then
            call utmess('F', 'POSTRELE_32')
        end if
        quant = 'TRACE_NORMALE'
        repere = 'LOCAL'
    else if ((ntd1 .ne. 0) .or. (ntd2 .ne. 0)) then
        quant = 'TRACE_DIRECTIONNELLE'
        repere = 'GLOBAL'
    else if (nso .ne. 0) then
        quant = 'SOMME'
        repere = 'GLOBAL'
    else if ((ntc .ne. 0) .or. (nnc .ne. 0)) then
        quant = 'ASTER'
    else
        quant = '&&INCONNU'
    end if
    iret = 1
!
    call rvopti(mcf, iocc, nch19, nomgd, typegd, &
                option)
!
    if (ntc .ne. 0) then
        nomob1 = '&&OP0051.NOMCMP.USER'
        call utncmp(nch19, nbc, nomob1)
        if (nbc .eq. 0) then
            call utmess('F', 'POSTRELE_54', si=iocc)
        end if
        call jeveuo(nomob1, 'L', ancpu)
        call wkvect(nomojb, 'V V K8', nbc, avk8)
        do i = 1, nbc, 1
            zk8(avk8+i-1) = zk8(ancpu+i-1)
        end do
        call jedetr(nomob1)
    else if ((nin .ne. 0) .or. (nep .ne. 0)) then
        if ((option .eq. 'SIGM_ELNO') .or. (option .eq. 'SIEF_ELNO') .or. &
            (option .eq. 'EPSI_ELNO') .or. (option .eq. 'EPSG_ELNO') .or. &
            (option .eq. 'EPME_ELNO') .or. (option .eq. 'EPMG_ELNO') .or. &
            (option .eq. 'SIGM_NOEU') .or. (option .eq. 'SIEF_NOEU') .or. &
            (option .eq. 'EPSI_NOEU') .or. (option .eq. 'EPSG_NOEU') .or. &
            (option .eq. 'EPME_NOEU_DEPL') .or. (option .eq. 'EPMG_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+i-1)
            end do
        else if ((option .eq. 'EFGE_ELNO') .or. (option .eq. &
                                                 'EFGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+13+i-1)
            end do
            iret = 0
        else if ((option .eq. 'DEGE_ELNO') .or. (option .eq. &
                                                 'DEGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+6+i-1)
            end do
            iret = 0
        else
            iret = 0
            valk(1) = option
            valk(2) = ' '
            valk(3) = 'SIGM_ELNO'
            valk(4) = 'SIEF_ELNO'
            valk(5) = 'EPSI_ELNO'
            valk(6) = 'EPSG_ELNO'
            valk(7) = 'EPME_ELNO'
            valk(8) = 'EPMG_ELNO'
            valk(9) = 'EPSI_ELNO_ELGA'
            valk(10) = 'DEGE_ELNO'
            valk(11) = 'EFGE_ELNO'
            valk(12) = 'SIGM_NOEU'
            valk(13) = 'SIEF_NOEU'
            valk(14) = 'EPSI_NOEU'
            valk(15) = 'EPSG_NOEU'
            valk(16) = 'EPME_NOEU_DEPL'
            valk(17) = 'EPMG_NOEU'
            valk(18) = 'EPSI_NOEU_ELGA'
            valk(19) = 'DEGE_NOEU'
            valk(20) = 'EFGE_NOEU'
            call utmess('F', 'POSTRELE_33', nk=20, valk=valk, si=iocc)
        end if
!
    else if ((ntn1 .ne. 0) .or. (ntn2 .ne. 0)) then
!
!      /* LA NORMALE N' EST CALCULEE QUE POUR (X,Y) */
!
        if ((option .eq. 'FLUX_ELNO') .or. (option .eq. 'FLUX_NOEU_DEPL') .or. &
            (option .eq. 'FLUX_NOEU')) then
            call wkvect(nomojb, 'V V K8', 3, avk8)
            do i = 1, 3, 1
                zk8(avk8+i-1) = zk8(acpgd+i-1)
            end do
        else if ((option .eq. 'SIGM_ELNO') .or. (option .eq. 'SIEF_ELNO') .or. &
                 (option .eq. 'EPSI_ELNO') .or. (option .eq. 'EPSG_ELNO') .or. &
                 (option .eq. 'EPME_ELNO') .or. (option .eq. 'EPMG_ELNO') .or. &
                 (option .eq. 'SIGM_NOEU') .or. (option .eq. 'SIEF_NOEU') .or. &
                 (option .eq. 'EPSI_NOEU') .or. (option .eq. 'EPSG_NOEU') .or. &
                 (option .eq. 'EPME_NOEU_DEPL') .or. (option .eq. 'EPMG_NOEU')) then
            call wkvect(nomojb, 'V V K8', 5, avk8)
            zk8(avk8+1-1) = zk8(acpgd+1-1)
            zk8(avk8+2-1) = zk8(acpgd+2-1)
            zk8(avk8+3-1) = zk8(acpgd+4-1)
            zk8(avk8+4-1) = zk8(acpgd+5-1)
            zk8(avk8+5-1) = zk8(acpgd+6-1)
        else if ((option .eq. 'DEGE_ELNO') .or. (option .eq. &
                                                 'DEGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+i+6-1)
            end do
            iret = 0
        else if ((option .eq. 'EFGE_ELNO') .or. (option .eq. &
                                                 'EFGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+i+13-1)
            end do
            iret = 0
        else
            iret = 0
            valk(1) = option
            valk(2) = ' '
            valk(3) = 'SIGM_ELNO'
            valk(4) = 'SIEF_ELNO'
            valk(5) = 'EPSI_ELNO'
            valk(6) = 'EPSG_ELNO'
            valk(7) = 'EPME_ELNO'
            valk(8) = 'EPMG_ELNO'
            valk(9) = 'DEGE_ELNO'
            valk(10) = 'EFGE_ELNO'
            valk(11) = 'SIGM_NOEU'
            valk(12) = 'SIEF_NOEU'
            valk(13) = 'EPSI_NOEU'
            valk(14) = 'EPSG_NOEU'
            valk(15) = 'EPME_NOEU_DEPL'
            valk(16) = 'EPMG_NOEU'
            valk(17) = 'DEGE_NOEU'
            valk(18) = 'EFGE_NOEU'
            valk(19) = ' '
            valk(20) = 'FLUX_R'
            call utmess('F', 'POSTRELE_34', nk=20, valk=valk, si=iocc)
        end if
!
    else if ((ntd1 .ne. 0) .or. (ntd2 .ne. 0)) then
!
!      /* LA DIRECTION EST AU PLUS SUIVANT (X,Y,Z) */
!
        call getvr8(mcf, 'DIRECTION', iocc=iocc, nbval=3, vect=dir, &
                    nbret=ibid)
        if (abs(dir(1)) .lt. 1.0d-6) then
            dirx = .false.
        else
            dirx = .true.
        end if
        if (abs(dir(2)) .lt. 1.0d-6) then
            diry = .false.
        else
            diry = .true.
        end if
        if (abs(dir(3)) .lt. 1.0d-6) then
            dirz = .false.
        else
            dirz = .true.
        end if
        call wkvect(nomaux, 'V V I', 6, avicp)
        if (dirx .and. diry .and. dirz) then
            codir = 7
        else if ((.not. dirx) .and. diry .and. dirz) then
            codir = 6
        else if ((.not. diry) .and. dirx .and. dirz) then
            codir = 5
        else if ((.not. dirz) .and. dirx .and. diry) then
            codir = 4
        else if ((.not. dirx) .and. (.not. diry) .and. dirz) then
            codir = 3
        else if ((.not. dirx) .and. (.not. dirz) .and. diry) then
            codir = 2
        else if ((.not. diry) .and. (.not. dirz) .and. dirx) then
            codir = 1
        else
        end if
        pt = 1
        if ((option .eq. 'DEPL_NOEU_DEPL') .or. (option .eq. 'FORC_NOEU_FORC')) then
            if (dirx) then
                zi(avicp+pt-1) = 1
                pt = pt+1
            end if
            if (diry) then
                zi(avicp+pt-1) = 2
                pt = pt+1
            end if
            if (dirz) then
                zi(avicp+pt-1) = 3
                pt = pt+1
            end if
        else if ((option .eq. 'SIGM_ELNO') .or. (option .eq. 'SIEF_ELNO') .or. &
                 (option .eq. 'EPSI_ELNO') .or. (option .eq. 'EPSG_ELNO') .or. &
                 (option .eq. 'EPME_ELNO') .or. (option .eq. 'EPMG_ELNO') .or. &
                 (option .eq. 'SIGM_NOEU') .or. (option .eq. 'SIEF_NOEU') .or. &
                 (option .eq. 'EPSI_NOEU') .or. (option .eq. 'EPSG_NOEU') .or. &
                 (option .eq. 'EPME_NOEU_DEPL') .or. (option .eq. 'EPMG_NOEU')) then
            call wkvect(nomnew, 'V V I', 3, avinew)
            if (dirx) then
                zi(avinew+1-1) = 1
                zi(avinew+2-1) = 4
                zi(avinew+3-1) = 5
                call i2trgi(zi(avicp), zi(avinew), 3, pt)
            end if
            if (diry) then
                zi(avinew+1-1) = 2
                zi(avinew+2-1) = 4
                zi(avinew+3-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 3, pt)
            end if
            if (dirz) then
                zi(avinew+1-1) = 3
                zi(avinew+2-1) = 5
                zi(avinew+3-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 3, pt)
            end if
            call jedetr(nomnew)
        else if ((option .eq. 'EFGE_ELNO') .or. (option .eq. &
                                                 'EFGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+13+i-1)
            end do
            iret = 0
            call wkvect(nomnew, 'V V I', 4, avinew)
            if (dirx) then
                zi(avinew+1-1) = 1
                zi(avinew+2-1) = 3
                zi(avinew+3-1) = 4
                zi(avinew+4-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 4, pt)
            end if
            if (diry) then
                zi(avinew+1-1) = 2
                zi(avinew+2-1) = 3
                zi(avinew+3-1) = 5
                zi(avinew+4-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 4, pt)
            end if
            call jedetr(nomnew)
        else if ((option .eq. 'DEGE_ELNO') .or. (option .eq. &
                                                 'DEGE_NOEU')) then
            call wkvect(nomojb, 'V V K8', 6, avk8)
            do i = 1, 6, 1
                zk8(avk8+i-1) = zk8(acpgd+6+i-1)
            end do
            iret = 0
            call wkvect(nomnew, 'V V I', 4, avinew)
            if (dirx) then
                zi(avinew+1-1) = 1
                zi(avinew+2-1) = 3
                zi(avinew+3-1) = 4
                zi(avinew+4-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 4, pt)
            end if
            if (diry) then
                zi(avinew+1-1) = 2
                zi(avinew+2-1) = 3
                zi(avinew+3-1) = 5
                zi(avinew+4-1) = 6
                call i2trgi(zi(avicp), zi(avinew), 4, pt)
            end if
            call jedetr(nomnew)
        else
            iret = 0
            valk(1) = option
            valk(2) = ' '
            valk(3) = 'SIGM_ELNO'
            valk(4) = 'SIEF_ELNO'
            valk(5) = 'EPSI_ELNO'
            valk(6) = 'EPSG_ELNO'
            valk(7) = 'EPME_ELNO'
            valk(8) = 'EPMG_ELNO'
            valk(9) = 'DEGE_ELNO'
            valk(10) = 'EFGE_ELNO'
            valk(11) = 'DEGE_ELNO'
            valk(12) = 'SIGM_NOEU'
            valk(13) = 'SIEF_NOEU'
            valk(14) = 'EPSI_NOEU'
            valk(15) = 'EPSG_NOEU'
            valk(16) = 'EPME_NOEU_DEPL'
            valk(17) = 'EPMG_NOEU'
            valk(18) = 'DEGE_NOEU'
            valk(19) = 'EFGE_NOEU'
            valk(20) = 'DEGE_NOEU'
            valk(21) = ' '
            valk(22) = 'DEPL_R'
            valk(23) = 'FORC_R'
            call utmess('F', 'POSTRELE_35', nk=23, valk=valk, si=iocc)
        end if
        if (iret .ne. 0) then
            pt = pt-1
            if (pt .gt. 0) then
                call wkvect(nomojb, 'V V K8', pt, avk8)
                do i = 1, pt, 1
                    zk8(avk8+i-1) = zk8(acpgd+zi(avicp+i-1)-1)
                end do
            else if (iret .ne. 0) then
                iret = 0
                call utmess('F', 'POSTRELE_36', si=iocc)
            else
            end if
        end if
        call jeexin(nomaux, i)
        if (i .ne. 0) then
            call jedetr(nomaux)
        end if
!
    else if (nso .ne. 0) then
!
!     /* CAS SOMME (REPERE GLOBAL) */
!
        call wkvect(nomojb, 'V V K8', nbcpc, avk8)
        do i = 1, nbcpc, 1
            zk8(avk8+i-1) = zk8(alscpc+i-1)
        end do
    else
!
!     /* CAS NOM_CMP */
!
        call jeexin(nch19//'.DESC', ibid)
        call dismoi("DOCU", nch19, "CHAMP", repk=docu)
!
        if (repere .eq. 'GLOBAL') then
            call wkvect(nomojb, 'V V K8', nbcpc, avk8)
            do i = 1, nbcpc, 1
                zk8(avk8+i-1) = zk8(alscpc+i-1)
            end do
        else
!JMP       CALL WKVECT('&&RVCPNC.LISTE.IS','V V I',12,ALSI)
            call wkvect('&&RVCPNC.LISTE.IS', 'V V I', nbcpgd, alsi)
            if ((nomgd(1:6) .eq. 'SIEF_R') .or. (nomgd .eq. 'EPSI_R')) then
!          /* CHGT DE REPERE POUR SIGMA, EPSI, (N,M) OU (E,K) */
                call jelira(nlscpc, 'LONMAX', nn)
                AS_ALLOCATE(vk8=tmp, size=nbcpgd)
                do i = 1, nbcpgd
                    tmp(i) = ' '
                end do
                do i = 1, nn
                    tmp(i) = zk8(alscpc+i-1)
                end do
                call num2k8(nomgd, tmp, zk8(acpgd), nbcpgd, zi(alsi))
                AS_DEALLOCATE(vk8=tmp)
                do i = 1, 6, 1
                    n1 = n1+zi(alsi+i-1)
                end do
                do i = 1, 3, 1
                    n2 = n2+zi(alsi+6+i-1)
                    n3 = n3+zi(alsi+9+i-1)
                end do
                if (n1 .ne. 0) then
                    option = nomgd
                    if (docu .eq. 'CHNO') then
                        option(5:14) = '_NOEU_DEPL'
                    else
                        option(5:14) = '_ELNO_DEPL'
                    end if
                    call jedetr(nlscpc)
                    call wkvect(nlscpc, 'V V K8', n1, alscpc)
                    call wkvect('&&RVCPNC.TABIS.1', 'V V I', 6, iadt1)
                    call wkvect('&&RVCPNC.TABIS.2', 'V V I', 3, iadt2)
                    pt = 1
                    ptnc = 1
                    do i = 1, 6, 1
                        if (zi(alsi+i-1) .ne. 0) then
                            zk8(alscpc+pt-1) = zk8(acpgd+i-1)
                            pt = pt+1
                            if ((i .eq. 1) .or. (i .eq. 2) .or. (i .eq. 4)) then
                                zi(iadt2+1-1) = 1
                                zi(iadt2+2-1) = 2
                                zi(iadt2+3-1) = 4
                                call i2trgi(zi(iadt1), zi(iadt2), 3, ptnc)
                            else if ((i .eq. 5) .or. (i .eq. 6)) &
                                then
                                zi(iadt2+1-1) = 5
                                zi(iadt2+2-1) = 6
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else
                                zi(iadt2+1-1) = 3
                                call i2trgi(zi(iadt1), zi(iadt2), 1, ptnc)
                            end if
                        end if
                    end do
                    call wkvect(nomojb, 'V V K8', ptnc-1, avk8)
                    do i = 1, ptnc-1, 1
                        zk8(avk8+i-1) = zk8(acpgd+zi(iadt1+i-1)-1)
                    end do
                    call jedetr('&&RVCPNC.TABIS.1')
                    call jedetr('&&RVCPNC.TABIS.2')
                    call utmess('I', 'POSTRELE_37', si=iocc)
                else
                    if ((n2+n3) .ne. 0) then
                        if (nomgd .eq. 'SIEF_R') then
                            option(1:4) = 'EFGE'
                        end if
                        if (docu .eq. 'CHNO') then
                            option(5:14) = '_NOEU_DEPL'
                        else
                            option(5:14) = '_ELNO_DEPL'
                        end if
                        call wkvect(nomojb, 'V V K8', 6, avk8)
                        call jedetr(nlscpc)
                        call wkvect(nlscpc, 'V V K8', n2+n3, alscpc)
                        pt = 1
                        do i = 1, 6, 1
                            zk8(avk8+i-1) = zk8(acpgd+i+6-1)
                            if (zi(alsi+i+6-1) .ne. 0) then
                                zk8(alscpc+pt-1) = zk8(acpgd+i+6-1)
                                pt = pt+1
                            end if
                        end do
                    else
                        iret = 0
                        call utmess('F', 'POSTRELE_38', si=iocc)
                    end if
                end if
            else if ((nomgd(1:6) .eq. 'DEPL_R') .or. (nomgd(1:6) &
                                                      .eq. 'FORC_R')) then
                call numek8(zk8(alscpc), zk8(acpgd), nbcpc, nbcpgd, zi(alsi))
                do i = 1, 6, 1
                    n1 = n1+zi(alsi+i-1)
                end do
                if (n1 .ne. 0) then
                    call jedetr(nlscpc)
                    call wkvect(nlscpc, 'V V K8', n1, alscpc)
                    call wkvect('&&RVCPNC.TABIS.1', 'V V I', 6, iadt1)
                    call wkvect('&&RVCPNC.TABIS.2', 'V V I', 2, iadt2)
                    pt = 1
                    ptnc = 1
                    do i = 1, 6, 1
                        if (zi(alsi+i-1) .ne. 0) then
                            zk8(alscpc+pt-1) = zk8(acpgd+i-1)
                            pt = pt+1
                            if (i .eq. 1) then
                                zi(iadt2+1-1) = 1
                                zi(iadt2+2-1) = 2
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else if (i .eq. 2) then
                                zi(iadt2+1-1) = 1
                                zi(iadt2+2-1) = 2
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else if (i .eq. 3) then
                                zi(iadt2+1-1) = 3
                                call i2trgi(zi(iadt1), zi(iadt2), 1, ptnc)
                            else if (i .eq. 4) then
                                zi(iadt2+1-1) = 4
                                zi(iadt2+2-1) = 5
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else if (i .eq. 5) then
                                zi(iadt2+1-1) = 4
                                zi(iadt2+2-1) = 5
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else
                                zi(iadt2+1-1) = 6
                                call i2trgi(zi(iadt1), zi(iadt2), 1, ptnc)
                            end if
                        end if
                    end do
                    call wkvect(nomojb, 'V V K8', ptnc-1, avk8)
                    do i = 1, ptnc-1, 1
                        zk8(avk8+i-1) = zk8(acpgd+zi(iadt1+i-1)-1)
                    end do
                    call jedetr('&&RVCPNC.TABIS.1')
                    call jedetr('&&RVCPNC.TABIS.2')
                    if (nomgd(1:6) .eq. 'DEPL_R') then
                        option = 'DEPL_NOEU_DEPL'
                    else
                        option = 'FORC_NOEU_FORC'
                    end if
                    if (docu .eq. 'CHNO') then
                        option(5:9) = '_NOEU'
                    else
                        option(5:9) = '_ELNO'
                    end if
                else
                    iret = 0
                    call utmess('F', 'POSTRELE_38', si=iocc)
                end if
            else if (nomgd(1:6) .eq. 'FLUX_R') then
                call numek8(zk8(alscpc), zk8(acpgd), nbcpc, nbcpgd, zi(alsi))
                do i = 1, 3, 1
                    n1 = n1+zi(alsi+i-1)
                end do
                if (n1 .ne. 0) then
                    call jedetr(nlscpc)
                    call wkvect(nlscpc, 'V V K8', n1, alscpc)
                    call wkvect('&&RVCPNC.TABIS.1', 'V V I', 6, iadt1)
                    call wkvect('&&RVCPNC.TABIS.2', 'V V I', 2, iadt2)
                    pt = 1
                    ptnc = 1
                    do i = 1, 3, 1
                        if (zi(alsi+i-1) .ne. 0) then
                            zk8(alscpc+pt-1) = zk8(acpgd+i-1)
                            pt = pt+1
                            if (i .eq. 1) then
                                zi(iadt2+1-1) = 1
                                zi(iadt2+2-1) = 2
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else if (i .eq. 2) then
                                zi(iadt2+1-1) = 1
                                zi(iadt2+2-1) = 2
                                call i2trgi(zi(iadt1), zi(iadt2), 2, ptnc)
                            else
                                zi(iadt2+1-1) = 3
                                call i2trgi(zi(iadt1), zi(iadt2), 1, ptnc)
                            end if
                        end if
                    end do
                    call wkvect(nomojb, 'V V K8', ptnc-1, avk8)
                    do i = 1, ptnc-1, 1
                        zk8(avk8+i-1) = zk8(acpgd+zi(iadt1+i-1)-1)
                    end do
                    call jedetr('&&RVCPNC.TABIS.1')
                    call jedetr('&&RVCPNC.TABIS.2')
                    option = 'FLUX_ELNO'
                    if (docu .eq. 'CHNO') then
                        option(5:9) = '_NOEU'
                    end if
                else
                    iret = 0
                    call utmess('F', 'POSTRELE_38', si=iocc)
                end if
            else if (nomgd(1:6) .eq. 'TEMP_R') then
                call wkvect(nomojb, 'V V K8', 1, avk8)
                zk8(avk8+1-1) = zk8(acpgd+1-1)
                repere = 'GLOBAL'
            else
                iret = 0
                valk(1) = nomgd
                valk(2) = 'SIEF_R '
                valk(3) = 'SIGM_ELNO'
                valk(4) = 'EPSI_R '
                valk(5) = 'EP.._ELNO_DEPL'
                valk(6) = 'DEGE_ELNO'
                valk(7) = 'DEPL_R'
                valk(8) = 'FORC_R'
                call utmess('F', 'POSTRELE_39', nk=8, valk=valk, si=iocc)
            end if
            call jedetr('&&RVCPNC.LISTE.IS')
        end if
    end if
    call jedema()
end subroutine
