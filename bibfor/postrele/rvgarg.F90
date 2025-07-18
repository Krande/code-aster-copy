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

subroutine rvgarg(nxdnom, nxdnum, nvchef, nvcodo, nxdvar)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/iunifi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeimpo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/numek8.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utcmp2.h"
#include "asterfort/utmess.h"
#include "asterfort/utncmp.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: nxdnom, nxdnum, nvchef, nvcodo, nxdvar
!
!      SAISIE ET PREPA VERIF COHERENCE DES ARGUMENTS D 'APPEL DE OP0051
!          1. VERIFICATION D' EXISTENCE DES OJB MIS EN JEU
!          2. VERIFICATION DE LEGALITE DES CMP MISES EN JEU
!          3. VERIFICATION DE COHERENCE SUR LES MAILLAGES SOUS-JACENTS
!     CES VERIFICATIONS SERONT EFFECTUEES PAR RVCOHE
!     SAISIE DES CMP MISES EN JEU ET DES OPERATION DEMANDEE
!     ------------------------------------------------------------------
! OUT NXDNOM : K : XD V K8 NOM DES CMP UTILISATEURS (OC(I) = OCCUR(I))
! OUT NXDNUM : K : XD V K8 NUM DES CMP UTILISATEURS (OC(I) = OCCUR(I))
! OUT NVCHEF : K : S V K24 NOM D' UN CHAMP EFFECTIFS(TYPE CHAM_GD)
! OUT NVCODO : K : S V I   CODE OPERATION POST
!     ------------------------------------------------------------------
!
    character(len=80) :: text80, text1
    character(len=24) :: naux24, kordre, nomobj, operation
    character(len=19) :: nchp19
    character(len=16) :: nchsym
    character(len=8) :: k8b, nresu, nchgd, granch, nomcp(50), nomail
    character(len=4) :: typech
    aster_logical :: existe, parMesh
    integer(kind=8) :: anomcp, anumcp, ancpu1, ancpu2, acpgd, avchef
    integer(kind=8) :: n1, n2, i, iocc, gd, n3, adr, nbelp, nbinv, ibid, avcodo
    integer(kind=8) :: nbpost, nbchgd, nbcpgd, nbcmp, nbresu, nbtcp, nbsom
    integer(kind=8) :: ifr, j, jordr, jxvar, n4, nbc, nbnc, numecp(50), tord(1)
    real(kind=8) :: rbid
    complex(kind=8) :: cbid
!
!=================== CORPS DE LA ROUTINE =============================
!
    call jemarq()
    ifr = iunifi('RESULTAT')
!
    rbid = 1.0d-6
    cbid = dcmplx(rbid, rbid)
    call getfac('ACTION', nbpost)
    call jecrec(nxdnom, 'V V K8', 'NU', 'DISPERSE', 'VARIABLE', &
                nbpost)
    call jecrec(nxdvar, 'V V I ', 'NU', 'DISPERSE', 'VARIABLE', &
                nbpost)
    call jecrec(nxdnum, 'V V I ', 'NU', 'DISPERSE', 'VARIABLE', &
                nbpost)
    call wkvect(nvchef, 'V V K24', nbpost, avchef)
    call wkvect(nvcodo, 'V V I', nbpost, avcodo)
    kordre = '&&RVGARG.NUMEORDR'
    do iocc = 1, nbpost, 1
        n1 = 0
        n2 = 0
        call getvtx('ACTION', 'OPERATION', iocc=iocc, nbval=0, nbret=n3)
        n3 = -n3
        call wkvect('&&RVGARG.NOM.OPERATION', 'V V K80', n3, adr)
        call getvtx('ACTION', 'OPERATION', iocc=iocc, nbval=n3, vect=zk80(adr), &
                    nbret=n4)
        operation = zk80(adr)
        if (n3 .eq. 1) then
            text1 = zk80(adr+1-1)
            if (text1(1:1) .eq. 'E') then
                zi(avcodo+iocc-1) = 1
            else
                zi(avcodo+iocc-1) = 3
            end if
        else
            zi(avcodo+iocc-1) = 2
        end if
        call jedetr('&&RVGARG.NOM.OPERATION')
        call getvid('ACTION', 'RESULTAT', iocc=iocc, nbval=0, nbret=nbresu)
        call getvid('ACTION', 'CHAM_GD', iocc=iocc, nbval=0, nbret=nbchgd)
        nbresu = -nbresu
        nbchgd = -nbchgd
        if (nbresu .ne. 0) then
!        /* CAS D' UN RESULTAT COMPOSE */
            call getvid('ACTION', 'RESULTAT', iocc=iocc, scal=nresu, nbret=n1)
            call getvtx('ACTION', 'NOM_CHAM', iocc=iocc, scal=text80, nbret=n1)
            call dismoi('NOM_MAILLA', nresu, 'RESULTAT', repk=nomail)
            parMesh = isParallelMesh(nomail)
            nchsym = text80(1:16)
            call jenonu(jexnom(nresu//'           .DESC', nchsym), n1)
            if (n1 .ne. 0) then
!           /* LE CHAMP SYMBOLIQUE EXISTE (POTENTIELLEMENT)*/
                call rsorac(nresu, 'LONUTI', 0, rbid, k8b, &
                            cbid, rbid, 'RELATIF', tord, 1, &
                            ibid)
                n3 = tord(1)
                if (n3 .gt. 0) then
                    call wkvect(kordre, 'V V I', n3, jordr)
                    call rsorac(nresu, 'TOUT_ORDRE', 0, rbid, k8b, &
                                cbid, rbid, 'RELATIF', zi(jordr), n3, &
                                ibid)
                    do j = 1, n3
                        call rsexch(' ', nresu, nchsym, zi(jordr+j-1), naux24, &
                                    n2)
                        if (n2 .eq. 0) exit
                    end do
                    call jedetr(kordre)
                else
                    n2 = 1
                end if
            else
!           /* LE CHAMP SYMBOLIQUE N' EXISTE PAS */
                n2 = 1
                write (ifr, *) 'CHAMP SYMBOLIQUE >', nchsym, '< NON '// &
                    'AUTORISE POUR LE RESULTAT >', nresu, '<'
                write (ifr, *) 'LES CHAMPS SYMBOLIQUES AUTORISES SONT :'
                call jeimpo(ifr, nresu//'           .DESC', ' ')
            end if
            if ((n1 .eq. 0) .or. (n2 .ne. 0)) then
!           /* ALTERNATIVE :                              */
!           /* LE CHAMPS SYMBOLIQUE EST ILEGAL OU         */
!           /* AUCUN CHAMP EFFECTIF ASSOCIE N' A ETE CREE */
                existe = .false.
            else
                existe = .true.
                nchp19 = naux24(1:19)
            end if
        else
!        /* CAS D'UN CHAMP_GD */
            call getvid('ACTION', 'CHAM_GD', iocc=iocc, scal=nchgd, nbret=n1)
            existe = .true.
            nchp19 = nchgd//'           '
            call dismoi('NOM_MAILLA', nchgd, 'CHAMP', repk=nomail)
            parMesh = isParallelMesh(nomail)
        end if
        if (parMesh .and. operation .ne. 'EXTRACTION') then
            call utmess('F', 'POSTRELE_24')
        end if
        call jecroc(jexnum(nxdnom, iocc))
        call jecroc(jexnum(nxdnum, iocc))
!
        if (.not. existe) then
!        /* LE CHAMPS SYMBOLIQUE EST ILEGAL, OU                 */
!        /* IL EST LEGAL, MAIS IL N' ADMET AUCUN CHAMP EFFECTIF */
            call jeecra(jexnum(nxdnom, iocc), 'LONMAX', 1)
            call jeveuo(jexnum(nxdnom, iocc), 'E', anomcp)
            zk8(anomcp) = '&NOEXIST'
            call jeecra(jexnum(nxdnum, iocc), 'LONMAX', 1)
            call jeveuo(jexnum(nxdnum, iocc), 'E', anumcp)
            zi(anumcp) = 0
            zk24(avchef+iocc-1) = '&NONEXISTEOUNONCREE     '
        else
            zk24(avchef+iocc-1) = nchp19//'     '
            call dismoi('TYPE_CHAMP', nchp19, 'CHAMP', repk=typech)
            call dismoi('NOM_GD', nchp19, 'CHAMP', repk=granch)
            call dismoi('NUM_GD', nchp19, 'CHAMP', repi=gd)
            call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', nbcpgd)
            call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', acpgd)
            call getvtx('ACTION', 'NOM_CMP', iocc=iocc, nbval=0, nbret=nbcmp)
            call getvtx('ACTION', 'TOUT_CMP', iocc=iocc, nbval=0, nbret=nbtcp)
            call getvtx('ACTION', 'INVARIANT', iocc=iocc, nbval=0, nbret=nbinv)
            call getvtx('ACTION', 'ELEM_PRINCIPAUX', iocc=iocc, nbval=0, nbret=nbelp)
            call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=0, nbret=nbsom)
            nbcmp = -nbcmp
            nbtcp = -nbtcp
            nbinv = -nbinv
            nbelp = -nbelp
            nbsom = -nbsom
            if ((nbcmp .ne. 0) .or. (nbsom .ne. 0)) then
!           /* PASSAGE D' UNE OU DEUX LISTE DE NOM DE CMPS    */
!           /* MOT-CLE (NOM_CMP) OU (RESULTANTE ET/OU MOMENT) */
                if (nbcmp .ne. 0) then
                    call wkvect('&&OP0051.NOMCMP.USER', 'V V K8', nbcmp, ancpu1)
                    call getvtx('ACTION', 'NOM_CMP', iocc=iocc, nbval=nbcmp, vect=zk8(ancpu1), &
                                nbret=n1)
                else
                    if (typech .eq. 'ELNO' .and. granch .eq. 'VARI_R') then
                        call utmess('F', 'POSTRELE_20')
                    end if
                    call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=0, nbret=n1)
                    call getvtx('ACTION', 'MOMENT', iocc=iocc, nbval=0, nbret=n2)
                    n1 = -n1
                    n2 = -n2
                    nbcmp = n1+n2
                    call wkvect('&&OP0051.NOMCMP.USER', 'V V K8', nbcmp, ancpu1)
                    call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=n1, vect=zk8(ancpu1))
                    call getvtx('ACTION', 'MOMENT', iocc=iocc, nbval=n2, vect=zk8(ancpu1+n1))
                end if
                if (typech .eq. 'ELNO' .and. granch .eq. 'VARI_R') then
                    call utcmp2(granch, 'ACTION', iocc, 50, nomcp, &
                                numecp, nbnc)
                    call jeecra(jexnum(nxdvar, iocc), 'LONMAX', nbnc)
                    call jeecra(jexnum(nxdvar, iocc), 'LONUTI', nbnc)
                    call jeveuo(jexnum(nxdvar, iocc), 'E', jxvar)
                    do i = 1, nbnc
                        zi(jxvar+i-1) = numecp(i)
                    end do
                    nbcmp = 1
                    zk8(ancpu1) = 'VARI'
                else
                    call jeecra(jexnum(nxdvar, iocc), 'LONMAX', nbcmp)
                    call jeecra(jexnum(nxdvar, iocc), 'LONUTI', 0)
                end if
                call jeecra(jexnum(nxdnom, iocc), 'LONMAX', nbcmp)
                call jeveuo(jexnum(nxdnom, iocc), 'E', anomcp)
                call jeecra(jexnum(nxdnum, iocc), 'LONMAX', nbcmp)
                call jeveuo(jexnum(nxdnum, iocc), 'E', anumcp)
                do i = 1, nbcmp, 1
                    zk8(anomcp+i-1) = zk8(ancpu1+i-1)
                end do
                call numek8(zk8(acpgd), zk8(anomcp), nbcpgd, nbcmp, zi(anumcp))
                call jedetr('&&OP0051.NOMCMP.USER')
            else if (nbtcp .ne. 0) then
!
                nomobj = '&&OP0051.NOMCMP.USER'
                call utncmp(nchp19, nbc, nomobj)
                if (nbc .eq. 0) then
                    call utmess('F', 'POSTRELE_54', si=iocc)
                end if
                call jeveuo(nomobj, 'L', ancpu2)
                call jeecra(jexnum(nxdnom, iocc), 'LONMAX', nbc)
                call jeveuo(jexnum(nxdnom, iocc), 'E', anomcp)
                call jeecra(jexnum(nxdnum, iocc), 'LONMAX', nbc)
                call jeveuo(jexnum(nxdnum, iocc), 'E', anumcp)
                do i = 1, nbc, 1
                    zk8(anomcp+i-1) = zk8(ancpu2+i-1)
                end do
                call numek8(zk8(acpgd), zk8(anomcp), nbcpgd, nbc, zi(anumcp))
                call jedetr(nomobj)
                if (typech .eq. 'ELNO' .and. granch .eq. 'VARI_R') then
                    call jeecra(jexnum(nxdvar, iocc), 'LONMAX', nbc)
                    call jeecra(jexnum(nxdvar, iocc), 'LONUTI', nbc)
                    call jeveuo(jexnum(nxdvar, iocc), 'E', jxvar)
                    zi(jxvar) = -1
                else
                    call jeecra(jexnum(nxdvar, iocc), 'LONUTI', 0)
                end if
            else
!           /* PASSAGE DE CMPS IMPLICITES */
                call jeecra(jexnum(nxdnom, iocc), 'LONMAX', 1)
                call jeveuo(jexnum(nxdnom, iocc), 'E', anomcp)
                call jeecra(jexnum(nxdnum, iocc), 'LONMAX', 1)
                call jeveuo(jexnum(nxdnum, iocc), 'E', anumcp)
                zk8(anomcp) = 'IMPLICIT'
                zi(anumcp) = -1
            end if
        end if
    end do
    call jedema()
end subroutine
