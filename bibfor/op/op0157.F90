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
!
subroutine op0157()
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/irgene.h"
#include "asterfort/irtitr.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsinfo.h"
#include "asterfort/rstran.h"
#include "asterfort/rsutnu.h"
#include "asterfort/ulexis.h"
#include "asterfort/ulopen.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! IMPR_GENE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=3) :: toucha, toucmp, toupar, interp
    character(len=4) :: motfac
    character(len=8) :: k8b, form, meshName
    character(len=16) :: nomcmd, typcon, crit, fich
    character(len=19) :: gene, knum, kdisc, krang
    character(len=80) :: titre
    aster_logical :: lhist, lResu, lField
    integer(kind=8) :: ifi, iocc, ire2, iret, isy, jcmpg
    integer(kind=8) :: jdisc, jnosy, jordr, jpara, jrang, n, n01
    integer(kind=8) :: n10, n11, n21, n22, nbcmpg, nbdisc, nbnosy
    integer(kind=8) :: nbordr, nbpara, nc, nocc, np, nr
    real(kind=8) :: prec
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call getres(k8b, k8b, nomcmd)
!
!     --- FORMAT ---
!
    call getvtx(' ', 'FORMAT', scal=form, nbret=n)
!
!     --- FICHIER ---
!
    ifi = 0
    fich = ' '
    call getvis(' ', 'UNITE', scal=ifi, nbret=n)
    if (.not. ulexis(ifi)) then
        call ulopen(ifi, ' ', fich, 'NEW', 'O')
    end if
!
    motfac = 'GENE'
    call getfac(motfac, nocc)
!
    do iocc = 1, nocc
!
!        --- SEPARATION DES DIFFERENTES OCCURENCES---
!
        if (form .eq. 'RESULTAT') write (ifi, '(/,1X,80(''-''))')
!
        call getvid(motfac, 'RESU_GENE', iocc=iocc, scal=gene, nbret=nr)
        call gettco(gene, typcon)
!
!        --- ECRITURE DU TITRE ---
!
        lResu = ASTER_TRUE
        lField = ASTER_FALSE
        meshName = ' '
        call irtitr(lResu, lField, gene, meshName, form, ifi, titre)
!
!        --- IMPRESSION DE LA STRUCTURE DU RESU_GENE ---
!
        if (typcon .eq. 'MODE_GENE') then
            call getvtx(motfac, 'INFO_GENE', iocc=iocc, scal=k8b, nbret=n01)
            if (k8b(1:3) .eq. 'OUI') call rsinfo(gene, ifi)
        end if
!
!        --- QUELS SONT LES NOM_CHAMP A IMPRIMER ---
!
        toucha = 'OUI'
        call getvtx(motfac, 'TOUT_CHAM', iocc=iocc, scal=toucha, nbret=n21)
        call getvtx(motfac, 'NOM_CHAM', iocc=iocc, nbval=0, nbret=n22)
        if (n22 .lt. 0) then
            nbnosy = -n22
            call wkvect('&&OP0157.NOM_SYMB', 'V V K16', nbnosy, jnosy)
            call getvtx(motfac, 'NOM_CHAM', iocc=iocc, nbval=nbnosy, vect=zk16(jnosy), &
                        nbret=n)
        else if (toucha .eq. 'OUI') then
            if (typcon .eq. 'MODE_GENE') then
                call jelira(gene//'.DESC', 'NOMUTI', nbnosy)
                call wkvect('&&OP0157.NOM_SYMB', 'V V K16', nbnosy, jnosy)
                do isy = 1, nbnosy
                    call jenuno(jexnum(gene//'.DESC', isy), zk16(jnosy-1+isy))
                end do
            else
                nbnosy = 3
                call wkvect('&&OP0157.NOM_SYMB', 'V V K16', nbnosy, jnosy)
                zk16(jnosy+1-1) = 'DEPL'
                zk16(jnosy+2-1) = 'VITE'
                zk16(jnosy+3-1) = 'ACCE'
            end if
        else if (toucha .eq. 'NON') then
            nbnosy = 0
            jnosy = 1
        end if
!
!        --- QUELS SONT LES CMP_GENE A IMPRIMER ---
!
        nbcmpg = -1
        jcmpg = 1
        toucmp = '   '
        call getvtx(motfac, 'TOUT_CMP_GENE', iocc=iocc, scal=toucmp, nbret=n21)
        call getvis(motfac, 'NUME_CMP_GENE', iocc=iocc, nbval=0, nbret=n22)
        if (toucmp .eq. 'NON') then
            nbcmpg = 0
        else if (n22 .lt. 0) then
            nbcmpg = -n22
            call wkvect('&&OP0157.NOM_CMPG', 'V V I', nbcmpg, jcmpg)
            call getvis(motfac, 'NUME_CMP_GENE', iocc=iocc, nbval=nbcmpg, vect=zi(jcmpg), &
                        nbret=n)
        end if
!
!        --- ON RECHERCHE LES PARAMETRES A ECRIRE ---
!
        nbpara = -1
        jpara = 1
        toupar = '   '
        call getvtx(motfac, 'TOUT_PARA', iocc=iocc, scal=toupar, nbret=n11)
        call getvtx(motfac, 'NOM_PARA', iocc=iocc, nbval=0, nbret=n10)
        if (toupar .eq. 'NON') then
            nbpara = 0
        else if (n10 .ne. 0) then
            nbpara = -n10
            call wkvect('&&OP0157.NOMUTI_PARA', 'V V K16', nbpara, jpara)
            call getvtx(motfac, 'NOM_PARA', iocc=iocc, nbval=nbpara, vect=zk16(jpara), &
                        nbret=n)
        end if
!
!        --- LES ACCES ---
!
        nbordr = 0
        jordr = 1
        nbdisc = 0
        jdisc = 1
        jrang = 1
        if (typcon .eq. 'MODE_GENE') then
            knum = '&&OP0157.NUME_ORDRE'
            call getvr8(motfac, 'PRECISION', iocc=iocc, scal=prec, nbret=np)
            call getvtx(motfac, 'CRITERE', iocc=iocc, scal=crit, nbret=nc)
            call rsutnu(gene, motfac, iocc, knum, nbordr, &
                        prec, crit, iret)
            if (iret .ne. 0) goto 12
            call jeveuo(knum, 'L', jordr)
        elseif ((typcon .eq. 'HARM_GENE') .or. (typcon .eq. 'TRAN_GENE')) &
            then
            kdisc = '&&OP0157.DISCRET'
            krang = '&&OP0157.NUME_ORDRE'
            interp = 'NON'
            call rstran(interp, gene, motfac, iocc, kdisc, &
                        krang, nbdisc, iret)
            if (iret .ne. 0) goto 12
            call jeexin(kdisc, ire2)
            if (ire2 .gt. 0) then
                call jeveuo(kdisc, 'E', jdisc)
                call jeveuo(krang, 'E', jrang)
            end if
        end if
!
!        --- HISTOR ---
!
        lhist = .true.
        call getvtx(motfac, 'INFO_CMP_GENE', iocc=iocc, scal=k8b, nbret=n)
        if (k8b(1:3) .eq. 'NON') lhist = .false.
!
        call irgene(iocc, gene, form, ifi, nbnosy, &
                    zk16(jnosy), nbcmpg, zi(jcmpg), nbpara, zk16(jpara), &
                    nbordr, zi(jordr), nbdisc, zr(jdisc), zi(jrang), &
                    lhist)
!
12      continue
        call jedetr('&&OP0157.NOM_SYMB')
        call jedetr('&&OP0157.NOM_CMPG')
        call jedetr('&&OP0157.NOMUTI_PARA')
        call jedetr('&&OP0157.NUME_ORDRE')
        call jedetr('&&OP0157.DISCRET')
    end do
    call ulopen(-ifi, ' ', fich, ' ', 'O')
!
    call jedema()
end subroutine
