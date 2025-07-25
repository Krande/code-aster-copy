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

subroutine rvouex(mcf, iocc, nchpt, lstcmp, lstmac, &
                  lstnac, iret)
    implicit none
#include "jeveux.h"
#include "asterfort/celcel.h"
#include "asterfort/celver.h"
#include "asterfort/cncinv.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/i2trgi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/rvgnoe.h"
#include "asterfort/utmach.h"
#include "asterfort/utmess.h"
#include "asterfort/utncmp.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: lstcmp, lstmac, lstnac
    character(len=*) :: mcf, nchpt
    integer(kind=8) :: iocc, iret
!
!  OPERATION REALISEE
!  ------------------
!     CONSTRUCTION DES VECTEURS DES MAILLES ET/OU NOEUDS ACTIFS
!
!  ARGUMENTS EN ENTREE
!  -------------------
!     IOCC   : NUMERO DE L' OCCURENCE TRAITEE
!     NCHPT  : NOM DU CHAM_GD A TRAITER
!     LSTCMP : NOM DU VECTEUR DES NUMEROS DE CMP MISES EN JEU
!
!  ARGUMENTS EN SORTIE
!  -------------------
!     LSTMAC : NOM DE L' OJB 'V V I' DES NUMEROS DE MAILLES ACTIVES
!     LSTNAC : NOM DE L' OJB 'V V I' DES NUMEROS DE NOEUDS ACTIFS
!     IRET   : CODE RETOUR 1 = OK, 0 = KO
!
!  REMARQUE
!  --------
!     SUIVANT LA NATURE DU CHAMP TRAITE UN SEUL DES OJB LSTMAC ET
!     LSTNAC EST CONSTRUIT
!
!**********************************************************************
!
!
!
!  VARIABLES LOCALES
!  -----------------
    integer(kind=8) :: adr, acncin, alsmac, alsnac, acmp, adrvlc, arepe
    integer(kind=8) :: nbtma, nbm, nbmac, nbnac, nbmalu
    integer(kind=8) :: i, in, n, m, libre, n1, ibid, igrel, jnuma, j
    integer(kind=8) :: ibib, imolo, n2, kk, ier, nbvari, nbr
    integer(kind=8) :: ii, jmmail, nbtrou, nbcmp, nbcmp1, nc, jcmp, jcmp1, ntc
    character(len=4) :: docu
    character(len=8) :: nmaila, nomgd, resuco, nomvar, num
    character(len=15) :: nconec
    character(len=16) :: motcle(2), typmcl(2), nchsym
    character(len=19) :: nchp19
    character(len=24) :: ncncin, nrepe, lismai, malist, nomobj, valk(3)
    integer(kind=8), pointer :: celd(:) => null()
    integer(kind=8), pointer :: entier(:) => null()
    data nbvari/100/
!**********************************************************************
!
    call jemarq()
!
    call jeveuo(jexnum(lstcmp, iocc), 'L', acmp)
    malist = '&&RVOUEX_MALIST'
!
!
    nchp19 = nchpt
    iret = 1
!
    if (nchp19(1:1) .ne. '&') then
!
        call dismoi('DOCU', nchp19, 'CHAMP', repk=docu)

        if (docu == "CHML") then
!
!          -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
!
            call celcel('NBVARI_CST', nchp19, 'V', '&&RVOUEX.CHAMEL1')
            nchp19 = '&&RVOUEX.CHAMEL1'
            call celver(nchp19, 'NBSPT_1', 'COOL', kk)
            if (kk .eq. 1) then
                call dismoi('NOM_GD', nchp19, 'CHAMP', repk=nomgd)
                call utmess('I', 'PREPOST_36', sk=nomgd)
                call celcel('PAS_DE_SP', nchp19, 'V', '&&RVOUEX.CHAMEL2')
                nchp19 = '&&RVOUEX.CHAMEL2'
            end if
            call jelira(nchp19//'.CELD', 'DOCU', cval=docu)
        end if
!
        call dismoi('NOM_MAILLA', nchp19, 'CHAMP', repk=nmaila)
        nconec = nmaila//'.CONNEX'
        ncncin = '&&OP0051.CONNECINVERSE  '
!
        call jelira(nconec, 'NMAXOC', nbtma)
!
        nbtrou = 0
        jmmail = 1
!
        call getvtx(mcf, 'NOM_CMP', iocc=iocc, nbval=0, nbret=nc)
        if (nc .lt. 0) then
            nbcmp = -nc
            call wkvect('&&RVOUEX.NOM_CMP', 'V V K8', nbcmp, jcmp)
            call getvtx(mcf, 'NOM_CMP', iocc=iocc, nbval=nbcmp, vect=zk8(jcmp), &
                        nbret=nc)
!
! VERIFICATION QUE LES COMPOSANTES DEMANDEES
! APPARTIENNENT BIEN AU CHAMP
!
            call getvid(mcf, 'RESULTAT', iocc=iocc, scal=resuco, nbret=ibib)
            if (ibib .ne. 0) then
                nomobj = '&&RVOUEX.NOM_CMP1'
                call jeexin(nomobj, ier)
                if (ier .ne. 0) call jedetr(nomobj)
                call utncmp(nchp19, nbcmp1, nomobj)
                call jeveuo(nomobj, 'L', jcmp1)
                call getvtx(mcf, 'NOM_CHAM', iocc=iocc, scal=nchsym, nbret=n1)
                nbr = nbcmp1
                if (zk8(jcmp1) .eq. 'VARI') nbr = nbvari
                do i = 1, nbcmp
                    do j = 1, nbr
                        if (zk8(jcmp1) .eq. 'VARI') then
                            call codent(j, 'G', num)
                            nomvar = 'V'//num(1:7)
                            if (zk8(jcmp-1+i) .eq. nomvar) goto 102
                        else
                            if (zk8(jcmp-1+i) .eq. zk8(jcmp1-1+j)) goto 102
                        end if
                    end do
                    valk(1) = zk8(jcmp-1+i)
                    valk(2) = nchsym
                    valk(3) = resuco
                    call utmess('F', 'POSTRELE_65', nk=3, valk=valk)
102                 continue
                end do
            end if
!
            call utmach(nchp19, nbcmp, zk8(jcmp), 'NU', malist, &
                        nbtrou)
            if (nbtrou .ne. 0) call jeveuo(malist, 'L', jmmail)
            call jedetr('&&RVOUEX.NOM_CMP')
        end if
!
        call getvtx(mcf, 'TOUT_CMP', iocc=iocc, nbval=0, nbret=ntc)
        if (ntc .lt. 0) then
            nomobj = '&&RVOUEX.NOMCMP.USER'
            call utncmp(nchp19, nbcmp, nomobj)
            call jeveuo(nomobj, 'L', jcmp)
            call utmach(nchp19, nbcmp, zk8(jcmp), 'NU', malist, &
                        nbtrou)
            if (nbtrou .ne. 0) call jeveuo(malist, 'L', jmmail)
            call jedetr(nomobj)
        end if
!
        if (docu .eq. 'CHML') then
!             ----------------
            call jeveuo(nchp19//'.CELK', 'L', adr)
            nrepe = zk24(adr) (1:19)//'.REPE'
            call jeveuo(nrepe, 'L', arepe)
!
            ibid = 0
            call rvgnoe(mcf, iocc, nmaila, lstnac, 0, &
                        [ibid])
!
            call getvtx(mcf, 'GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
            call getvtx(mcf, 'MAILLE', iocc=iocc, nbval=0, nbret=n2)
            if ((n1+n2) .eq. 0) then
                nbmalu = 0
            else
                lismai = '&&RVOUEX.NUME_MAIL'
                motcle(1) = 'GROUP_MA'
                motcle(2) = 'MAILLE'
                typmcl(1) = 'GROUP_MA'
                typmcl(2) = 'MAILLE'
                call reliem(' ', nmaila, 'NU_MAILLE', mcf, iocc, &
                            2, motcle, typmcl, lismai, nbmalu)
                call jeveuo(lismai, 'L', jnuma)
            end if
!
            call jeexin(ncncin, n2)
            if (n2 .eq. 0) call cncinv(nmaila, [ibid], 0, 'V', ncncin)
!
            call jelira(lstnac, 'LONMAX', nbnac)
            call jeveuo(lstnac, 'L', alsnac)
!
            call jecreo('&&RVOUEX.LISTE.ENTIER', 'V V I')
            call jeecra('&&RVOUEX.LISTE.ENTIER', 'LONMAX', nbtma)
            call jeveuo('&&RVOUEX.LISTE.ENTIER', 'E', vi=entier)
!
            libre = 1
            call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', adrvlc)
            call jeveuo(jexnum(ncncin, 1), 'L', acncin)
!
            do in = 1, nbnac, 1
                n = zi(alsnac+in-1)
                nbm = zi(adrvlc+n+1-1)-zi(adrvlc+n-1)
                adr = zi(adrvlc+n-1)
!
                call i2trgi(entier, zi(acncin+adr-1), nbm, libre)
!
            end do
!
            nbmac = libre-1
            libre = 1
!
            call jeveuo(nchp19//'.CELD', 'L', vi=celd)
!
            do i = 1, nbmac, 1
                m = entier(i)
                if (nbtrou .ne. 0) then
                    do ii = 1, nbtrou
                        if (m .eq. zi(jmmail+ii-1)) goto 114
                    end do
                    goto 110
114                 continue
                end if
                if (m .ne. 0) then
                    if (nbmalu .ne. 0) then
                        do j = 1, nbmalu, 1
                            if (m .eq. zi(jnuma+j-1)) goto 404
                        end do
                        goto 110
404                     continue
                    end if
                    igrel = zi(arepe+2*(m-1))
                    if (igrel .gt. 0) then
                        imolo = celd(celd(4+igrel)+2)
                        if (imolo .gt. 0) then
                            entier(libre) = entier(i)
                            libre = libre+1
                        end if
                    end if
                end if
110             continue
            end do
!
            nbmac = libre-1
!
            if (nbmac .gt. 0) then
!
                call wkvect(lstmac, 'V V I', nbmac, alsmac)
!
                do i = 1, nbmac, 1
                    zi(alsmac+i-1) = entier(i)
                end do
!
            else
!
                iret = 0
!
            end if
!
            call jedetr('&&RVOUEX.LISTE.ENTIER')
            call jedetr('&&RVOUEX.NUME_MAIL')
!
        else
!             ----------------
!
!
            call rvgnoe(mcf, iocc, nmaila, lstnac, nbtrou, &
                        zi(jmmail))
!
        end if
!
    end if
!
    call jedetr(malist)
    call detrsd('CHAM_ELEM', '&&RVOUEX.CHAMEL1')
    call detrsd('CHAM_ELEM', '&&RVOUEX.CHAMEL2')
    call jedema()
end subroutine
