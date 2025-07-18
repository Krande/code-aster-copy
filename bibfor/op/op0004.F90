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
subroutine op0004()
    implicit none
!     OPERATEUR DEFI_NAPPE
!     STOCKAGE DANS UN OBJET DE TYPE FONCTION
!     -----------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/foimpr.h"
#include "asterfort/foston.h"
#include "asterfort/foverf.h"
#include "asterfort/fovern.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ordon1.h"
#include "asterfort/ordonn.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/uttrif.h"
#include "asterfort/wkvect.h"
    character(len=4) :: interp(2)
    character(len=24) :: valk
    character(len=8) :: k8b
    character(len=16) :: nomcmd, typfon, verif
    character(len=19) :: nomfon
    aster_logical :: defonc
    integer(kind=8) :: iret, iret2
    integer(kind=8) :: vali(2)
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifm, ifonc, iocc, ival, jval, l
    integer(kind=8) :: l1, ladrf, lnomf, lpar, lpar2, lpara, lparc
    integer(kind=8) :: lpro, lval, mxva, n, n1, n2, n3
    integer(kind=8) :: nbcoup, nbfonc, nbpara, nbval, niv, nprol, nv
!
!-----------------------------------------------------------------------
    call jemarq()
!
    verif = ' '
    defonc = .false.
    call getres(nomfon, typfon, nomcmd)
    call getvtx(' ', 'VERIF', scal=verif, nbret=n1)
    call getvr8(' ', 'PARA', nbval=0, nbret=n2)
    call getvid(' ', 'FONCTION', nbval=0, nbret=n3)
    nbpara = abs(n2)
    if (n3 .ne. 0) then
        nbfonc = -n3
    else
        defonc = .true.
        call getfac('DEFI_FONCTION', nbfonc)
    end if
!
    if (nbpara .ne. nbfonc) then
        vali(1) = nbpara
        vali(2) = nbfonc
        call utmess('F', 'UTILITAI8_3', ni=2, vali=vali)
    end if
!
!     --- VERIFICATION DE LA CROISSANCE DES PARAMETRES ---
    if (verif .eq. 'CROISSANT') then
        call wkvect('&&OP0004.TEMP.PARA', 'V V R', nbpara, lparc)
        call getvr8(' ', 'PARA', nbval=nbpara, vect=zr(lparc), nbret=n)
!        VERIF QUE LES PARA SONT STRICT CROISSANTS
        iret = 2
        call foverf(zr(lparc), nbpara, iret)
        if (iret .ne. 2) then
            call utmess('F', 'UTILITAI2_72')
        end if
        call jedetr('&&OP0004.TEMP.PARA')
    end if
!
    if (defonc) then
        do iocc = 1, nbfonc
            call getvr8('DEFI_FONCTION', 'VALE', iocc=iocc, nbval=0, nbret=nv)
            nv = -nv
            if (mod(nv, 2) .ne. 0) then
                vali(1) = iocc
                call utmess('F', 'UTILITAI8_4', si=vali(1))
            end if
            if (verif .eq. 'CROISSANT') then
                nbcoup = nv/2
                call wkvect('&&OP0004.TEMP.PARA', 'V V R', nv, lpara)
                call wkvect('&&OP0004.TEMP.PAR2', 'V V R', nbcoup, lpar2)
                call getvr8('DEFI_FONCTION', 'VALE', iocc=iocc, nbval=nv, vect=zr(lpara), &
                            nbret=nbval)
                do i = 0, nbcoup-1
                    zr(lpar2+i) = zr(lpara+2*i)
                end do
!              VERIF QUE LES PARA SONT STRICT CROISSANTS
                iret = 2
                call foverf(zr(lpar2), nbcoup, iret)
                if (iret .ne. 2) then
                    call utmess('F', 'UTILITAI2_72')
                end if
                call jedetr('&&OP0004.TEMP.PARA')
                call jedetr('&&OP0004.TEMP.PAR2')
            end if
        end do
    end if
!
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infmaj()
    call infniv(ifm, niv)
!
!   --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON.PROL ---
    nprol = 7+2*nbfonc
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', nprol, lpro)
    zk24(lpro) = 'NAPPE   '
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=l1)
    if (l1 .eq. 1) interp(2) = interp(1)
    zk24(lpro+1) = interp(1)//interp(2)
    call getvtx(' ', 'NOM_PARA', scal=zk24(lpro+2), nbret=l)
    call getvtx(' ', 'NOM_RESU', scal=zk24(lpro+3), nbret=l)
    call getvtx(' ', 'PROL_GAUCHE', scal=zk24(lpro+4) (1:1), nbret=l)
    call getvtx(' ', 'PROL_DROITE', scal=zk24(lpro+4) (2:2), nbret=l)
    zk24(lpro+5) = nomfon
!
!   --- CREATION ET REMPLISSAGE DE L'OBJET NOMFON.PARA ---
    call wkvect(nomfon//'.PARA', 'G V R', nbpara, lpar)
    call getvr8(' ', 'PARA', nbval=nbpara, vect=zr(lpar), nbret=n)
!
    call wkvect('&&OP0004.NOM.FONCTIONS', 'V V K24', nbfonc, lnomf)
    if (defonc) then
        call getvtx(' ', 'NOM_PARA_FONC', scal=zk24(lpro+6), nbret=l)
        mxva = 0
        do ifonc = 1, nbfonc
            call getvr8('DEFI_FONCTION', 'VALE', iocc=ifonc, nbval=0, nbret=nbval)
            mxva = max(mxva, -nbval)
        end do
        call wkvect('&&OP0004.VALEURS.LUES', 'V V R', mxva, jval)
        call wkvect('&&OP0004.POINTEURS.F', 'V V I', nbfonc, ladrf)
        do ifonc = 1, nbfonc
            zk24(lnomf+ifonc-1) = '&&OP0004.F'
            call codent(ifonc, 'G', zk24(lnomf+ifonc-1) (11:19))
            zk24(lnomf+ifonc-1) (20:24) = '.VALE'
            call getvr8('DEFI_FONCTION', 'VALE', iocc=ifonc, nbval=mxva, vect=zr(jval), &
                        nbret=nbval)
            call wkvect(zk24(lnomf+ifonc-1), 'V V R', nbval, lval)
            zi(ladrf+ifonc-1) = lval
            nbcoup = nbval/2
            do ival = 1, nbcoup
                zr(lval-1+ival) = zr(jval-1+2*ival-1)
                zr(lval-1+nbcoup+ival) = zr(jval-1+2*ival)
            end do
!
!           --- VERIFICATION QU'ON A BIEN CREER UNE FONCTION ---
!               ET REMISE DES ABSCISSES EN ORDRE CROISSANT
!           CE N'EST PAS LA PEINE SI LA CROISSANTE STRICTE A ETE IMPOSEE
            if (verif .ne. 'CROISSANT') then
                iret2 = 0
                call foverf(zr(lval), nbcoup, iret2)
                if (iret2 .eq. 0) then
                    typfon = 'FONCTION'
                    call uttrif(zr(lval), nbcoup, typfon)
                    valk = nomfon
                    call utmess('A', 'UTILITAI8_5', sk=valk)
                else if (iret2 .lt. 0) then
                    call ordon1(zr(lval), nbcoup)
                    vali(1) = ifonc
                    call utmess('A', 'UTILITAI8_6', si=vali(1))
                end if
            end if
!
            call getvtx('DEFI_FONCTION', 'INTERPOL', iocc=ifonc, nbval=2, vect=interp, &
                        nbret=l1)
            if (l1 .eq. 1) interp(2) = interp(1)
            zk24(lpro+6+2*ifonc-1) = interp(1)//interp(2)
            call getvtx('DEFI_FONCTION', 'PROL_GAUCHE', iocc=ifonc, &
                        scal=zk24(lpro+6+2*ifonc) (1:1), nbret=l)
            call getvtx('DEFI_FONCTION', 'PROL_DROITE', iocc=ifonc, &
                        scal=zk24(lpro+6+2*ifonc) (2:2), nbret=l)
        end do
    else
        call getvid(' ', 'FONCTION', nbval=nbfonc, vect=zk24(lnomf), nbret=n)
        call fovern(zk24(lnomf), nbfonc, zk24(lpro), iret)
    end if
!
!   --- CREATION ET REMPLISSAGE DE LA COLLECTION NOMFON.VALE ---
    call jecrec(nomfon//'.VALE', 'G V R', 'NU', 'CONTIG', 'VARIABLE', &
                nbfonc)
    call foston(nomfon//'.VALE', zk24(lnomf), nbfonc)
!
!   --- ON ORDONNE LA NAPPE SUIVANT LES PARAMETRES CROISSANTS ---
!       CE N'EST PAS LA PEINE SI LA CROISSANTE STRICTE A ETE IMPOSEE
    if (verif .ne. 'CROISSANT') then
        call ordonn(nomfon, 0)
    end if
!
!   --- CREATION D'UN TITRE ---
    call titre()
!
!   --- IMPRESSIONS ---
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, k8b)
!
    call jedema()
end subroutine
