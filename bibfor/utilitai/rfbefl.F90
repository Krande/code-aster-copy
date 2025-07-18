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
subroutine rfbefl(base)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/foattr.h"
#include "asterfort/foimpr.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ordonn.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: base
!
!     OPERATEUR "RECU_FONCTION"  MOT CLE "BASE_ELAS_FLUI"
!     ------------------------------------------------------------------
    integer(kind=8) :: ifm, niv
    character(len=4) :: interp(2)
    character(len=8) :: k8b, basefl, ttordr, typflu
    character(len=16) :: nomcmd, typcon, parax, paray
    character(len=19) :: nomfon
    character(len=24) :: vite, numeo, numo, freq
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, id, imod, ind, ind1, ind2
    integer(kind=8) :: inumeo, j, lfon, lfreq, lnumo, lpro, lvar
    integer(kind=8) :: lvite, min, n1, n2, n3, n4, n5
    integer(kind=8) :: nbm, nbno, nbv, npv, nummod, iven
    integer(kind=8) :: ifsvr, ifsvi, lremf, nbzex, ivcn, nbconn, dec
    real(kind=8) :: pas, epsi
!-----------------------------------------------------------------------
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
    epsi = r8prem()
!
    call getres(nomfon, typcon, nomcmd)
    basefl = base
    nbv = 0
    interp(1) = 'LIN '
    interp(2) = 'LIN '
!
!     --- RECUPERATION DES ENTREES ---
!
    call getvtx(' ', 'PARA_X', scal=parax, nbret=n2)
    call getvtx(' ', 'PARA_Y', scal=paray, nbret=n3)
    call getvis(' ', 'NUME_MODE', scal=nummod, nbret=n4)
    call getvtx(' ', 'INTERPOL', nbval=2, vect=interp, nbret=n5)
    if (n5 .eq. 1) interp(2) = interp(1)
!
!     --- REMPLISSAGE DU .PROL ---
!
    ASSERT(lxlgut(nomfon) .le. 24)
    call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
    zk24(lpro) = 'FONCTION'
    zk24(lpro+1) = interp(1)//interp(2)
    zk24(lpro+2) = parax
    zk24(lpro+3) = paray
    zk24(lpro+4) = 'EE      '
    zk24(lpro+5) = nomfon
!
    if (parax(1:10) .eq. 'NB_CONNORS') then
!ASSERT(paray(1:9).eq.'VITE_CRIT')
!
        call jeveuo(basefl//'           .REMF', 'L', lremf)
        typflu = zk8(lremf+0)
        call jeveuo(typflu//'           .FSVR', 'L', ifsvr)
        call jeveuo(typflu//'           .FSVI', 'L', ifsvi)
        nbzex = zi(ifsvi+1)
        nbconn = zi(ifsvi+nbzex+2)
        pas = (zr(ifsvr+4)-zr(ifsvr+3))/(nbconn-1)
!
!       --- CAS D'UN NOMBRE DE CONNORS UNIQUE (PAS DE DISCRETISATION EN X)
!       --- FORCER UNE VALEUR MINIMUM DU PAS EN NOMBRE DE CONNORS
        if (abs(pas) .lt. epsi) then
            pas = epsi*1.d3
        end if
!
!------------- REMPLISSAGE DE L'ABSCISSE
        call wkvect(nomfon//'.VALE', 'G V R', 2*nbconn, lvar)
        lfon = lvar+nbconn
        do i = 1, nbconn
            zr(lvar+i-1) = zr(ifsvr+3)+(i-1)*pas
        end do
!
        call jeveuo(basefl//'.VCN', 'L', ivcn)
!------------- REMPLISSAGE DE L'ORDONNEE
        if (paray(1:7) .eq. 'INSTAB_') then
            call jeveuo(basefl//'.VEN', 'L', iven)
            dec = 0
!           ---- POUR LA METHODE TOUTES COMPOSANTS, DECALER DE NBMODES DANS .VEN
            if (paray(1:15) .eq. 'INSTAB_TOUT_CMP') then
                call jelira(basefl//'.VEN', 'LONMAX', dec)
                dec = dec/2
            end if
            do i = 1, nbconn
                if (abs(zr(ivcn+nbconn*(nummod-1)+i-1)) .gt. epsi) then
                    zr(lfon+i-1) = zr(iven+dec+nummod-1)/zr(ivcn+nbconn*(nummod-1)+i-1)
                else
                    call utmess('A', 'UTILITAI4_12', nr=1, valr=[zr(ifsvr+3)+(i-1)*pas])
                    zr(lfon+i-1) = 0.d0
                end if
            end do
        else if (paray(1:9) .eq. 'VITE_CRIT') then
            do i = 1, nbconn
                write (*, *) "i=", i, "vc=", zr(ivcn+nbconn*(nummod-1)+i-1)
                zr(lfon+i-1) = zr(ivcn+nbconn*(nummod-1)+i-1)
            end do
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(parax(1:8) .eq. 'VITE_FLU')
!
        call getvtx(' ', 'TOUT_ORDRE', scal=ttordr, nbret=n3)
!
!     --- RECUPERATION DES OJB ---
        vite = basefl//'           .VITE'
        call jelira(vite, 'LONUTI', npv)
        call jeveuo(vite, 'L', lvite)
        freq = basefl//'           .FREQ'
        call jeveuo(freq, 'L', lfreq)
        numo = basefl//'           .NUMO'
        call jelira(numo, 'LONUTI', nbm)
        call jeveuo(numo, 'L', lnumo)
!
!   --- VERIFICATION DE LA VALIDITE DES NUMEROS D'ORDRE DES VITESSES -
!
        if (ttordr .ne. 'OUI') then
            numeo = '&&RFBEFL.NUME_ORDRE'
            call getvis(' ', 'NUME_ORDRE', nbval=0, nbret=nbno)
            nbno = -nbno
            call wkvect(numeo, 'V V I', nbno, inumeo)
            call getvis(' ', 'NUME_ORDRE', nbval=nbno, vect=zi(inumeo), nbret=n1)
            min = zi(inumeo)
            do i = 1, nbno
                id = min-zi(inumeo+i-1)
                if (id .gt. 0) min = zi(inumeo+i-1)
            end do
            if (min .gt. npv) then
                call utmess('F', 'UTILITAI4_9')
            end if
        end if
!
!         --- DETERMINATION DU NUMERO D'ORDRE DU MODE VOULU ---
!
        do imod = 1, nbm
            id = nummod-zi(lnumo+imod-1)
            if (id .eq. 0) goto 30
        end do
        call utmess('F', 'UTILITAI4_10')
30      continue
!
!         --- CAS 1 : REMPLISSAGE POUR TOUS LES NUMEROS D'ORDRE ---
!
        if (ttordr .eq. 'OUI') then
            call wkvect(nomfon//'.VALE', 'G V R', 2*npv, lvar)
            lfon = lvar+npv
            do i = 1, npv
                zr(lvar+i-1) = zr(lvite+i-1)
                if (paray(1:4) .eq. 'FREQ') then
                    ind = 2*nbm*(i-1)+2*(imod-1)
                    zr(lfon+i-1) = zr(lfreq+ind)
                else
                    ind = 2*nbm*(i-1)+2*(imod-1)+1
                    zr(lfon+i-1) = zr(lfreq+ind)
                end if
            end do
        else
!
!         --- CAS 2 : REMPLISSAGE POUR UNE LISTE DE NUMEROS D'ORDRE ---
!
!-------------2.1 ON ORDONNE LA LISTE DES NUMEROS D'ORDRE
            if (nbno .gt. 1) then
                do i = 1, nbno
                    ind = i
                    min = zi(inumeo+i-1)
                    do j = i+1, nbno
                        id = min-zi(inumeo+j-1)
                        if (id .gt. 0) then
                            ind = j
                            min = zi(inumeo+j-1)
                        end if
                    end do
                    zi(inumeo+ind-1) = zi(inumeo+i-1)
                    zi(inumeo+i-1) = min
                end do
            end if
!
!-------------2.2 DETERMINATION DU NOMBRE DE NUMEROS D'ORDRE VALIDES
            if (nbno .gt. 1) then
                do i = 1, nbno
                    if (zi(inumeo+i-1) .gt. npv) goto 44
                    nbv = nbv+1
                end do
44              continue
            else
                nbv = 1
            end if
!
!-------------2.3 REMPLISSAGE
            call wkvect(nomfon//'.VALE', 'G V R', 2*nbv, lvar)
            lfon = lvar+nbv
            do i = 1, nbv
                ind1 = zi(inumeo+i-1)
                zr(lvar+i-1) = zr(lvite+ind1-1)
                if (paray(1:4) .eq. 'FREQ') then
                    ind2 = 2*nbm*(ind1-1)+2*(imod-1)
                    zr(lfon+i-1) = zr(lfreq+ind2)
                else
                    ind2 = 2*nbm*(ind1-1)+2*(imod-1)+1
                    zr(lfon+i-1) = zr(lfreq+ind2)
                end if
            end do
            call jedetr(numeo)
!
        end if
    end if
!
    call foattr(' ', 1, nomfon)
!
!     --- VERIFICATION QU'ON A BIEN CREER UNE FONCTION ---
!         ET REMISE DES ABSCISSES EN ORDRE CROISSANT
    call ordonn(nomfon, 0)
!
    call titre()
    if (niv .gt. 1) call foimpr(nomfon, niv, ifm, 0, k8b)
!
    call jedema()
end subroutine
