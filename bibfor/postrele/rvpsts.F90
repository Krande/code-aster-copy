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
subroutine rvpsts(iocc, sdlieu, sdeval, sdmoye)
    implicit none
!
!
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rvchve.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: sdeval
    character(len=24) :: sdmoye, sdlieu
    integer(kind=8) :: iocc
!
!**********************************************************************
!
!  OPERATION REALISEE
!  ------------------
!
!     OPERATION MOYENNE DU POST-TRAITEMENT POUR UN LIEU
!
!  ARGUMENTS EN ENTREE
!  -------------------
!
!     SDLIEU : SD DU LIEU TRAITEE
!     SDEVAL : SD DE L' EVALUATION DE LA QUANTITE SUR CE LIEU
!
!  ARGUMENTS EN SORTIE
!  -------------------
!
!     SDMOYE : NOM DE LA SD CONSERVANT LA SOMME ::= RECORD
!
!      .VALE : XD V R, UN OC PAR OC DE .ABSC DU LIEU
!                      DIM(V) = NB_CMP_SOMME*NB_COUCHE*NB_SS_PT
!      .NOCP : S V K8  NOM DES CMP
!
!**********************************************************************
!
!  FONCTIONS EXTERNES
!  ------------------
!
!
!  -----------------------------------------
!
!
!  ---------------------------------
!
!  VARIABLES LOCALES
!  -----------------
!
    character(len=24) :: nabsc, ntab
    character(len=24) :: nsocp, nsova
    character(len=4) :: docul, docu
!
    integer(kind=8) :: avale, amoye, aabsc, atab, adr1, adr2
    integer(kind=8) :: deb, fin, lmoye, nbcp, nbco, nbsp, nboc, nbsgt, nres, nmom
    integer(kind=8) :: l1, l2, l3, l5, ioc, i, j, k, l, n, nbpt
!
    real(kind=8) :: t1, t(3), x, y, z, xpi, ypi, zpi, rx, ry, rz, mx, my, mz
    real(kind=8) :: zero
!
!==================== CORPS DE LA ROUTINE =============================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: lll
    integer(kind=8), pointer :: pnbn(:) => null()
    integer(kind=8), pointer :: padr(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    zero = 0.0d0
    xpi = 0.d0
    ypi = 0.d0
    zpi = 0.d0
!
    ntab = '&&RVPSTM.VECT.INTER'
    nabsc = sdlieu(1:19)//'.ABSC'
    nsova = sdmoye(1:19)//'.VALE'
    nsocp = sdmoye(1:19)//'.NOCP'
!
    call getvtx('ACTION', 'RESULTANTE', iocc=iocc, nbval=0, nbret=nres)
    call getvtx('ACTION', 'MOMENT', iocc=iocc, nbval=0, nbret=nmom)
    nres = -nres
    nmom = -nmom
!
!   -- si on utilise RESULTANTE ou MOMENT, on verifie que REPERE est correct :
    if (nres .gt. 0 .or. nmom .gt. 0) then
        call rvchve(iocc, xpi, ypi, zpi)
    end if
!
    call jelira(sdlieu(1:19)//'.REFE', 'DOCU', cval=docul)
    call jelira(sdeval//'.VALE', 'DOCU', cval=docu)
    call jelira(nabsc, 'NMAXOC', nboc)
    call jelira(sdeval//'.NOCP', 'LONMAX', nbcp)
    call jecrec(nsova, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nboc)
!
    call jeveuo(sdeval//'.PNCO', 'L', i)
    nbco = zi(i)
    call jeveuo(sdeval//'.PNSP', 'L', i)
    nbsp = zi(i)
!
    call jeveuo(sdeval//'.VALE', 'L', avale)
    call jeveuo(sdeval//'.PADR', 'L', vi=padr)
    call jeveuo(sdeval//'.PNBN', 'L', vi=pnbn)
!
    l2 = nbsp*nbcp
    l1 = nbco*l2
    l3 = 2*l1
!
    if (nmom .eq. 0) then
        call wkvect(nsocp, 'V V K8', nbcp, deb)
        call jeveuo(sdeval//'.NOCP', 'L', fin)
        do ioc = 1, nbcp, 1
            zk8(deb+ioc-1) = zk8(fin+ioc-1)
        end do
        lmoye = nbcp*nbco*nbsp
    else
        call wkvect(nsocp, 'V V K8', 6, deb)
        zk8(deb-1+1) = 'RESULT_X'
        zk8(deb-1+2) = 'RESULT_Y'
        zk8(deb-1+3) = 'RESULT_Z'
        zk8(deb-1+4) = 'MOMENT_X'
        zk8(deb-1+5) = 'MOMENT_Y'
        zk8(deb-1+6) = 'MOMENT_Z'
        lmoye = 6*nbco*nbsp
    end if
!
    nbsgt = 0
    do ioc = 1, nboc, 1
        call jelira(jexnum(nabsc, ioc), 'LONMAX', fin)
        nbsgt = max(nbsgt, fin)
    end do
    call wkvect(ntab, 'V V R', l3*(nbsgt+1), atab)
!
    fin = 0
!
    do ioc = 1, nboc, 1
!
        call jecroc(jexnum(nsova, ioc))
        call jeecra(jexnum(nsova, ioc), 'LONMAX', lmoye)
        call jeveuo(jexnum(nsova, ioc), 'E', amoye)
        call jelira(jexnum(nabsc, ioc), 'LONMAX', nbpt)
        call jeveuo(jexnum(nabsc, ioc), 'L', aabsc)
!
        nbsgt = nbpt-1
        deb = fin+1
        fin = deb+nbsgt
!
        if ((docu .eq. 'CHLM') .or. (docul .ne. 'LSTN')) then
!
            fin = fin-1
!
        end if
!
!     /* VECTEUR INTER IE : PASSAGE A UN SS_CHAM_NO */
!
        if ((docul .eq. 'LSTN') .or. (docu .eq. 'CHNO')) then
!
            do i = 1, nbpt, 1
!
                adr1 = padr(1+deb+i-2)
                n = pnbn(1+deb+i-2)
!
                do j = 1, nbco, 1
!
                    l5 = (j-1)*n*l2
!
                    do k = 1, l2, 1
!
                        t1 = 0.0d0
!
                        lll = 0
                        do l = 1, n, 1
!
                            if (zr(avale-1+adr1+l5+(l-1)*l2+k-1) .eq. r8vide()) goto 230
                            lll = lll+1
                            t1 = t1+zr(avale-1+adr1+l5+(l-1)*l2+k-1)
!
230                         continue
                        end do
!
                        if (lll .eq. 0) then
                            t1 = 0.d0
                        else
                            t1 = t1/lll
                        end if
!
                        adr2 = (i-1)*l1+(j-1)*l2+k-1
!
                        zr(atab+adr2) = t1
!
                    end do
!
                end do
!
            end do
!
        else
!
            adr1 = padr(deb)
            do i = 1, nbco, 1
                l5 = (i-1)*l2
                do j = 1, l2, 1
                    if (zr(avale+adr1+2*l5+j-2) .eq. r8vide()) then
                        zr(atab+l5+j-1) = 0.d0
                    else
                        zr(atab+l5+j-1) = zr(avale+adr1+2*l5+j-2)
                    end if
                end do
            end do
!
            do i = 1, nbsgt-1, 1
!
                adr1 = padr(1+deb+i-2)
!
                do j = 1, nbco, 1
!
                    l5 = (j-1)*l2+i*l1
                    adr2 = avale+adr1-1+(j-1)*l3+l2
!
                    do k = 1, l2, 1
!
                        if (zr(adr2+k-1) .eq. r8vide() .and. zr(adr2+l2*(nbco*2-1)+k-1) &
                            .eq. r8vide()) then
                            zr(atab+l5+k-1) = 0.d0
                        else if (zr(adr2+k-1) .eq. r8vide()) then
                            zr(atab+l5+k-1) = zr(adr2+l2*(nbco*2-1)+k-1)
                        elseif (zr(adr2+l2*(nbco*2-1)+k-1) .eq. r8vide() &
                                ) then
                            zr(atab+l5+k-1) = zr(adr2+k-1)
                        else
                            zr(atab+l5+k-1) = 0.5d0*(zr(adr2+k-1)+zr(adr2+l2*(nbco*2-1)+k-1))
                        end if
!
                    end do
!
                end do
!
            end do
!
            adr1 = avale+padr(1+deb+nbsgt-2)-1
            adr2 = atab+nbsgt*l1
            do j = 1, nbco, 1
                l5 = (j-1)*l2
                do k = 1, l2, 1
                    if (zr(adr1+(i-1)*l3+k-1) .eq. r8vide()) then
                        zr(adr2+l5+k-1) = 0.d0
                    else
                        zr(adr2+l5+k-1) = zr(adr1+(i-1)*l3+k-1)
                    end if
                end do
            end do
!
        end if
!
!
!     /* CALCUL DES SOMMES SUR LE SS_CHAM_NO */
!
        if (nmom .eq. 0) then
            do i = 1, l1, 1
                t1 = zero
                do j = 1, nbpt, 1
                    t1 = t1+zr(atab+(j-1)*l1+i-1)
                end do
                zr(amoye+i-1) = t1
            end do
        else
            do i = 1, l1, 1
                zr(amoye+i-1) = zero
            end do
            call getvr8('ACTION', 'POINT', iocc=iocc, nbval=0, nbret=n)
            n = -n
            call getvr8('ACTION', 'POINT', iocc=iocc, nbval=n, vect=t, &
                        nbret=i)
!
!       -- point P : (x,y,z)
            x = t(1)
            y = t(2)
            if (n .eq. 2) then
                z = zero
            else
                z = t(3)
            end if
            call jeveuo(jexnum(sdlieu(1:19)//'.COOR', ioc), 'L', k)
            do i = 1, nbco*nbsp, 1
                adr1 = (i-1)*nbcp
                adr2 = (i-1)*6
                do j = 1, nbpt, 1
                    l = (j-1)*3+k-1
                    l5 = (j-1)*l1
!
!           -- calcul du vecteur PM :
                    xpi = zr(l+1)-x
                    ypi = zr(l+2)-y
                    zpi = zr(l+3)-z
!
!           -- changement de repere pour le vecteur PM :
                    call rvchve(iocc, xpi, ypi, zpi)
!
!
                    if (nmom .eq. 3) then
                        rx = zr(atab+l5+adr1+1-1)
                        ry = zr(atab+l5+adr1+2-1)
                        rz = zr(atab+l5+adr1+3-1)
                        mx = zr(atab+l5+adr1+4-1)
                        my = zr(atab+l5+adr1+5-1)
                        mz = zr(atab+l5+adr1+6-1)
                    else
                        rx = zr(atab+l5+adr1+1-1)
                        ry = zr(atab+l5+adr1+2-1)
                        rz = zero
                        mx = zr(atab+l5+adr1+3-1)
                        my = zr(atab+l5+adr1+4-1)
                        mz = zero
                    end if
!
                    mx = ypi*rz-zpi*ry+mx
                    my = zpi*rx-xpi*rz+my
                    mz = xpi*ry-ypi*rx+mz
!
                    zr(amoye+adr2+1-1) = zr(amoye+adr2+1-1)+rx
                    zr(amoye+adr2+2-1) = zr(amoye+adr2+2-1)+ry
                    zr(amoye+adr2+3-1) = zr(amoye+adr2+3-1)+rz
                    zr(amoye+adr2+4-1) = zr(amoye+adr2+4-1)+mx
                    zr(amoye+adr2+5-1) = zr(amoye+adr2+5-1)+my
                    zr(amoye+adr2+6-1) = zr(amoye+adr2+6-1)+mz
                end do
            end do
        end if
!
    end do
!
    call jedetr(ntab)
!
    call jedema()
end subroutine
