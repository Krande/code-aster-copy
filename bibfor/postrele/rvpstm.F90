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
subroutine rvpstm(sdlieu, sdeval, sdmoye)
    implicit none
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: sdeval
    character(len=24) :: sdmoye, sdlieu
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
!     SDLIEU : SD DU LIEU TRAITE
!     SDEVAL : SD DE L' EVALUATION DE LA QUANTITE SUR CE LIEU
!
!  ARGUMENTS EN SORTIE
!  -------------------
!
!     SDMOYE : NOM DE LA SD CONSERVANT LES MOYENNES
!
!              XD V R, UN OC PAR OC DE .ABSC DU LIEU
!                      DIM(V) = 6*NB_CMP*NB_COUCHE*NB_SS_PT
!
!                      1 --> MOYENNE DE TYPE 1
!                      2 --> MOYENNE DE TYPE 2
!                      3 --> MINIMUM
!                      4 --> MAXIMUM
!                      5 --> MINIMUM LINEAIRE
!                      6 --> MAXIMUM LINEAIRE
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
    integer(kind=8) :: avale, apnbn, apadr, amoye, aabsc, atab, adr1, adr2
    integer(kind=8) :: deb, fin, lmoye, nbcp, nbco, nbsp, nboc, nbsgt
    integer(kind=8) :: l1, l2, l3, l5, l6, l7, ioc, ico, isgt, isp, k, i, n, inoe
    real(kind=8) :: m1, m2, ma, mi, s1, s2, t1, t2, s12, xl, t12, smil
    aster_logical :: deja
    character(len=1) :: bl
    character(len=4) :: docul, docu
    character(len=24) :: nvale, npnbn, npadr, nabsc, nnocp, ntab
!
!==================== CORPS DE LA ROUTINE =============================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icmp, iret, lll
!-----------------------------------------------------------------------
    call jemarq()
    bl = ' '
!
    ntab = '&&RVPSTM.VECT.INTER'
!
    nvale = sdeval//'.VALE'
    npnbn = sdeval//'.PNBN'
    nnocp = sdeval//'.NOCP'
    npadr = sdeval//'.PADR'
!
    nabsc = sdlieu(1:19)//'.ABSC'
!
    call jelira(sdlieu(1:19)//'.REFE', 'DOCU', cval=docul)
    call jelira(nvale, 'DOCU', cval=docu)
    call jelira(nabsc, 'NMAXOC', nboc)
    call jeexin(nnocp, iret)
    if (iret .eq. 0) then
        call utmess('F', 'POSTRELE_5')
    end if
    call jelira(nnocp, 'LONMAX', nbcp)
    call jecrec(sdmoye, 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                nboc)
!
    call jeveuo(sdeval//'.PNCO', 'L', avale)
!
    nbco = zi(avale)
!
    call jeveuo(sdeval//'.PNSP', 'L', avale)
!
    nbsp = zi(avale)
!
    call jeveuo(nvale, 'L', avale)
    call jeveuo(npadr, 'L', apadr)
    call jeveuo(npnbn, 'L', apnbn)
!
    l2 = nbsp*nbcp
    l1 = nbco*l2
    l3 = 2*l1
!
    if (l2 .gt. 6) then
        call utmess('F', 'POSTRELE_7')
    end if
    lmoye = 6*nbcp*nbco*nbsp
    fin = 0
!
    do ioc = 1, nboc, 1
!
        call jecroc(jexnum(sdmoye, ioc))
        call jeecra(jexnum(sdmoye, ioc), 'LONMAX', lmoye)
        call jeveuo(jexnum(sdmoye, ioc), 'E', amoye)
        call jelira(jexnum(nabsc, ioc), 'LONMAX', nbsgt)
        call jeveuo(jexnum(nabsc, ioc), 'L', aabsc)
        deja = .false.
!
        nbsgt = nbsgt-1
        if (nbsgt .eq. 0) then
            call utmess('F', 'POSTRELE_6')
        end if
        deb = fin+1
        fin = deb+nbsgt
!
        if ((docu .eq. 'CHLM') .or. (docul .ne. 'LSTN')) fin = fin-1
!
!     /* VECTEUR INTER */
!
        call wkvect(ntab, 'V V R', l3*(nbsgt+1), atab)
!
        if ((docul .eq. 'LSTN') .or. (docu .eq. 'CHNO')) then
!
            do isgt = 1, nbsgt+1, 1
!
                adr1 = zi(apadr+deb+isgt-2)
                n = zi(apnbn+deb+isgt-2)
!
                do ico = 1, nbco, 1
!
                    l5 = (ico-1)*n*l2
!
                    do k = 1, l2, 1
!
                        t1 = 0.0d0
!
                        lll = 0
                        do i = 1, n, 1
!
                            if (zr(avale-1+adr1+l5+(i-1)*l2+k-1) .eq. r8vide()) goto 230
                            lll = lll+1
                            t1 = t1+zr(avale-1+adr1+l5+(i-1)*l2+k-1)
!
230                         continue
                        end do
!
                        if (lll .eq. 0) then
                            t1 = r8vide()
                        else
                            t1 = t1/lll
                        end if
!
                        adr2 = (isgt-1)*l3+(ico-1)*l2+k
!
                        zr(atab+adr2-1) = t1
                        zr(atab+adr2+l1-1) = t1
!
                    end do
!
                end do
!
            end do
!
        else
!
            do isgt = 1, nbsgt, 1
!
                adr1 = zi(apadr+deb+isgt-2)
!
                do ico = 1, nbco, 1
!
                    l5 = (ico-1)*l2
!
                    do k = 1, l2, 1
!
                        adr2 = (isgt-1)*l3+l5+l1+k
!
                        zr(atab+adr2-1) = zr(avale+adr1+2*l5+k-2)
                        zr(atab+adr2+l1-1) = zr(avale+adr1+2*l5+l2+k-2)
!
                    end do
!
                end do
!
            end do
!
        end if
!
!     /* CONTRIBUTION ELEMENTAIRE */
!
        do icmp = 1, nbcp, 1
!
            xl = 0.d0
!
            do ico = 1, nbco, 1
!
                l5 = l2*(ico-1)
!
                do isp = 1, nbsp, 1
!
                    l6 = nbcp*(isp-1)
                    m1 = 0.0d0
                    m2 = 0.0d0
                    ma = -1.0d50
                    mi = 1.0d50
                    inoe = 0
!
                    do isgt = 1, nbsgt, 1
!
                        adr1 = l3*(isgt-1)+l5+l6+icmp
                        adr2 = adr1+l1
!
                        t1 = zr(atab-1+l1+adr1)
                        t2 = zr(atab-1+l1+adr2)
!
                        if (t1 .eq. r8vide()) then
                            inoe = inoe+1
                            goto 140
                        end if
                        if (t2 .eq. r8vide()) then
                            if (isgt .eq. nbsgt) inoe = inoe+1
                            goto 140
                        end if
!
                        s1 = zr(aabsc+isgt-1)-zr(aabsc)
                        s2 = zr(aabsc+isgt+1-1)-zr(aabsc)
                        s12 = s2-s1
                        xl = xl+s12
                        t12 = (t1+t2)/2.0d0
                        smil = (s1+s2)/2.0d0
                        m1 = m1+s12*(t1+t2)
                        m2 = m2+s12/3.0d0*(t1*s1+4.0d0*t12*smil+t2*s2)
                        ma = max(ma, t1, t2)
                        mi = min(mi, t1, t2)
!
140                     continue
                    end do
!
                    if (inoe .ne. 0) then
                        if (.not. deja) then
                            if (inoe .eq. 1) then
                                call utmess('A', 'POSTRELE_62', si=inoe)
                            else
                                call utmess('A', 'POSTRELE_63', si=inoe)
                            end if
                            deja = .true.
                        end if
                    end if
!
                    m1 = m1/xl
                    m2 = m2/(xl*xl)
!
                    m1 = 0.5d0*m1
                    m2 = 6.0d0*(m2-m1)
!
                    l7 = 6*(nbsp*(nbco*(icmp-1)+ico-1)+isp-1)
!
                    zr(amoye+l7+1-1) = m1
                    zr(amoye+l7+2-1) = m2
                    zr(amoye+l7+3-1) = mi
                    zr(amoye+l7+4-1) = ma
                    zr(amoye+l7+5-1) = m1-0.5d0*m2
                    zr(amoye+l7+6-1) = m1+0.5d0*m2
!
                end do
!
            end do
!
        end do
!
        call jedetr(ntab)
!
    end do
!
    call jedema()
end subroutine
