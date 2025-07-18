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
subroutine calcsp(casint, nomu, table, freq, masg, &
                  nbm, nbmr, imod1, nuor, ivite)
    implicit none
!     CALCUL POUR CHAQUE VITESSES DES INTERSPECTRES DE REPONSES
!-----------------------------------------------------------------------
! IN  : CASINT: BOOLEEN CARACTERISANT L'OPTION DE CALCUL
!       CASINT= TRUE   => ON CALCULE TOUS LES INTERSPECTRES
!       CASINT= FALSE  => ON CALCULE LES AUTOSPECTRES UNIQUEMENT
! IN  : NOMU  : NOM UTILISATEUR DU CONCEPT TABL_INTSP CORRESPONDANT
!               AUX INTERSPECTRES DE REPONSES : A PRODUIRE
! IN  : TABLE : NOM UTILISATEUR DU CONCEPT TABL_INTSP CORRESPONDANT
!               AUX INTERSPECTRES D'EXCITATIONS : DONNEE DU CALCUL
! IN  : FREQ  : TABLEAU DES FREQUENCES ET AMORTISSEMENTS
! IN  : MASG  : TABLEAU DES MASSES GENERALISEES
! IN  : NBM   : NOMBRE DE MODES DE LA BASE DE CONCEPT MELASFLU
! IN  : NBMR  : NOMBRE DE MODES PRIS EN COMPTE
! IN  : IMOD1 : INDICE DU PREMIER MODE PRIS EN COMPTE DANS LA BASE DE
!               CONCEPT MELASFLU
! IN  : NUOR  : LISTE DES NUMEROS D'ORDRE DES MODES PRIS EN COMPTE
! IN  : IVITE : NUMERO VITESSE DU FLUIDE
!     ----------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: casint
    character(len=8) :: nomu, table
    integer(kind=8) :: nbm, nbmr, imod1, nuor(*), ivite
    real(kind=8) :: freq(*), masg(*)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ideb, ifonc, ihi, ihi1, ihr
    integer(kind=8) :: ihr1, il, im, im1, im2, imb, imodf
    integer(kind=8) :: ip, iv, lvale, nbpf
!
    real(kind=8) :: fr, fri, hhi, hhr, hii1, hii2, hir1
    real(kind=8) :: hir2, pi
!-----------------------------------------------------------------------
    integer(kind=8) :: ival(3), vali(2)
    integer(kind=8) :: lnumi, lnumj, lfreq, i1, nbabs
    integer(kind=8) :: lrnumi, lrnumj, lrfreq, mxval, mrxval, ipf
    real(kind=8) :: mgi, ksi, hdenom
    character(len=24) :: chnumi, chnumj, chfreq, chvale
    character(len=24) :: crnumi, crnumj, crfreq, crvale
!
!-----------------------------------------------------------------------
    call jemarq()
!
    pi = r8pi()
    imodf = imod1+nbmr-1
    if (imodf .gt. nbm) then
        call utmess('F', 'MODELISA2_76', ni=2, vali=[nbm, nbmr])
    end if
!
    chnumi = table//'.NUMI'
    chnumj = table//'.NUMJ'
    chfreq = table//'.DISC'
    chvale = table//'.VALE'
    call jeveuo(chnumi, 'L', lnumi)
    call jeveuo(chnumj, 'L', lnumj)
    call jeveuo(chfreq, 'L', lfreq)
    call jelira(chnumi, 'LONMAX', mxval)
    call jelira(chfreq, 'LONMAX', nbpf)
!
!
    crnumi = nomu//'.NUMI'
    crnumj = nomu//'.NUMJ'
    crfreq = nomu//'.DISC'
    crvale = nomu//'.VALE'
    call wkvect(crfreq, 'G V R', nbpf, lrfreq)
    do ip = 1, nbpf
        zr(lrfreq+ip-1) = zr(lfreq+ip-1)
    end do
!
    mrxval = 0
    do im2 = 1, nbmr
        ideb = im2
        if (casint) ideb = 1
        do im1 = ideb, im2
            mrxval = mrxval+1
        end do
    end do
!
    call wkvect(crnumi, 'G V I', mrxval, lrnumi)
    call wkvect(crnumj, 'G V I', mrxval, lrnumj)
!
    call jecrec(crvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mrxval)
!
! --- CREATION DE VECTEURS DE TRAVAIL ---
!
    call wkvect('&&CALCSP.TEMP.HR  ', 'V V R8', nbmr*nbpf, ihr)
    call wkvect('&&CALCSP.TEMP.HI  ', 'V V R8', nbmr*nbpf, ihi)
!
    iv = ivite
    do im = imod1, imodf
        fri = freq(2*nbm*(iv-1)+2*(im-1)+1)
        if (fri .lt. 0.d0) then
            vali(1) = nuor(im)
            vali(2) = iv
            call utmess('F', 'MODELISA2_90', ni=2, vali=vali)
            goto 20
        end if
    end do
!
    do im = imod1, imodf
!
        mgi = masg(im)*4.d0*pi*pi
        fri = freq(2*nbm*(iv-1)+2*(im-1)+1)
        ksi = freq(2*nbm*(iv-1)+2*(im-1)+2)
        if (ksi .le. 0.d0) ksi = 1.d-06
!
        imb = im-imod1+1
!
        do ip = 1, nbpf
            fr = zr(lfreq+ip-1)
            ihr1 = ihr+nbpf*(imb-1)+ip-1
            ihi1 = ihi+nbpf*(imb-1)+ip-1
            zr(ihr1) = (mgi*(fri*fri-fr*fr))
            zr(ihi1) = (mgi*ksi*fr*fri*2.d0)
!
        end do
    end do
!
    ipf = 1
    do im2 = 1, nbmr
!
        ival(2) = nuor(im2)
!
        ideb = im2
        if (casint) ideb = 1
!
        do im1 = ideb, im2
!
            ival(3) = nuor(im1)
!
            do i1 = 1, mxval
                if ((zi(lnumi-1+i1) .eq. ival(2)) .and. (zi(lnumj-1+i1) .eq. ival(3))) then
                    call jeveuo(jexnum(chvale, i1), 'L', ifonc)
                end if
            end do
!
            call jecroc(jexnum(crvale, ipf))
            zi(lrnumi-1+ipf) = ival(2)
            zi(lrnumj-1+ipf) = ival(3)
            if (ival(2) .eq. ival(3)) then
                nbabs = nbpf
            else
                nbabs = 2*nbpf
            end if
            call jeecra(jexnum(crvale, ipf), 'LONMAX', nbabs)
            call jeecra(jexnum(crvale, ipf), 'LONUTI', nbabs)
            call jeveuo(jexnum(crvale, ipf), 'E', lvale)
            ipf = ipf+1
!
            do il = 1, nbpf
                hir1 = zr(ihr+nbpf*(im1-1)+il-1)
                hii1 = zr(ihi+nbpf*(im1-1)+il-1)
                hir2 = zr(ihr+nbpf*(im2-1)+il-1)
                hii2 = zr(ihi+nbpf*(im2-1)+il-1)
                hdenom = (hir1*hir1+hii1*hii1)*(hir2*hir2+hii2*hii2)
                hhr = (hir1*hir2+hii1*hii2)/hdenom
                hhi = (hir2*hii1-hir1*hii2)/hdenom
!
                if (ival(2) .eq. ival(3)) then
                    zr(lvale+il-1) = hhr*zr(ifonc+il-1)
                else
                    zr(lvale+2*(il-1)) = hhr*zr(ifonc+2*(il-1))-hhi*zr(ifonc+2*(il-1)+1)
                    zr(lvale+2*(il-1)+1) = hhr*zr(ifonc+2*(il-1)+1)+ &
                                           hhi*zr(ifonc+2*(il-1))
                end if
            end do
!
        end do
    end do
20  continue
!
!
! --- MENAGE
!
    call jedetr('&&CALCSP.TEMP.HR  ')
    call jedetr('&&CALCSP.TEMP.HI  ')
!
    call jedema()
end subroutine
