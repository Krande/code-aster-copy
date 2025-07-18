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

subroutine refode(nbcmb, angle, nomch, nuharm, tyharm, &
                  coef, basz, chpres)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/celver.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nbelem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbcmb, nuharm(*)
    character(len=*) :: nomch(*), basz, tyharm(*), chpres
    real(kind=8) :: angle, coef(*)
!     RECOMBINAISON DE FOURIER
!     -----------------------------------------------------------------
!     IN  NBCMB  : I   : NOMBRE DE CHAMPS A RECOMBINER
!     IN  ANGLE  : R8  : SECTION OU A LIEU LA RECOMBINAISON ( EN RD )
!     IN  NOMCH  : K8  : NOM DES CHAMPS A RECOMBINER
!     IN  NUHARM : I   : NUMERO DE L'HARMONIQUE
!     IN  TYHARM : K4  : TYPE DE L'HARMONIQUE (SYME OU ANTI)
!     IN  COEF   : R8  : COEF MULTIPLICATEUR ASSOCIE A L'HARMONIQUE
!     IN  CHPRES : K19 : NOM DU CHAMP RESULTAT
!     -----------------------------------------------------------------
    integer(kind=8) :: ibid, mode
    character(len=1) :: base
    character(len=4) :: docu
    character(len=5) :: refe, desc, vale
    character(len=8) :: noma, nomgd
    character(len=19) :: ch19, ligrel
    aster_logical :: lmeca, lther
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, ic, icoef, idecgr
    integer(kind=8) :: igrel, im, ino, ip, iret, ival
    integer(kind=8) :: jdesc, jrefe, jvale, k
    integer(kind=8) :: kdesc, krefe, kvale, lvale, nbdesc, nbec, nbelgr
    integer(kind=8) :: nbgr, nbnoeu, nbpt, nbrefe, nbscal, nbvale
    real(kind=8) :: ang
    integer(kind=8), pointer :: celd(:) => null()
    real(kind=8), pointer :: celv(:) => null()
    character(len=24), pointer :: celk(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    base = basz
    ch19 = nomch(1)
!
    call dismoi("DOCU", ch19, "CHAMP", repk=docu)

!
    if (docu .eq. 'CHNO') then
        desc = ' '
        refe = '.REFE'
        vale = '.VALE'
        nbdesc = 0
    else if (docu .eq. 'CHML') then
        desc = '.CELD'
        refe = '.CELK'
        vale = '.CELV'
        call jelira(ch19//desc, 'LONMAX', nbdesc)
        call jeveuo(ch19//desc, 'L', jdesc)
    else
        call utmess('F', 'UTILITAI_21')
    end if
!
    lmeca = .false.
    lther = .false.
    call dismoi('NOM_GD', ch19, 'CHAMP', repk=nomgd)
    if (nomgd .eq. 'DEPL_R' .or. nomgd .eq. 'SIEF_R' .or. nomgd .eq. 'EPSI_R') then
        lmeca = .true.
    else if (nomgd .eq. 'TEMP_R' .or. nomgd .eq. 'FLUX_R') then
        lther = .true.
    end if
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nbec)
!
!     --- CONSTRUCTION D'UN CHAMP RESULTAT SUR LE MODELE DE NOMCH(1)
!
    call jelira(ch19//vale, 'LONMAX', nbvale)
    call jelira(ch19//refe, 'LONMAX', nbrefe)
    call jeveuo(ch19//refe, 'L', jrefe)
!
    ch19 = chpres
    call jeexin(ch19//vale, iret)
    if (iret .eq. 0) then
        if (docu .eq. 'CHML') then
            call wkvect(ch19//desc, base//' V I', nbdesc, kdesc)
        end if
        call wkvect(ch19//vale, base//' V R', nbvale, kvale)
        call wkvect(ch19//refe, base//' V K24', nbrefe, krefe)
    else
        if (docu .eq. 'CHML') then
            call jeveuo(ch19//desc, 'E', kdesc)
        end if
        call jeveuo(ch19//vale, 'E', kvale)
        call jeveuo(ch19//refe, 'E', krefe)
    end if
!

    if (docu .eq. 'CHNO') then
        call jeecra(ch19//refe, 'DOCU', cval=docu)
    else if (docu .eq. 'CHML') then
        call jeecra(ch19//desc, 'DOCU', cval=docu)
        do i = 0, nbdesc-1
            zi(kdesc+i) = zi(jdesc+i)
        end do
    end if

!
!
    call wkvect('&&REFODE.VALE', 'V V R', nbvale, lvale)
!
    if (docu .eq. 'CHNO') then
        call jeveuo(zk24(jrefe+1) (1:19)//'.PRNO', 'L', vi=prno)
        do i = 0, nbrefe-1
            zk24(krefe+i) = zk24(jrefe+i)
        end do
        call dismoi('NOM_MAILLA', nomch(1), 'CHAMP', repk=noma)
        call jelira(noma//'.COORDO    .VALE', 'LONMAX', nbnoeu)
        nbnoeu = nbnoeu/3
!
!        --- BOUCLE SUR LES CHAMPS A RECOMBINER ---
!
        do im = 1, nbcmb
            ang = angle*dble(nuharm(im))
            ch19 = nomch(im)
            call jeveuo(ch19//vale, 'L', jvale)
!
            if (lmeca) then
!
                if (tyharm(im) (1:4) .eq. 'SYME') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*cos(ang)*zr(jvale+i+1)
                            zr(lvale+i+2) = zr(lvale+i+2)+coef(im)*cos(ang)*zr(jvale+i+2)
                            zr(lvale+i+3) = zr(lvale+i+3)-coef(im)*sin(ang)*zr(jvale+i+3)
                        end if
                    end do
!
                else if (tyharm(im) (1:4) .eq. 'ANTI') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*sin(ang)*zr(jvale+i+1)
                            zr(lvale+i+2) = zr(lvale+i+2)+coef(im)*sin(ang)*zr(jvale+i+2)
                            zr(lvale+i+3) = zr(lvale+i+3)+coef(im)*cos(ang)*zr(jvale+i+3)
                        end if
                    end do
!
                else if (tyharm(im) (1:4) .eq. 'TOUS') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*sin(ang)*zr(jvale+i+1)+co&
                                            &ef(im)*cos(ang)*zr(jvale+i+1)
                            zr(lvale+i+2) = zr(lvale+i+2)+coef(im)*sin(ang)*zr(jvale+i+2)+co&
                                            &ef(im)*cos(ang)*zr(jvale+i+2)
                            zr(lvale+i+3) = zr(lvale+i+3)+coef(im)*cos(ang)*zr(jvale+i+3)-co&
                                            &ef(im)*sin(ang)*zr(jvale+i+3)
                        end if
                    end do
!
                end if
!
            else if (lther) then
!
                if (tyharm(im) (1:4) .eq. 'SYME') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*cos(ang)*zr(jvale+i+1)
                        end if
                    end do
!
                else if (tyharm(im) (1:4) .eq. 'ANTI') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*sin(ang)*zr(jvale+i+1)
                        end if
                    end do
!
                else if (tyharm(im) (1:4) .eq. 'TOUS') then
!
                    do ino = 1, nbnoeu
                        i = prno(1+(ino-1)*(nbec+2))-2
                        if (i .ne. -2) then
                            zr(lvale+i+1) = zr(lvale+i+1)+coef(im)*sin(ang)*zr(jvale+i+1)+co&
                                            &ef(im)*cos(ang)*zr(jvale+i+1)
                        end if
                    end do
!
                end if
            end if
        end do
!
    else if (docu .eq. 'CHML') then
        do i = 0, nbrefe-1
            zk24(krefe+i) = zk24(jrefe+i)
        end do
!
        do im = 1, nbcmb
            i1 = -1
            ang = angle*dble(nuharm(im))
            ch19 = nomch(im)
!
!           -- ON VERIFIE QUE LE CHAM_ELEM N'EST PAS TROP DYNAMIQUE :
            call celver(ch19, 'NBVARI_CST', 'STOP', ibid)
            call celver(ch19, 'NBSPT_1', 'STOP', ibid)
!
            call jeveuo(ch19//'.CELD', 'L', vi=celd)
            call jeveuo(ch19//'.CELK', 'L', vk24=celk)
            call jeveuo(ch19//'.CELV', 'L', vr=celv)
            nbgr = celd(2)
            ligrel = celk(1) (1:19)
!
            do igrel = 1, nbgr
                mode = celd(celd(4+igrel)+2)
                if (mode .eq. 0) goto 210
                nbscal = digdel(mode)
                icoef = max(1, celd(4))
                if (icoef .ne. 1) then
                    call utmess('F', 'ALGELINE3_33')
                end if
                nbelgr = nbelem(ligrel, igrel)
                idecgr = celd(celd(4+igrel)+8)
!
                if (lmeca) then
!              --- ON EST EN AXIS, IL Y A 6 COMPOSANTES PAR POINT ---
                    nbpt = nbscal/6
!
                    if (tyharm(im) (1:4) .eq. 'SYME') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)-coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)-coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                            end do
                        end do
!
                    else if (tyharm(im) (1:4) .eq. 'ANTI') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                            end do
                        end do
!
                    else if (tyharm(im) (1:4) .eq. 'TOUS') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)-coef(im)*sin(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)-coef(im)*sin(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                            end do
                        end do
                    end if
!
                else if (lther) then
!              --- ON EST EN AXIS, IL Y A 3 COMPOSANTES PAR POINT ---
                    nbpt = nbscal/3
!
                    if (tyharm(im) (1:4) .eq. 'SYME') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)-coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                            end do
                        end do
!
                    else if (tyharm(im) (1:4) .eq. 'ANTI') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-&
                                               &1)*nbscal+ic)
                            end do
                        end do
!
                    else if (tyharm(im) (1:4) .eq. 'TOUS') then
!
                        do k = 1, nbelgr
                            ic = -1
                            do ip = 1, nbpt
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*sin(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)+coef(im)*cos(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                                i1 = i1+1
                                ic = ic+1
                                zr(lvale+i1) = zr(lvale+i1)+coef(im)*cos(ang)*celv(idecgr+(k-1&
                                               &)*nbscal+ic)-coef(im)*sin(ang)*celv(1-1+idecgr&
                                               &+(k-1)*nbscal+ic)
                            end do
                        end do
                    end if
                end if
210             continue
            end do
!
        end do
    end if
!
    do ival = 0, nbvale-1
        zr(kvale+ival) = zr(lvale+ival)
    end do
    ch19 = chpres
    call jedetr('&&REFODE.VALE')
!
    call jedema()
end subroutine
