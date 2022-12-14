! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine asexc2(motfac, nbocc, nbmode, momec, amort,&
                  corfre, noma, ndir, nomsup, nomspe,&
                  dirspe, echspe, nature, nbsupm, nsupp,&
                  knoeu, kvspe, kaspe, nopara, nordr)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/fointe.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer :: nbocc, nbmode, ndir(*), nature(3, *), nsupp(*), nordr(*)
    real(kind=8) :: amort(*), dirspe(3, *), echspe(3, *)
    character(len=8) :: nomsup(3, *), nomspe(3, *), noma
    character(len=*) :: motfac, kvspe, kaspe, knoeu, momec
    character(len=24) :: nopara(*)
    aster_logical :: corfre
!
!     COMMANDE : COMB_SISM_MODAL
!                TRAITEMENT DU MOT-CLE "EXCIT" POUR LE MONO-APPUI
!
!
! IN  : MOTFAC : MOT CLE FACTEUR
! IN  : NBOCC  : NOMBRE D'OCCURENCE DU MOT CLE FACTEUR
! IN  : NBMODE : NOMBRE DE MODES
! IN  : AMORT  : AMORTISSEMENTS MODAUX
! IN  : MOMEC  : MODES MECANIQUES
! IN  : CORFRE : CORRECTION FREQUENCE SI .TRUE.
! OUT : NDIR   : DIRECTION DU SEISME A ETUDIER
! OUT : VALSPE : VALEURS DU SPECTRE
! OUT : ASYSPE : VALEURS ASYMPTOTIQUES DU SPECTRE
!
!
!
!
!
    integer :: i, id, ier, igr, ii, iii, im, inat, ino, ioc
    integer :: iret, is, j, jaspe, jdgn, jkno
    integer :: jvspe, n1, nbsupm, ngr, nimpr
    integer :: nno
!
    real(kind=8) :: amor, coef, deuxpi, echel, epsi, freq, dirsp0(3)
    real(kind=8) :: echsp0(3), valpu(2), omega, omega2
    real(kind=8) :: resu, un, uns2pi, xnorm, zero, fcoup
!
    character(len=1) :: dir(3)
    character(len=4) :: knat
    character(len=8) :: k8b, spect, noeu, nomsp0(3), nompu(2), dircopy
    character(len=9) :: niveau
    character(len=24) :: obj1, obj2, valk(2), grnoeu
    integer :: ival
    character(len=24), pointer :: group_no(:) => null()
    character(len=8), pointer :: noeud(:) => null()
    real(kind=8) :: correc
!
!     ------------------------------------------------------------------
!
    data  nompu / 'AMOR' , 'FREQ'    /
    data   dir  / 'X' , 'Y' , 'Z' /
!
!     ------------------------------------------------------------------
!
    call jemarq()
    ier = 0
    epsi = 1.d-03
    zero = 0.d0
    un = 1.d0
    deuxpi = r8depi()
    uns2pi = un / deuxpi
    nsupp(1) = 0
    nsupp(2) = 0
    nsupp(3) = 0
    obj1 = noma//'.GROUPENO'
    obj2 = noma//'.NOMNOE'
!
!     --- LECTURE MOT-CLE FACTEUR IMPRESSION ---
!
    call getvtx('IMPRESSION', 'NIVEAU', iocc=1, scal=niveau, nbret=nimpr)
    if (nimpr .eq. 0) niveau='TOUT     '
!
    call getvr8(' ', 'FREQ_COUP', iocc=1, scal=fcoup, nbret=n1)
    if (n1 .eq. 0) then
        call rsadpa(momec, 'L', 1, nopara(2), nordr(nbmode),&
                    0, sjv=ival, istop=0)
        fcoup = uns2pi * sqrt(zr(ival))
    endif
!
!     --- NOMBRE DE SUPPORTS PAR DIRECTION ---
    do ioc = 1, nbocc
!
        echsp0(1) = un
        echsp0(2) = un
        echsp0(3) = un
        dirsp0(1) = un
        dirsp0(2) = un
        dirsp0(3) = un
        xnorm = un
!
!        --- UN SPECTRE SUIVANT UN AXE ---
        call getvr8(motfac, 'AXE', iocc=ioc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvr8(motfac, 'AXE', iocc=ioc, nbval=3, vect=dirsp0,&
                        nbret=n1)
            xnorm = zero
            do id = 1, 3
                xnorm = xnorm + dirsp0(id) * dirsp0(id)
            end do
            if (xnorm .lt. epsi) then
                ier = ier + 1
                call utmess('E', 'SEISME_4')
                goto 10
            endif
            xnorm = un / sqrt(xnorm)
            call getvid(motfac, 'SPEC_OSCI', iocc=ioc, scal=spect, nbret=n1)
            nomsp0(1) = spect
            nomsp0(2) = spect
            nomsp0(3) = spect
            call getvr8(motfac, 'ECHELLE', iocc=ioc, scal=echel, nbret=n1)
            if (n1 .ne. 0) then
                echsp0(1) = echel
                echsp0(2) = echel
                echsp0(3) = echel
            endif
!
!        --- UN SPECTRE DANS LES 3 DIRECTIONS ---
        else
            call getvr8(motfac, 'TRI_AXE', iocc=ioc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvr8(motfac, 'TRI_AXE', iocc=ioc, nbval=3, vect=dirsp0,&
                            nbret=n1)
                call getvid(motfac, 'SPEC_OSCI', iocc=ioc, scal=spect, nbret=n1)
                nomsp0(1) = spect
                nomsp0(2) = spect
                nomsp0(3) = spect
                call getvr8(motfac, 'ECHELLE', iocc=ioc, scal=echel, nbret=n1)
                if (n1 .ne. 0) then
                    echsp0(1) = echel
                    echsp0(2) = echel
                    echsp0(3) = echel
                endif
!
!        --- 3 SPECTRES DANS LES 3 DIRECTIONS ---
            else
!
                call getvid(motfac, 'SPEC_OSCI', iocc=ioc, nbval=3, vect=nomsp0,&
                            nbret=n1)
                call getvr8(motfac, 'ECHELLE', iocc=ioc, nbval=3, vect=echsp0,&
                            nbret=n1)
!
            endif
        endif
!
        call getvtx(motfac, 'NATURE', iocc=ioc, scal=knat, nbret=n1)
        if (knat .eq. 'ACCE') inat = 1
        if (knat .eq. 'VITE') inat = 2
        if (knat .eq. 'DEPL') inat = 3
!
        do id = 1, 3
            dirsp0(id) = xnorm * dirsp0(id)
            if (abs(dirsp0(id)) .gt. epsi) then
                ndir(id) = 1
!
                call getvem(noma, 'NOEUD', motfac, 'NOEUD', ioc,&
                            0, noeu, n1)
                if (n1 .ne. 0) then
                    nno = -n1
                    AS_ALLOCATE(vk8=noeud, size=nno)
                    call getvem(noma, 'NOEUD', motfac, 'NOEUD', ioc,&
                                nno, noeud, n1)
                    do ino = 1, nno
                        noeu = noeud(ino)
                        call jenonu(jexnom(obj2, noeu), iret)
                        if (iret .eq. 0) then
                            ier = ier + 1
                            valk(1) = noeu
                            valk(2) = noma
                            call utmess('E', 'SEISME_1', nk=2, valk=valk)
                            goto 20
                        endif
                        do is = 1, nsupp(id)
                            if (nomsup(id,is) .eq. noeu) then
                                ier = ier + 1
                                call utmess('E', 'SEISME_7', sk=noeu)
                                goto 20
                            endif
                        end do
                        nsupp(id) = nsupp(id) + 1
                        nomsup(id,nsupp(id)) = noeu
                        nomspe(id,nsupp(id)) = nomsp0(id)
                        dirspe(id,nsupp(id)) = dirsp0(id)
                        echspe(id,nsupp(id)) = echsp0(id)
                        nature(id,nsupp(id)) = inat
 20                     continue
                    end do
                    AS_DEALLOCATE(vk8=noeud)
!
                else
                    call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO', ioc,&
                                0, k8b, n1)
                    ngr = -n1
                    AS_ALLOCATE(vk24=group_no, size=ngr)
                    call getvem(noma, 'GROUP_NO', motfac, 'GROUP_NO', ioc,&
                                ngr, group_no, n1)
                    do igr = 1, ngr
                        grnoeu = group_no(igr)
                        call jeexin(jexnom(obj1, grnoeu), iret)
                        if (iret .eq. 0) then
                            ier = ier + 1
                            valk(1) = grnoeu
                            valk(2) = noma
                            call utmess('E', 'SEISME_2', nk=2, valk=valk)
                            goto 30
                        endif
                        call jelira(jexnom(obj1, grnoeu), 'LONUTI', nno)
                        call jeveuo(jexnom(obj1, grnoeu), 'L', jdgn)
                        do ino = 1, nno
                            call jenuno(jexnum(obj2, zi(jdgn+ino-1)), noeu)
                            do is = 1, nsupp(id)
                                if (nomsup(id,is) .eq. noeu) then
                                    ier = ier + 1
                                    call utmess('E', 'SEISME_7', sk=noeu)
                                    goto 32
                                endif
                            end do
                            nsupp(id) = nsupp(id) + 1
                            nomsup(id,nsupp(id)) = noeu
                            nomspe(id,nsupp(id)) = nomsp0(id)
                            dirspe(id,nsupp(id)) = dirsp0(id)
                            echspe(id,nsupp(id)) = echsp0(id)
                            nature(id,nsupp(id)) = inat
 32                         continue
                        end do
 30                     continue
                    end do
                    AS_DEALLOCATE(vk24=group_no)
                endif
            endif
        end do
!
 10     continue
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'SEISME_6')
    endif
!
!     --- NOM DES SUPPORTS PAR DIRECTION ---
    nbsupm = max(nsupp(1),nsupp(2),nsupp(3))
    call wkvect(knoeu, 'V V K8', 3*nbsupm, jkno)
    do is = 1, 3*nbsupm
        zk8(jkno+is-1) = '        '
    end do
    do id = 1, 3
        i = nbsupm*(id-1)
        do is = 1, nsupp(id)
            i = i + 1
            zk8(jkno+i-1) = nomsup(id,is)
        end do
    end do
!
!     --- INTERPOLATION DES SPECTRES ---
    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
        call utmess('I', 'SEISME_59')
    endif
    call wkvect(kvspe, 'V V R', 3*nbsupm*nbmode, jvspe)
    do im = 1, nbmode
        ii = 0
        amor = amort(im)
        call rsadpa(momec, 'L', 1, nopara(2), nordr(im),&
                    0, sjv=ival, istop=0)
        omega2 = zr(ival)
        omega = sqrt( omega2 )
        freq = uns2pi * omega
        valpu(1) = amor
        valpu(2) = freq
        if (corfre) then
            correc = sqrt( un - amor*amor )
        else
            correc =1.
        endif
        do id = 1, 3
            iii = 0
            if (ndir(id) .eq. 1) then
                do is = 1, nsupp(id)
                    call fointe('F ', nomspe(id, is), 2, nompu, valpu,&
                                resu, ier)
                    coef = dirspe(id,is)*echspe(id,is)
                    resu = resu * coef * correc
                    if (nature(id,is) .eq. 2) resu = resu * omega
                    if (nature(id,is) .eq. 3) resu = resu * omega2
                    j = id + 3*(im-1) + 3*nbmode*(is-1)
                    zr(jvspe+j-1) = resu
                    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
                        if (ii .eq. 0) then
                            ii = 1
                            iii = 1
                            dircopy = dir(id)
                            call utmess('I', 'SEISME_60', si=im, nk=2,&
                                        valk=[dircopy, nomsup(id, is)], nr=3,&
                                        valr=[freq, amor, resu])
                        else
                            if (iii .eq. 0) then
                                iii = 1
                                dircopy = dir(id)
                                call utmess('I', 'SEISME_73', nk=2,&
                                            valk=[dircopy, nomsup(id, is)], sr=resu)
                            else
                                call utmess('I', 'SEISME_61', sk=nomsup(id, is), sr=resu)
                            endif
                        endif
                    endif
                end do
            endif
        end do
    end do
!
!     --- VALEURS ASYMPTOTIQUES DES SPECTRES ---
    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
        call utmess('I', 'SEISME_62')
    endif
    call wkvect(kaspe, 'V V R', 3*nbsupm, jaspe)
    do id = 1, 3
        j = nbsupm*(id-1)
        if (ndir(id) .eq. 1) then
            iii = 0
            do is = 1, nsupp(id)
                coef = dirspe(id,is)*echspe(id,is)
                valpu(1) = amort(nbmode)
                valpu(2) = fcoup
                if (corfre) then
                    correc = sqrt( un - amor*amor )
                else
                    correc = 1.
                endif
                omega = deuxpi * fcoup
                call fointe('F ', nomspe(id, is), 2, nompu, valpu,&
                            resu, ier)
                if (nature(id,is) .eq. 2) resu = resu * omega
                if (nature(id,is) .eq. 3) resu = resu * omega * omega
                j = j + 1
                resu = resu * coef * correc
                zr(jaspe+j-1) = resu
                if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
                    if (iii .eq. 0) then
                        iii = 1
                        call utmess('I', 'SEISME_63', nk=2, valk=[dir(id), nomsup(id, is)],&
                                    sr=resu)
                    else
                        call utmess('I', 'SEISME_64', sk=nomsup(id, is), sr=resu)
                    endif
                endif
            end do
        endif
    end do
!
    call jedema()
end subroutine
