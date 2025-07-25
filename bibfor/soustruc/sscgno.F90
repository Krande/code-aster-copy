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

subroutine sscgno(ma, nbgnin)
    implicit none
!     BUT: TRAITER LE MOT CLEF CREA_GROUP_NO
!          DE L'OPERATEUR: DEFI_GROUP
!
!     IN:  MA    : NOM DU MAILLAGE
!          NBGNP : NOMBRE DE GROUP_NO A CREER
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/cgnoec.h"
#include "asterfort/cgnoes.h"
#include "asterfort/cgnofu.h"
#include "asterfort/cgnoin.h"
#include "asterfort/cgnoiv.h"
#include "asterfort/cgnoor.h"
#include "asterfort/cgnopl.h"
#include "asterfort/cgnoso.h"
#include "asterfort/cgnoxf.h"
#include "asterfort/cgrcbp.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ornofd.h"
#include "asterfort/ssgngm.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    real(kind=8) :: vecori(3)
!
    integer(kind=8) :: nalar, nbma
    character(len=8) :: alarm, typm, ndorig, ndextr
    character(len=8) :: ma, nono, k8b, kpos, nom1
    character(len=16) :: concep, cmd, option, motcle, typmcl, motfac
    character(len=24) :: grpnoe, cooval, lisno, mafour
    character(len=24) :: valk(2), nogno, nogno2
    character(len=80) :: card
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iagm1, iagm2, iagma, iagn2, iagno, ialii1
    integer(kind=8) :: ialii2, idlino, ifm, ign, ign1
    integer(kind=8) :: ign2, ii, iii, ili1, ili2, im1
    integer(kind=8) :: ind1, ind2, ino, iocc, ireste, iret, jjj
    integer(kind=8) :: jnoeu, jvale, kkk, maxcol, n, n1
    integer(kind=8) :: n2, n3, n4, n5, n6, n6a, n6b
    integer(kind=8) :: n7, n8, n9, nb, nbcol, nbgna2, nbgnaj
    integer(kind=8) :: nbgnin, nbgrmn, nbid, nbis, nbk8, nbline, nbno
    integer(kind=8) :: nbnot, nbocc, niv, ntrou, num, ier
    aster_logical :: l_write
    aster_logical :: l_parallel_mesh, lcolle
    character(len=24), pointer :: lik8(:) => null()
    character(len=8), pointer :: l_noeud(:) => null()
    integer(kind=8), pointer :: noeud2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
!     RECUPERATION DU NIVEAU D'IMPRESSION
!     -----------------------------------
    call infniv(ifm, niv)
!
    call getres(k8b, concep, cmd)
    lisno = '&&SSCGNO.LISTE_NOEUDS'
    grpnoe = ma//'.GROUPENO       '
    cooval = ma//'.COORDO    .VALE'
    call jeveuo(cooval, 'L', jvale)
    call jelira(grpnoe, 'NMAXOC', nbgrmn)
    nbis = nbgrmn
    nbk8 = nbgrmn
    call wkvect('&&SSCGNO.LII1', 'V V I', nbis, ialii1)
    call wkvect('&&SSCGNO.LII2', 'V V I', nbis, ialii2)
!
    motcle = 'GROUP_MA'
    typmcl = 'GROUP_MA'
!
    motfac = 'CREA_GROUP_NO'
    call getfac(motfac, nbocc)
    call getvtx(' ', 'ALARME', scal=alarm, nbret=nalar)
    nbgnaj = 0
!
    l_parallel_mesh = isParallelMesh(ma)
!
! ----------------------------------------------------------------------
!
    do iocc = 1, nbocc
!
        call getvtx(motfac, 'NOEUD', iocc=iocc, nbval=0, nbret=n2)
        call getvtx(motfac, 'INTERSEC', iocc=iocc, nbval=0, nbret=n3)
        call getvtx(motfac, 'UNION', iocc=iocc, nbval=0, nbret=n4)
        call getvtx(motfac, 'DIFFE', iocc=iocc, nbval=0, nbret=n5)
        call getvtx(motfac, 'GROUP_MA', iocc=iocc, nbval=0, nbret=n6)
        call getvtx(motfac, 'TOUT_GROUP_MA', iocc=iocc, nbval=0, nbret=n7)
        call getvtx(motfac, 'GROUP_NO', iocc=iocc, nbval=0, nbret=n8)
        call getvtx(motfac, 'OPTION', iocc=iocc, nbval=0, nbret=n9)
!
        if (l_parallel_mesh .and. n2 < 0) then
            call utmess('F', 'MODELISA7_86')
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "TOUT_GROUP_MA" :
!       --------------------------
        if (n7 .lt. 0) then
            call ssgngm(ma, iocc, nbgna2)
            nbgnaj = nbgnaj+nbgna2
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "GROUP_MA" :
!       ---------------------
        if (n6 .lt. 0 .and. n9 .eq. 0) then
            call ssgngm(ma, iocc, nbgna2)
            nbgnaj = nbgnaj+nbgna2
            goto 100
        end if
!
! ----------------------------------------------------------------------
!
        nogno = ''
        call getvtx(motfac, 'NOM', iocc=iocc, scal=nogno, nbret=n1)
        call jenonu(jexnom(grpnoe, nogno), iret)
        if (iret .gt. 0) then
            call utmess('F', 'SOUSTRUC_37', sk=nogno)
        end if
!
        n2 = -n2
        n3 = -n3
        n4 = -n4
        n5 = -n5
        n8 = -n8
        n9 = -n9
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "INTERSEC" :
!       ---------------------
        if (n3 .gt. 0) then
            AS_ALLOCATE(vk24=lik8, size=n3)
            call getvem(ma, 'GROUP_NO', motfac, 'INTERSEC', iocc, &
                        n3, lik8, nbid)
            n3 = nbid
            if (n3 .eq. 0) then
                AS_DEALLOCATE(vk24=lik8)
            end if
        end if
        if (n3 .gt. 0) then
            call jenonu(jexnom(grpnoe, lik8(1)), ign1)
            call jelira(jexnum(grpnoe, ign1), 'LONUTI', ili1)
            call jeveuo(jexnum(grpnoe, ign1), 'L', iagm1)
            if (ili1 .gt. nbis) then
                nbis = 2*ili1
                call jedetr('&&SSCGNO.LII1')
                call jedetr('&&SSCGNO.LII2')
                call wkvect('&&SSCGNO.LII1', 'V V I', nbis, ialii1)
                call wkvect('&&SSCGNO.LII2', 'V V I', nbis, ialii2)
            end if
            n = ili1
            do ii = 1, n
                zi(ialii1-1+ii) = zi(iagm1-1+ii)
            end do
!
            do ign = 2, n3
                call jenonu(jexnom(grpnoe, lik8(ign)), ign2)
                call jelira(jexnum(grpnoe, ign2), 'LONUTI', ili2)
                call jeveuo(jexnum(grpnoe, ign2), 'L', iagm2)
                call utlisi('INTER', zi(ialii1), n, zi(iagm2), ili2, &
                            zi(ialii2), nbis, ntrou)
                n = ntrou
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(ialii2-1+ii)
                end do
            end do
            AS_DEALLOCATE(vk24=lik8)
!
            if (n .eq. 0) then
                if (alarm .eq. 'OUI') then
                    call utmess('A', 'SOUSTRUC_38', sk=nogno)
                end if
                goto 100
            end if
            call jecroc(jexnom(grpnoe, nogno))
            call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, n))
            call jeecra(jexnom(grpnoe, nogno), 'LONUTI', n)
            call jeveuo(jexnom(grpnoe, nogno), 'E', iagma)
            do ii = 1, n
                zi(iagma-1+ii) = zi(ialii1-1+ii)
            end do
            nbgnaj = nbgnaj+1
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "UNION" :
!       ------------------
        if (n4 .gt. 0) then
            AS_ALLOCATE(vk24=lik8, size=n4)
            call getvem(ma, 'GROUP_NO', motfac, 'UNION', iocc, &
                        n4, lik8, nbid)
            n4 = nbid
            if (n4 .eq. 0) then
                AS_DEALLOCATE(vk24=lik8)
            end if
        end if
        if (n4 .gt. 0) then
            call jenonu(jexnom(grpnoe, lik8(1)), ign1)
            call jelira(jexnum(grpnoe, ign1), 'LONUTI', ili1)
            call jeveuo(jexnum(grpnoe, ign1), 'L', iagm1)
            if (ili1 .gt. nbis) then
                nbis = 2*ili1
                call jedetr('&&SSCGNO.LII1')
                call jedetr('&&SSCGNO.LII2')
                call wkvect('&&SSCGNO.LII1', 'V V I', nbis, ialii1)
                call wkvect('&&SSCGNO.LII2', 'V V I', nbis, ialii2)
            end if
            n = ili1
            do ii = 1, n
                zi(ialii1-1+ii) = zi(iagm1-1+ii)
            end do
!
            do ign = 2, n4
                call jenonu(jexnom(grpnoe, lik8(ign)), ign2)
                call jelira(jexnum(grpnoe, ign2), 'LONUTI', ili2)
                call jeveuo(jexnum(grpnoe, ign2), 'L', iagm2)
                call utlisi('UNION', zi(ialii1), n, zi(iagm2), ili2, &
                            zi(ialii2), nbis, ntrou)
!
                if (ntrou .lt. 0) then
                    nbis = -2*ntrou
                    call jedetr('&&SSCGNO.LII2')
                    call wkvect('&&SSCGNO.LII2', 'V V I', nbis, ialii2)
                    call utlisi('UNION', zi(ialii1), n, zi(iagm2), ili2, &
                                zi(ialii2), nbis, ntrou)
                    call jedetr('&&SSCGNO.LII1')
                    call wkvect('&&SSCGNO.LII1', 'V V I', nbis, ialii1)
                end if
                n = ntrou
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(ialii2-1+ii)
                end do
            end do
            AS_DEALLOCATE(vk24=lik8)
!
            if (n .eq. 0) then
                if (alarm .eq. 'OUI') then
                    call utmess('A', 'SOUSTRUC_38', sk=nogno)
                end if
            else
                call jecroc(jexnom(grpnoe, nogno))
                call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, n))
                call jeecra(jexnom(grpnoe, nogno), 'LONUTI', n)
                call jeveuo(jexnom(grpnoe, nogno), 'E', iagma)
                do ii = 1, n
                    zi(iagma-1+ii) = zi(ialii1-1+ii)
                end do
                nbgnaj = nbgnaj+1
            end if
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "DIFFE" :
!       ------------------
        if (n5 .gt. 0) then
            AS_ALLOCATE(vk24=lik8, size=n5)
            call getvem(ma, 'GROUP_NO', motfac, 'DIFFE', iocc, &
                        n5, lik8, nbid)
            n5 = nbid
            if (n5 .eq. 0) then
                AS_DEALLOCATE(vk24=lik8)
            end if
        end if
        if (n5 .gt. 0) then
            call jenonu(jexnom(grpnoe, lik8(1)), ign1)
            call jelira(jexnum(grpnoe, ign1), 'LONUTI', ili1)
            call jeveuo(jexnum(grpnoe, ign1), 'L', iagm1)
            if (ili1 .gt. nbis) then
                nbis = 2*ili1
                call jedetr('&&SSCGNO.LII1')
                call jedetr('&&SSCGNO.LII2')
                call wkvect('&&SSCGNO.LII1', 'V V I', nbis, ialii1)
                call wkvect('&&SSCGNO.LII2', 'V V I', nbis, ialii2)
            end if
            n = ili1
            do ii = 1, n
                zi(ialii1-1+ii) = zi(iagm1-1+ii)
            end do
!
            do ign = 2, n5
                call jenonu(jexnom(grpnoe, lik8(ign)), ign2)
                call jelira(jexnum(grpnoe, ign2), 'LONUTI', ili2)
                call jeveuo(jexnum(grpnoe, ign2), 'L', iagm2)
                call utlisi('DIFFE', zi(ialii1), n, zi(iagm2), ili2, &
                            zi(ialii2), nbis, ntrou)
                n = ntrou
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(ialii2-1+ii)
                end do
            end do
            AS_DEALLOCATE(vk24=lik8)
!
            if (n .eq. 0) then
                if (alarm .eq. 'OUI') then
                    call utmess('A', 'SOUSTRUC_38', sk=nogno)
                end if
            else
                call jecroc(jexnom(grpnoe, nogno))
                call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, n))
                call jeecra(jexnom(grpnoe, nogno), 'LONUTI', n)
                call jeveuo(jexnom(grpnoe, nogno), 'E', iagma)
                do ii = 1, n
                    zi(iagma-1+ii) = zi(ialii1-1+ii)
                end do
                nbgnaj = nbgnaj+1
            end if
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "OPTION" :
!       -------------------
        if (n9 .gt. 0) then
            call getvtx(motfac, 'OPTION', iocc=iocc, scal=option, nbret=n9)
!
!         -- TRAITEMENT DE L'OPTION "ENV_SPHERE" :
!         ----------------------------------------
            if (option(1:10) .eq. 'ENV_SPHERE') then
                call cgnoes(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION "ENV_CYLINDRE" :
!         ------------------------------------------
            else if (option(1:12) .eq. 'ENV_CYLINDRE') then
                call cgnoec(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION "PLAN" :
!         ----------------------------------
            else if (option(1:4) .eq. 'PLAN') then
                call cgnopl(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION "SEGM_DROI_ORDO" :
!         --------------------------------------------
            else if (option(1:14) .eq. 'SEGM_DROI_ORDO') then
                call cgnoso(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION "TUNNEL" :
!         ------------------------------------
            else if (option(1:6) .eq. 'TUNNEL') then
                call cgnofu(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION INCLUSION :
!         ------------------------------------
            else if (option .eq. 'INCLUSION') then
                call cgnoin(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION INTERVALLE_VALE :
!         --------------------------------------------
            else if (option .eq. 'INTERVALLE_VALE') then
                call cgnoiv(iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION "NOEUD_ORDO" :
!         ----------------------------------------
            else if (option(1:10) .eq. 'NOEUD_ORDO') then
                mafour = '&&SSCGNO.MALIGNE'
                call cgnoor(mafour, ma, motfac, iocc, 1, &
                            motcle, typmcl, ' ', nbma, ndorig, &
                            ndextr, typm, vecori)
                call ornofd(mafour, ma, nbma, lisno, ndorig, &
                            ndextr, 'V', vecori)
                call jedetr(mafour)
                call jelira(lisno, 'LONUTI', nbno)
!
!         -- TRAITEMENT DE L'OPTION FISS_XFEM :
!         ----------------------------------------
            else if (option .eq. 'FISS_XFEM') then
                call cgnoxf(motfac, iocc, ma, lisno, nbno)
!
!         -- TRAITEMENT DE L'OPTION RELA_CINE_BP :
!         ----------------------------------------
            else if (option .eq. 'RELA_CINE_BP') then
                l_write = .true.
                call cgrcbp(motfac, iocc, ma, l_write, nbgna2)
                nbgnaj = nbgnaj+nbgna2
                goto 100
!
            else
                call utmess('F', 'CALCULEL6_10', sk=option)
            end if
!
!         -- CREATION ET AFFECTATION DU GROUP_NO :
!         ----------------------------------------
            if (nbno .eq. 0) then
                if (alarm .eq. 'OUI') then
                    call utmess('A', 'SOUSTRUC_38', sk=nogno)
                end if
            else
                call jeveuo(lisno, 'L', idlino)
!
                call jecroc(jexnom(ma//'.GROUPENO', nogno))
                call jeecra(jexnom(ma//'.GROUPENO', nogno), 'LONMAX', max(1, nbno))
                call jeecra(jexnom(ma//'.GROUPENO', nogno), 'LONUTI', nbno)
                call jeveuo(jexnom(ma//'.GROUPENO', nogno), 'E', iagma)
!
                do ii = 1, nbno
                    zi(iagma-1+ii) = zi(idlino-1+ii)
                end do
                nbgnaj = nbgnaj+1
!
            end if
            call jedetr(lisno)
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "NOEUD" :
!       ------------------
        if (n2 .gt. 0) then
            if (l_parallel_mesh) then
                call utmess('F', 'MODELISA7_86')
            end if
            AS_ALLOCATE(vk8=l_noeud, size=n2)
            call getvem(ma, 'NOEUD', motfac, 'NOEUD', iocc, &
                        n2, l_noeud, nb)
            n2 = nb
            if (n2 .eq. 0) then
                AS_DEALLOCATE(vk8=l_noeud)
            end if
        end if
        if (n2 .gt. 0) then
            call wkvect('&&SSCGNO.NOEUD', 'V V I', n2, jnoeu)
            call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbnot)
            AS_ALLOCATE(vi=noeud2, size=nbnot)
!         --- ON VERIFIE QUE TOUS LES NOEUDS SONT DISTINCTS ---
            nbno = 0
            lcolle = .false.
            call jeexin(ma//'.NOMNOE', ier)
            if (ier .ne. 0) then
                lcolle = .true.
            end if
            do im1 = 1, n2
                nom1 = l_noeud(im1)
                num = char8_to_int(nom1, lcolle, ma, "NOEUD")
                noeud2(num) = noeud2(num)+1
                if (noeud2(num) .eq. 2) then
                    valk(1) = nom1
                    valk(2) = nogno
                    call utmess('A', 'SOUSTRUC_39', nk=2, valk=valk)
                    goto 20
                end if
                nbno = nbno+1
                zi(jnoeu+nbno-1) = num
20              continue
            end do
!
            call jecroc(jexnom(grpnoe, nogno))
            call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, nbno))
            call jeecra(jexnom(grpnoe, nogno), 'LONUTI', nbno)
            call jeveuo(jexnom(grpnoe, nogno), 'E', iagma)
            do ino = 0, nbno-1
                zi(iagma+ino) = zi(jnoeu+ino)
            end do
            nbgnaj = nbgnaj+1
            call jedetr('&&SSCGNO.NOEUD')
            AS_DEALLOCATE(vi=noeud2)
            AS_DEALLOCATE(vk8=l_noeud)
            goto 100
        end if
!
! ----------------------------------------------------------------------
! ----- MOT CLEF "GROUP_NO" :
!       ---------------------
        if (n8 .gt. 0) then
            call getvem(ma, 'GROUP_NO', motfac, 'GROUP_NO', iocc, &
                        1, nogno2, nbid)
            n8 = nbid
        end if
        if (n8 .gt. 0) then
            call getvem(ma, 'GROUP_NO', motfac, 'GROUP_NO', iocc, &
                        1, nogno2, nbid)
            call jenonu(jexnom(grpnoe, nogno2), ign2)
            call jelira(jexnum(grpnoe, ign2), 'LONUTI', ili2)
            call jeveuo(jexnum(grpnoe, ign2), 'L', iagn2)
!
            call getvtx(motfac, 'POSITION', iocc=iocc, nbval=0, nbret=n6b)
            ind1 = 0
            ind2 = 0
            if (n6b .eq. 0) then
                call getvis(motfac, 'NUME_INIT', iocc=iocc, scal=ind1, nbret=n6a)
                if (n6a .eq. 0) ind1 = 1
                call getvis(motfac, 'NUME_FIN', iocc=iocc, scal=ind2, nbret=n6a)
                if (n6a .eq. 0) ind2 = ili2
                if (ind2 .lt. ind1) then
                    call utmess('F', 'SOUSTRUC_33')
                end if
                if (ili2 .lt. ind2) then
                    call utmess('F', 'SOUSTRUC_34')
                end if
                n6a = ind2-ind1+1
            else
                n6a = 1
            end if
!
            call jecroc(jexnom(grpnoe, nogno))
            call jeecra(jexnom(grpnoe, nogno), 'LONMAX', max(1, n6a))
            call jeecra(jexnom(grpnoe, nogno), 'LONUTI', n6a)
            call jeveuo(jexnom(grpnoe, nogno), 'E', iagno)
            nbgnaj = nbgnaj+1
            if (n6b .ne. 0) goto 80
            n = ind2-ind1+1
            do ii = 1, n
                zi(iagno-1+ii) = zi(iagn2-2+ind1+ii)
            end do
            goto 100
80          continue
            call getvtx(motfac, 'POSITION', iocc=iocc, scal=kpos, nbret=n6b)
            if (kpos .eq. 'INIT') then
                zi(iagno) = zi(iagn2)
            else if (kpos .eq. 'FIN') then
                ii = ili2
                zi(iagno) = zi(iagn2+ii-1)
            else if (kpos .eq. 'MILIEU') then
                ii = (ili2+1)/2
                zi(iagno) = zi(iagn2+ii-1)
            end if
            goto 100
        end if
!
! ----------------------------------------------------------------------
!
100     continue
    end do
!
! ----------------------------------------------------------------------
! --- IMPRESSIONS NIVEAUX 1 ET 2 :
!     --------------------------
    if (niv .ge. 1 .and. nbgnaj .ne. 0) then
        write (ifm, '(/,/,A,I6,/,39(''=''))')&
     &    'NOMBRE  DE GROUPES DE NOEUDS CREES : ', nbgnaj
!
        if (nbocc .ge. 1) then
            write (ifm, '(/,15X,54(''-''),2(/,15X,A),/,15X,54(''-''))')&
     &      '!         NOM DU GROUPE         ! NBRE DE NOEUDS DU  !',&
     &      '!            NOEUDS             !      GROUPE_NO     !'
!
            do i = 1, nbgnaj
                ii = nbgnin+i
                call jenuno(jexnum(grpnoe, ii), nogno)
                call jelira(jexnum(grpnoe, ii), 'LONUTI', nbno)
                write (ifm, '(15X,A,2X,A24,5X,A,2X,I8,10X,A)') '!', &
                    nogno, '!', nbno, '!'
            end do
            write (ifm, '(15X,54(''-''),/)')
        end if
    end if
!
! ----------------------------------------------------------------------
! --- IMPRESSIONS NIVEAU 2 :
!     --------------------
    if (niv .eq. 2 .and. nbgnaj .ne. 0) then
        maxcol = 8
        lcolle = .false.
        call jeexin(ma//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if
        do i = 1, nbgnaj
            ii = nbgnin+i
            call jeveuo(jexnum(grpnoe, ii), 'L', iagno)
            call jenuno(jexnum(grpnoe, ii), nogno)
            call jelira(jexnum(grpnoe, ii), 'LONUTI', nbno)
            write (ifm, '(/,3A,/,27(''-''))') 'NOEUDS DU GROUPE ', &
                nogno, ' :'
            nbline = nbno/maxcol
            ireste = mod(nbno, maxcol)
            if (ireste .ne. 0) nbline = nbline+1
            nbcol = maxcol
            kkk = 0
            do jjj = 1, nbline
                if (ireste .ne. 0 .and. jjj .eq. nbline) nbcol = ireste
                do iii = 1, nbcol
                    kkk = kkk+1
                    nono = int_to_char8(zi(iagno-1+kkk), lcolle, ma, "NOEUD")
                    card((iii-1)*10+1:) = ' '//nono//' '
                end do
                write (ifm, '(A)') card(:10*nbcol)
            end do
        end do
        write (ifm, '(/,/)')
    end if
!
! ----------------------------------------------------------------------
! --- MENAGE
    call jedetr(lisno)
    call jedetr('&&SSCGNO.LII1')
    call jedetr('&&SSCGNO.LII2')
!
    call jedema()
!
end subroutine
