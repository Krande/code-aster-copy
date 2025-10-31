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
subroutine crcore()
    implicit none
!
!     COMMANDE:  CREA_RESU /CONV_RESU
!     CREE UNE STRUCTURE DE DONNEE DE TYPE
!           "DYNA_TRANS"   "EVOL_CHAR"
!
! --- ------------------------------------------------------------------
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/fointe.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/refdaj.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsmxno.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rssepa.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
!
    integer(kind=8) :: ibid, ier, icompt, iret, numini, numfin
    integer(kind=8) :: n1, nis, nbinst, nbval, nume, j
    integer(kind=8) :: iad, jinst, jchin, jchout, iddl, icmp, aprno
    integer(kind=8) :: nbv(1), jrefe
    integer(kind=8) :: jcpt, nbr, ivmx, k, iocc, nboini, inoe, ncmp, nbr0
    integer(kind=8) :: tnum(1)
    integer(kind=8) :: nbordr1, nbordr2, numei, neq, nbnoeu, gd, nec, iec
    integer(kind=8) :: ngn, nbno, ino, in, nbd, jchi1, ldgn, jdist
    integer(kind=8) :: nfx, nfy, nfz, jncmp, tabec(11), ncmpmx, iad2, ncmp0, j2
    integer(kind=8) :: nfrx, nfry, nfrz
!
    real(kind=8) :: rbid, tps, prec, coefr
    real(kind=8) :: dist, dire(3), coorre(3), vs, delay, dt, tp2
    real(kind=8) :: valpar(7)
    complex(kind=8) :: cbid
    aster_logical :: lddlex
!
    character(len=1) :: typmat
    character(len=4) :: typabs
    character(len=8) :: k8b, resu, criter, resui, matr, nomcmp
    character(len=8) :: modele, materi, carele, blan8, noma
    character(len=8) :: nompar(7)
    character(len=14) :: numedd
    character(len=16) :: type, oper, fonx, fony, fonz
    character(len=16) :: fonrx, fonry, fonrz
    character(len=19) :: nomch, listr8, list_load, resu19, profch
    character(len=19) :: chamno, chamn2, chamn1
    character(len=24) :: linst, nsymb, nsymb0, typres, lcpt, o1, o2
    character(len=24) :: matric(3), nprno
    character(len=24) :: nomgr, magrno, ldist
    character(len=8) :: nomgd
    character(len=24), pointer :: refn(:) => null()
    real(kind=8), pointer :: val(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: vicmp(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    data linst, listr8, lcpt, ldist/'&&CRCORE_LINST', '&&CRCORE_LISR8',&
     &     '&&CPT_CRCORE', '&&CRCORE_LDIST'/
! --- ------------------------------------------------------------------
    call jemarq()
!
    blan8 = ' '
    list_load = ' '
    nboini = 10
    modele = ' '
    carele = ' '
    materi = ' '
!
    call getres(resu, type, oper)
    resu19 = resu
    call getvtx(' ', 'TYPE_RESU', scal=typres, nbret=n1)
    iocc = 1
    call getvid('CONV_RESU', 'RESU_INIT', iocc=iocc, scal=resui, nbret=n1)
!
    call rscrsd('G', resu, typres, nboini)
!
    call jelira(resu//'           .ORDR', 'LONUTI', nbordr1)
    call getvtx('CONV_RESU', 'NOM_CHAM_INIT', iocc=iocc, scal=nsymb0, nbret=ibid)
    call getvr8('CONV_RESU', 'COEF', iocc=iocc, scal=coefr, nbret=ibid)
!
    numini = -1
    icompt = -1
    profch = ' '
!
!        MOT CLE INST PRESENT :
    nis = 0
    nbinst = 0
    call getvr8('CONV_RESU', 'INST', iocc=iocc, nbval=0, nbret=nis)
    if (nis .ne. 0) then
        typabs = 'INST'
        nbinst = -nis
    end if
!
    if (nis .ne. 0) then
        call wkvect(lcpt, 'V V I', nbinst, jcpt)
        call wkvect(linst, 'V V R', nbinst, jinst)
        call getvr8('CONV_RESU', typabs, iocc=iocc, nbval=nbinst, vect=zr(jinst), &
                    nbret=n1)
        call getvr8('CONV_RESU', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
        call getvtx('CONV_RESU', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
        call rsorac(resu, 'LONUTI', 0, rbid, k8b, &
                    cbid, rbid, k8b, nbv, 1, &
                    ibid)
!
        ivmx = rsmxno(resu)
        do k = 1, nbinst
            if (nbv(1) .gt. 0) then
                call rsorac(resu, typabs, ibid, zr(jinst+k-1), k8b, &
                            cbid, prec, criter, tnum, 1, &
                            nbr)
                nume = tnum(1)
            else
                nbr = 0
            end if
            if (nbr .lt. 0) then
                call utmess('F', 'ALGORITH2_48')
            else if (nbr .eq. 0) then
                zi(jcpt+k-1) = ivmx+1
                ivmx = ivmx+1
            else
                zi(jcpt+k-1) = nume
            end if
        end do
    else
!        MOT CLE LIST_INST PRESENT :
        n1 = 0
        call getvid('CONV_RESU', 'LIST_INST', iocc=iocc, scal=listr8, nbret=n1)
        if (n1 .ne. 0) then
            typabs = 'INST'
        end if
!
        call getvr8('CONV_RESU', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
        call getvtx('CONV_RESU', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
        call jelira(listr8//'.VALE', 'LONMAX', nbval)
!
        nbinst = nbval
        numini = 1
        numfin = nbinst
        nbinst = min(nbinst, nbval)
!
        call wkvect(linst, 'V V R', nbinst, jinst)
        call jeveuo(listr8//'.VALE', 'L', vr=val)
        call rsorac(resu, 'LONUTI', 0, rbid, k8b, &
                    cbid, rbid, k8b, nbv, 1, &
                    ibid)
        call wkvect(lcpt, 'V V I', nbinst, jcpt)
        ivmx = rsmxno(resu)
        j = 0
        do k = 1, nbval
            if (k .lt. numini) goto 40
            if (k .gt. numfin) goto 40
            j = j+1
            zr(jinst-1+j) = val(k)
            if (nbv(1) .gt. 0) then
                call rsorac(resu, typabs, ibid, val(k), k8b, &
                            cbid, prec, criter, tnum, 1, &
                            nbr)
                nume = tnum(1)
            else
                nbr = 0
            end if
            if (nbr .lt. 0) then
                call utmess('F', 'SEISME_78')
            else if (nbr .eq. 0) then
                zi(jcpt+j-1) = ivmx+1
                ivmx = ivmx+1
            else
                zi(jcpt+j-1) = nume
            end if
40          continue
        end do
    end if
    dt = zr(jinst+1)-zr(jinst)
    numedd = ' '
    call getvid('CONV_RESU', 'MATR_RIGI', iocc=iocc, scal=matr, nbret=n1)
    if (n1 .eq. 1) then
        call dismoi('NOM_NUME_DDL', matr, 'MATR_ASSE', repk=numedd)
    else
        call getvid('CONV_RESU', 'NUME_DDL', iocc=iocc, scal=numedd, nbret=n1)
        matr = ' '
    end if
    call dismoi('NOM_MODELE', numedd, 'NUME_DDL', repk=modele)
    call dismoi('NUME_EQUA', numedd, 'NUME_DDL', repk=profch)
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoeu)
    call dismoi('NUM_GD_SI', numedd, 'NUME_DDL', repi=gd)
    nprno = numedd//'.NUME.PRNO'
    call jenonu(jexnom(nprno(1:19)//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(nprno, ibid), 'L', aprno)
    nec = nbec(gd)
    typmat = 'R'
    call jeveuo(numedd//'.NUME.REFN', 'L', vk24=refn)
    nomgd = refn(2) (1:8)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jncmp)
    ASSERT(nec .le. 11)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iad2)
!    AS_ALLOCATE(vi=vicmp, size=ncmpmx*nbnoeu)
    call getvid('CONV_RESU', 'FONC_DX', iocc=iocc, scal=fonx, nbret=nfx)
    call getvid('CONV_RESU', 'FONC_DY', iocc=iocc, scal=fony, nbret=nfy)
    call getvid('CONV_RESU', 'FONC_DZ', iocc=iocc, scal=fonz, nbret=nfz)
    call getvid('CONV_RESU', 'FONC_DRX', iocc=iocc, scal=fonrx, nbret=nfrx)
    call getvid('CONV_RESU', 'FONC_DRY', iocc=iocc, scal=fonry, nbret=nfry)
    call getvid('CONV_RESU', 'FONC_DRZ', iocc=iocc, scal=fonrz, nbret=nfrz)
    if (nfx .ne. 0) then
        AS_ALLOCATE(vi=vicmp, size=ncmpmx*nbnoeu)
    end if
    if (typres(1:10) .eq. 'DYNA_TRANS') then
        nsymb = 'DEPL'
        if (nfx .ne. 0) then
            nsymb = nsymb0
        end if
    else
        nsymb = 'FORC_NODA'
    end if
    chamn2 = '&&CRCORE.CHAM_NO'
    call vtcreb(chamn2, 'V', 'R', nume_ddlz=numedd)
    call rsagsd(resu, nbinst)
    lddlex = .false.
    call getvtx('CONV_RESU', 'DDL_EXCLUS', iocc=iocc, scal=nomcmp, nbret=ibid)
    if (ibid .ne. 0) then
        lddlex = .true.
    end if
    if (lddlex) then
        if (nomcmp .eq. 'DX') icmp = 1
        if (nomcmp .eq. 'DY') icmp = 2
        if (nomcmp .eq. 'DZ') icmp = 3
        if (nomcmp .eq. 'DRX') icmp = 4
        if (nomcmp .eq. 'DRY') icmp = 5
        if (nomcmp .eq. 'DRZ') icmp = 6
    end if
!
    call getvem(noma, 'GROUP_NO', 'CONV_RESU', 'GROUP_NO_INTERF', 1, &
                1, nomgr, ngn)
    if (ngn .ne. 0) then
        magrno = noma//'.GROUPENO'
        call jeveuo(jexnom(magrno, nomgr), 'L', ldgn)
        call jelira(jexnom(magrno, nomgr), 'LONUTI', nbno)
        call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
        call getvr8('CONV_RESU', 'DIRECTION', iocc=iocc, nbval=3, vect=dire, &
                    nbret=nbd)
        dist = 0.d0
        do in = 1, 3
            dist = dist+dire(in)**2.0
        end do
        dist = dist**0.5
        do in = 1, 3
            dire(in) = dire(in)/dist
        end do
        call getvr8('CONV_RESU', 'COOR_REFE', iocc=iocc, nbval=3, vect=coorre, &
                    nbret=nbd)
        call getvr8('CONV_RESU', 'VITE_ONDE', iocc=iocc, scal=vs, nbret=nbd)
        call wkvect(ldist, 'V V R', nbno, jdist)
        do ino = 1, nbno
            inoe = zi(ldgn+ino-1)
            dist = 0.d0
            do in = 1, 3
                dist = dist+(vale(3*(inoe-1)+in)-coorre(in))*dire(in)
            end do
            if (dist .lt. 0.d0) then
                call utmess('A', 'SEISME_78')
            end if
            zr(jdist+ino-1) = dist
        end do
    end if
!
    do j = 1, nbinst
        if (j .ge. 2) call jemarq()
        call jerecu('V')
        icompt = zi(jcpt+j-1)
        tps = zr(jinst+j-1)
        call rsexch(' ', resu, nsymb, icompt, nomch, &
                    iret)
        if (iret .eq. 0) then
            call rsadpa(resu, 'L', 1, typabs, icompt, &
                        0, sjv=iad, styp=k8b)
        else if (iret .eq. 110) then
            call rsagsd(resu, 0)
            call rsexch(' ', resu, nsymb, icompt, nomch, iret)
        else if (iret .eq. 100) then
            call vtcreb(nomch, 'G', 'R', nume_ddlz=numedd)
        end if
        call jeveuo(nomch//'.VALE', 'E', jchout)
        call rsorac(resui, typabs, ibid, tps, k8b, &
                    cbid, prec, criter, tnum, 1, &
                    nbr)
        numei = tnum(1)
        nbr0 = nbr
        call rsexch(' ', resui, nsymb0, numei, chamno, iret)
        call vtcopy(chamno, chamn2, ier)
        if (ier .ne. 0) then
            call utmess("A", "FIELD0_3", sk=nsymb0)
        end if
!
        call jeveuo(chamn2//'.VALE', 'L', jchin)
        if (lddlex) then
            do inoe = 1, nbnoeu
                iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)
                ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
                if (iddl .ne. 0) then
                    if (icmp .gt. ncmp) goto 50
                    zr(jchin-1+iddl+icmp-1) = 0.d0
                end if
50              continue
            end do
        end if
        if (ngn .ne. 0) then
            do ino = 1, nbno
                inoe = zi(ldgn+ino-1)
                iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)
                ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
                delay = zr(jdist+ino-1)/vs
                tp2 = tps-dt*int(delay/dt)
                call rsorac(resui, typabs, ibid, tp2, k8b, &
                            cbid, prec, criter, tnum, 1, &
                            nbr)
                if (iddl .ne. 0) then
                    if (nbr .ne. 0) then
                        numei = tnum(1)
                        call rsexch(' ', resui, nsymb0, numei, chamn1, &
                                    iret)
                        call jeveuo(chamn1//'.VALE', 'L', jchi1)
                        do icmp = 1, ncmp
                            zr(jchin-1+iddl+icmp-1) = zr(jchi1-1+iddl+icmp-1)
                        end do
                    else
                        do icmp = 1, ncmp
                            zr(jchin-1+iddl+icmp-1) = 0.d0
                        end do
                    end if
                end if
            end do
        end if
        if (nfx .ne. 0) then
            nbr = nbr0
            if (nbr .ne. 0) then
                call jeveuo(chamno//'.VALE', 'L', jchi1)
            end if
            nompar(1) = 'DX'
            nompar(2) = 'DY'
            nompar(3) = 'DZ'
            nompar(4) = 'DRX'
            nompar(5) = 'DRY'
            nompar(6) = 'DRZ'
            nompar(7) = 'INST'
            valpar(7) = tps
            do inoe = 1, nbnoeu
                iddl = zi(aprno+(nec+2)*(inoe-1)+1-1)
                ncmp = zi(aprno+(nec+2)*(inoe-1)+2-1)
                if (j .gt. 1) goto 25
                do iec = 1, nec
                    tabec(iec) = zi(aprno-1+(inoe-1)*(nec+2)+2+iec)
                end do
                ncmp0 = 0
                do icmp = 1, ncmpmx
                    if (exisdg(tabec, icmp)) then
                        do j2 = 1, ncmp0
                            if (vicmp(ncmpmx*(inoe-1)+j2) .eq. icmp) goto 20
                        end do
                        ncmp0 = ncmp0+1
                        vicmp(ncmpmx*(inoe-1)+ncmp0) = icmp
                    end if
20                  continue
                end do
25              continue
                if (iddl .ne. 0) then
                    if (nbr .ne. 0) then
                        valpar(1) = zr(jchi1-1+iddl)
                        if (nfy .ne. 0) then
                            valpar(2) = zr(jchi1-1+iddl+1)
                        else
                            valpar(2) = 0.d0
                        end if
                        if (nfz .ne. 0) then
                            valpar(3) = zr(jchi1-1+iddl+2)
                        else
                            valpar(3) = 0.d0
                        end if
                        if (nfrx .ne. 0) then
                            valpar(4) = zr(jchi1-1+iddl+3)
                        else
                            valpar(4) = 0.d0
                        end if
                        if (nfry .ne. 0) then
                            valpar(5) = zr(jchi1-1+iddl+4)
                        else
                            valpar(5) = 0.d0
                        end if
                        if (nfrz .ne. 0) then
                            valpar(6) = zr(jchi1-1+iddl+5)
                        else
                            valpar(6) = 0.d0
                        end if
                        do icmp = 1, ncmp
                            zr(jchin-1+iddl+icmp-1) = zr(jchi1-1+iddl+icmp-1)
                            nomcmp = zk8(iad2-1+vicmp(ncmpmx*(inoe-1)+icmp))
!  write(6,*) 'icmp nomcmp vicmp =',icmp, nomcmp, vicmp(ncmpmx*(inoe-1)+icmp)
                            if (vicmp(ncmpmx*(inoe-1)+icmp) .eq. 1) then
                                call fointe('F ', fonx, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                            if (nfy .ne. 0 .and. vicmp(ncmpmx*(inoe-1)+icmp) .eq. 2) then
                                call fointe('F ', fony, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                            if (nfz .ne. 0 .and. vicmp(ncmpmx*(inoe-1)+icmp) .eq. 3) then
                                call fointe('F ', fonz, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                            if (nfrx .ne. 0 .and. vicmp(ncmpmx*(inoe-1)+icmp) .eq. 4) then
                                call fointe('F ', fonrx, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                            if (nfry .ne. 0 .and. vicmp(ncmpmx*(inoe-1)+icmp) .eq. 5) then
                                call fointe('F ', fonry, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                            if (nfrz .ne. 0 .and. vicmp(ncmpmx*(inoe-1)+icmp) .eq. 6) then
                                call fointe('F ', fonrz, 7, nompar, valpar, &
                                            zr(jchin-1+iddl+icmp-1), iret)
                            end if
                        end do
                    else
                        do icmp = 1, ncmp
                            zr(jchin-1+iddl+icmp-1) = 0.d0
                        end do
                    end if
                end if
            end do
        end if
!
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, coefr, zr(jchin), b_incx, zr(jchout), &
                   b_incy)
!
        o1 = chamno//'.DESC'
        o2 = nomch//'.DESC'
        call jedupo(o1, 'G', o2, .false._1)
!
        o1 = chamno//'.REFE'
        o2 = nomch//'.REFE'
        call jedupo(o1, 'G', o2, .false._1)
!
        call jeveuo(nomch//'.REFE', 'E', jrefe)
        zk24(jrefe+1) = profch
!
        call rsnoch(resu, nsymb, icompt)
        call rsadpa(resu, 'E', 1, typabs, icompt, &
                    0, sjv=iad, styp=k8b)
        zr(iad) = tps
        call rssepa(resu, icompt, modele, materi, carele, &
                    list_load)
        if (j .ge. 2) call jedema()
!
    end do
    call jedetr(linst)
    call jedetr(lcpt)
!
!     REMPLISSAGE DE .REFD POUR DYNA_*:
    call jelira(resu//'           .ORDR', 'LONUTI', nbordr2)
    if (nbordr2 .gt. nbordr1) then
        if (typres(1:10) .eq. 'DYNA_TRANS') then
            matric(1) = matr
            matric(2) = ' '
            matric(3) = ' '
            call refdaj('F', resu19, (nbordr2-nbordr1), numedd, 'DYNAMIQUE', &
                        matric, ier)
        end if
    end if
!    AS_DEALLOCATE(vi=vicmp)
    if (nfx .ne. 0) then
        AS_DEALLOCATE(vi=vicmp)
    end if
!
    call jedema()
end subroutine
