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

subroutine crkucv()
    implicit none
!
!     COMMANDE:  CREA_RESU /KUCV
!     CREE UNE STRUCTURE DE DONNEE DE TYPE
!           "DYNA_TRANS"   "EVOL_CHAR"
!
! --- ------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
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
!
    integer(kind=8) :: ibid, ier, icompt, iret, numini, numfin
    integer(kind=8) :: n1, nis, nbinst, nbval, nume, j
    integer(kind=8) :: iad, jinst, jchin, jchout
    integer(kind=8) :: nbv(1), jrefe
    integer(kind=8) :: jcpt, nbr, ivmx, k, iocc, nboini
    integer(kind=8) :: tnum(1)
    integer(kind=8) :: nbordr1, nbordr2, numei
    integer(kind=8) :: lmat, lma2, nr, neq
!
    real(kind=8) :: rbid, tps, prec
    complex(kind=8) :: cbid
!
    character(len=1) :: typmat
    character(len=4) :: typabs
    character(len=8) :: k8b, resu, criter, matr, resui
    character(len=8) :: modele, materi, carele, blan8
    character(len=14) :: numem
    character(len=16) :: type, oper
    character(len=19) :: nomch, listr8, list_load, resu19, profch
    character(len=19) :: chamno, chamn2
    character(len=19) :: mrigi, mamor
    character(len=24) :: linst, nsymb, typres, lcpt, o1, o2
    character(len=24) :: matric(3)
    real(kind=8), pointer :: val(:) => null()
!
    data linst, listr8, lcpt/'&&CRKUCV_LINST', '&&CRKUCV_LISR8',&
     &     '&&CPT_CRKUCV'/
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
    call getvid('KUCV', 'RESU_INIT', iocc=iocc, scal=resui, nbret=n1)
!
    call rscrsd('G', resu, typres, nboini)
!
    call jelira(resu//'           .ORDR', 'LONUTI', nbordr1)
!
    numini = -1
    icompt = -1
    profch = ' '
!
!        MOT CLE INST PRESENT :
    nis = 0
    nbinst = 0
    call getvr8('KUCV', 'INST', iocc=iocc, nbval=0, nbret=nis)
    if (nis .ne. 0) then
        typabs = 'INST'
        nbinst = -nis
    end if

    if (nis .ne. 0) then
        call wkvect(lcpt, 'V V I', nbinst, jcpt)
        call wkvect(linst, 'V V R', nbinst, jinst)
        call getvr8('KUCV', typabs, iocc=iocc, nbval=nbinst, vect=zr(jinst), &
                    nbret=n1)
        call getvr8('KUCV', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
        call getvtx('KUCV', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
        call rsorac(resu, 'LONUTI', 0, rbid, k8b, cbid, rbid, k8b, nbv, 1, ibid)
!
        ivmx = rsmxno(resu)
        do k = 1, nbinst
            if (nbv(1) .gt. 0) then
                call rsorac(resu, typabs, ibid, zr(jinst+k-1), k8b, &
                            cbid, prec, criter, tnum, 1, nbr)
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
        call getvid('KUCV', 'LIST_INST', iocc=iocc, scal=listr8, nbret=n1)
        if (n1 .ne. 0) then
            typabs = 'INST'
        end if
!
        call getvr8('KUCV', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
        call getvtx('KUCV', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
        call jelira(listr8//'.VALE', 'LONMAX', nbval)
!
        nbinst = nbval
        numini = 1
        numfin = nbinst

        nbinst = min(nbinst, nbval)
!
        call wkvect(linst, 'V V R', nbinst, jinst)
        call jeveuo(listr8//'.VALE', 'L', vr=val)
        call rsorac(resu, 'LONUTI', 0, rbid, k8b, cbid, rbid, k8b, nbv, &
                    1, ibid)
        call wkvect(lcpt, 'V V I', nbinst, jcpt)
        ivmx = rsmxno(resu)
        j = 0
        do k = 1, nbval
            if (k .lt. numini) goto 40
            if (k .gt. numfin) goto 40
            j = j+1
            zr(jinst-1+j) = val(k)
            if (nbv(1) .gt. 0) then
                call rsorac(resu, typabs, ibid, val(k), k8b, cbid, prec, criter, tnum, &
                            1, nbr)
                nume = tnum(1)
            else
                nbr = 0
            end if
            if (nbr .lt. 0) then
                call utmess('F', 'ALGORITH2_48')
            else if (nbr .eq. 0) then
                zi(jcpt+j-1) = ivmx+1
                ivmx = ivmx+1
            else
                zi(jcpt+j-1) = nume
            end if
40          continue
        end do
    end if
!
!
!     0. MATRICE :
!     -------------
    call getvid('KUCV', 'MATR_AMOR', iocc=iocc, scal=mamor, nbret=n1)

    call mtdscr(mamor)
    call jeveuo(mamor(1:19)//'.&INT', 'E', lmat)
    call dismoi('NOM_NUME_DDL', mamor, 'MATR_ASSE', repk=numem)
    call dismoi('NOM_MODELE', mamor, 'MATR_ASSE', repk=modele)
    call dismoi('NUME_EQUA', mamor, 'MATR_ASSE', repk=profch)
    call dismoi('NB_EQUA', numem, 'NUME_DDL', repi=neq)
    call getvid('KUCV', 'MATR_RIGI', iocc=iocc, scal=mrigi, nbret=nr)
    if (nr .ne. 0) then
        call mtdscr(mrigi)
        call jeveuo(mrigi(1:19)//'.&INT', 'E', lma2)
    end if
    typmat = 'R'
    if (typres(1:10) .eq. 'DYNA_TRANS') then
        nsymb = 'DEPL'
    else
        nsymb = 'FORC_NODA'
    end if
    chamn2 = '&&CRKUCV.CHAM_NO'
    call vtcreb(chamn2, 'V', 'R', nume_ddlz=numem)
    call rsagsd(resu, nbinst)
!
    do j = 1, nbinst
        if (j .ge. 2) call jemarq()
        call jerecu('V')
        icompt = zi(jcpt+j-1)
        tps = zr(jinst+j-1)
        call rsexch(' ', resu, nsymb, icompt, nomch, iret)
        if (iret .eq. 0) then
            call rsadpa(resu, 'L', 1, typabs, icompt, &
                        0, sjv=iad, styp=k8b)
        else if (iret .eq. 110) then
            call rsagsd(resu, 0)
            call rsexch(' ', resu, nsymb, icompt, nomch, iret)
        else if (iret .eq. 100) then
            call vtcreb(nomch, 'G', 'R', nume_ddlz=numem)
        end if
        call jeveuo(nomch//'.VALE', 'E', jchout)
        call rsorac(resui, typabs, ibid, tps, k8b, cbid, prec, criter, tnum, 1, nbr)
        numei = tnum(1)
        call rsexch(' ', resui, 'VITE', numei, chamno, iret)
        call vtcopy(chamno, chamn2, ier)
        if (ier .ne. 0) then
            call utmess("A", "FIELD0_3", sk='VITE')
        end if
        call jeveuo(chamn2//'.VALE', 'L', jchin)
        call mrmult('ZERO', lmat, zr(jchin), zr(jchout), 1, .true._1)
        if (nr .ne. 0) then
            call rsexch(' ', resui, 'DEPL', numei, chamno, iret)
            call vtcopy(chamno, chamn2, ier)
            if (ier .ne. 0) then
                call utmess("A", "FIELD0_3", sk='DEPL')
            end if
            call jeveuo(chamn2//'.VALE', 'L', jchin)
            call mrmult('CUMU', lma2, zr(jchin), zr(jchout), 1, .true._1)
        end if
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

        call rsnoch(resu, nsymb, icompt)
        call rsadpa(resu, 'E', 1, typabs, icompt, 0, sjv=iad, styp=k8b)
        zr(iad) = tps
        call rssepa(resu, icompt, modele, materi, carele, list_load)
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
            matric(1) = ' '
            matric(2) = ' '
            matric(3) = mamor
            call getvid('KUCV', 'MATR_RIGI', iocc=iocc, scal=matr, nbret=n1)
            if (n1 .eq. 1) then
                matric(1) = matr
            end if
            call refdaj('F', resu19, (nbordr2-nbordr1), numem, 'DYNAMIQUE', &
                        matric, ier)
        end if
    end if
!
    call jedema()
end subroutine
