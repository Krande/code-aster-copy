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
subroutine crcoch()
!
    use listLoad_module
!
    implicit none
!
!     COMMANDE:  CREA_RESU /CONV_CHAR
!     CREE UNE STRUCTURE DE DONNEE DE TYPE
!           "DYNA_TRANS"   "EVOL_CHAR"
!
! --- ------------------------------------------------------------------
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asasve.h"
#include "asterfort/ascova.h"
#include "asterfort/copisd.h"
#include "asterfort/crcoch_getloads.h"
#include "asterfort/dismoi.h"
#include "asterfort/fondpl.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcmfmc.h"
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
#include "asterfort/vechme.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: ibid, ier, icompt, iret, numini, numfin
    integer(kind=8) :: n1, nis, nbinst, nbval, nume, j, ie
    integer(kind=8) :: iad, jinst, jchout
    integer(kind=8) :: nbv(1), jrefe
    integer(kind=8) :: jcpt, nbr, ivmx, k, iocc, nboini
    integer(kind=8) :: nb_ondp, nbLoad, tnum(1)
    integer(kind=8) :: nbordr1, nbordr2
    integer(kind=8) :: neq, jchou1, jchou2
    real(kind=8) :: rbid, tps, prec, partps(3)
    complex(kind=8) :: cbid
    character(len=1) :: typmat
    character(len=4) :: typabs
    character(len=8) :: k8b, resu, criter, matr
    character(len=8) :: materField, caraElem, blan8
    character(len=16) :: type, oper
    character(len=19) :: nomch, listr8, resu19, profch
    character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
    character(len=24) :: listLoad, listLoadResu
    character(len=19) :: nomch1, nomch2
    character(len=24) :: linst, nsymb, typres, lcpt, londp
    character(len=24) :: matric(3)
    character(len=24) :: model, mate, numedd, vecond
    character(len=24) :: veonde, vaonde
    character(len=24) :: vechmp, vachmp, cnchmp, chargt
    real(kind=8), pointer :: val(:) => null()
    character(len=8), pointer :: v_ondp(:) => null()
!
    data linst, listr8, lcpt, londp/'&&CRCOCH_LINST', '&&CRCOCH_LISR8',&
     &     '&&CPT_CRCOCH', '&&CRCOCH_LONDP'/
! --- ------------------------------------------------------------------
    call jemarq()
!
    blan8 = ' '
    nboini = 10
    model = ' '
    caraElem = ' '
    materField = ' '
    vecond = '&&CRCOCH_VECOND'
    veonde = '&&CRCOCH_VEONDE'
    vaonde = '&&CRCOCH_VAONDE'
    listLoad = '&&CRCOCH.LISCHA'
    listLoadResu = ' '
    vechmp = '&&CRCOCH_VECHMP'
    vachmp = '&&CRCOCH_VACHMP'
    cnchmp = '&&CRCOCH_CNCHMP'
    chargt = '&&CRCOCH_CHARGT'
    nomch1 = '&&CRCOCH_NOMCH1'
    nomch2 = '&&CRCOCH_NOMCH2'
    numini = -1
    icompt = -1
    profch = ' '
    iocc = 1
    nbLoad = 0
!
    call getres(resu, type, oper)
    resu19 = resu
    call getvtx(' ', 'TYPE_RESU', scal=typres)

! - Get loads
    call getvid('CONV_CHAR', 'CHARGE', iocc=iocc, nbval=0, nbret=n1)
    if (n1 .ne. 0) then
! ----- Create datastructure for list of loads
        call crcoch_getloads(listLoad, nbLoad, nb_ondp, v_ondp)

! ----- Generate name of list of loads to save in results datastructure
        call nameListLoad(listLoadResu)
    end if
!
! - Create output datastructure
!
    call rscrsd('G', resu, typres, nboini)
    call jelira(resu//'           .ORDR', 'LONUTI', nbordr1)
!
!        MOT CLE INST PRESENT :
    nis = 0
    nbinst = 0
    call getvr8('CONV_CHAR', 'INST', iocc=iocc, nbval=0, nbret=nis)
    if (nis .ne. 0) then
        typabs = 'INST'
        nbinst = -nis
    end if
!
    if (nis .ne. 0) then
        call wkvect(lcpt, 'V V I', nbinst, jcpt)
        call wkvect(linst, 'V V R', nbinst, jinst)
        call getvr8('CONV_CHAR', typabs, iocc=iocc, nbval=nbinst, vect=zr(jinst))
        call getvr8('CONV_CHAR', 'PRECISION', iocc=iocc, scal=prec)
        call getvtx('CONV_CHAR', 'CRITERE', iocc=iocc, scal=criter)
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
                call utmess('F', 'CREARESU1_48')
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
        call getvid('CONV_CHAR', 'LIST_INST', iocc=iocc, scal=listr8, nbret=n1)
        if (n1 .ne. 0) then
            typabs = 'INST'
        end if
!
        call getvr8('CONV_CHAR', 'PRECISION', iocc=iocc, scal=prec, nbret=ibid)
        call getvtx('CONV_CHAR', 'CRITERE', iocc=iocc, scal=criter, nbret=ibid)
        call jelira(listr8//'.VALE', 'LONMAX', nbval)
!
        nbinst = nbval
        numini = 1
        numfin = nbinst
!
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
            if (k .lt. numini) cycle
            if (k .gt. numfin) cycle
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
                call utmess('F', 'CREARESU1_48')
            else if (nbr .eq. 0) then
                zi(jcpt+j-1) = ivmx+1
                ivmx = ivmx+1
            else
                zi(jcpt+j-1) = nume
            end if
        end do
    end if
!
    numedd = ' '
    call getvid('CONV_CHAR', 'MATR_RIGI', iocc=iocc, scal=matr, nbret=n1)
    if (n1 .eq. 1) then
        call dismoi('NOM_NUME_DDL', matr, 'MATR_ASSE', repk=numedd)
    end if
!
    call dismoi('NOM_MODELE', matr, 'MATR_ASSE', repk=model)
    call dismoi('NUME_EQUA', matr, 'MATR_ASSE', repk=profch)
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    call getvid('CONV_CHAR', 'CHAM_MATER', iocc=iocc, scal=materField, nbret=n1)
    caraElem = '        '
    call getvid('CONV_CHAR', 'CARA_ELEM', iocc=iocc, scal=caraElem, nbret=n1)
    call rcmfmc(materField, mate, l_ther_=ASTER_FALSE)
    typmat = 'R'
    if (typres(1:10) .eq. 'DYNA_TRANS') then
        nsymb = 'DEPL'
    else
        nsymb = 'FORC_NODA'
    end if
    call rsagsd(resu, nbinst)
!
    do j = 1, nbinst
        if (j .ge. 2) call jemarq()
        call jerecu('V')
        icompt = zi(jcpt+j-1)
        tps = zr(jinst+j-1)
        partps(1) = tps
        partps(2) = r8vide()
        partps(3) = r8vide()
        call rsexch(' ', resu, nsymb, icompt, nomch, &
                    iret)
        if (iret .eq. 0) then
            call rsadpa(resu, 'L', 1, typabs, icompt, &
                        0, sjv=iad)
        else if (iret .eq. 110) then
            call rsagsd(resu, 0)
            call rsexch(' ', resu, nsymb, icompt, nomch, &
                        iret)
        else if (iret .eq. 100) then
            call vtcreb(nomch, 'G', 'R', nume_ddlz=numedd)
        end if
        call jeveuo(nomch//'.VALE', 'E', jchout)
        do ie = 1, neq
            zr(jchout+ie-1) = 0.d0
        end do
        if (nb_ondp .ne. 0) then
            call jeexin(nomch1//'.VALE', iret)
            if (iret .eq. 0) then
                call vtcreb(nomch1, 'V', 'R', nume_ddlz=numedd)
            end if
            call jeveuo(nomch1//'.VALE', 'E', jchou1)
            do ie = 1, neq
                zr(jchou1+ie-1) = 0.d0
            end do
            call fondpl(model, materField, mate, numedd, neq, &
                        v_ondp, nb_ondp, vecond, veonde, vaonde, &
                        tps, zr(jchou1))
            do ie = 1, neq
                zr(jchout+ie-1) = zr(jchout+ie-1)-zr(jchou1+ie-1)
            end do
        end if
        if (nbLoad .ne. 0) then
            loadNameJv = listLoad(1:19)//'.LCHA'
            loadInfoJv = listLoad(1:19)//'.INFC'
            loadFuncJv = listLoad(1:19)//'.FCHA'
            call jeexin(nomch2//'.VALE', iret)
            if (iret .eq. 0) then
                call vtcreb(nomch2, 'V', 'R', nume_ddlz=numedd)
            end if
            call jeveuo(nomch2//'.VALE', 'E', jchou2)
            do ie = 1, neq
                zr(jchou2+ie-1) = 0.d0
            end do
            call vechme('S', &
                        model, caraElem, materField, mate, &
                        loadNameJv, loadInfoJv, &
                        partps, &
                        vechmp)
            call asasve(vechmp, numedd, 'R', vachmp)
            call ascova('D', vachmp, loadFuncJv, 'INST', tps, &
                        'R', nomch2)
            call jeveuo(nomch2//'.VALE', 'L', jchou2)
            do ie = 1, neq
                zr(jchout+ie-1) = zr(jchout+ie-1)+zr(jchou2+ie-1)
            end do
        end if
!
        call jeveuo(nomch//'.REFE', 'E', jrefe)
        zk24(jrefe-1+1) = " "
        zk24(jrefe-1+2) = profch
!
        call rsnoch(resu, nsymb, icompt)
        call rsadpa(resu, 'E', 1, typabs, icompt, &
                    0, sjv=iad)
        zr(iad) = tps

! ----- Copy of list of loads for save in results datastructure
        call copisd('LISTE_CHARGES', 'G', listLoad, listLoadResu)

! ----- Save parameters in results datastructure
        call rssepa(resu, icompt, model, materField, caraElem, listLoadResu)
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
!
    AS_DEALLOCATE(vk8=v_ondp)
!
    call jedema()
end subroutine
