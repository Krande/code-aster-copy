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
subroutine rstran(interp, resu, motcle, iocc, kdisc,&
                  krang, nbdisc, ier)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsindi.h"
#include "asterfort/rslipa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: interp, motcle
    character(len=*) :: resu, kdisc, krang
! person_in_charge: jacques.pellet at edf.fr
!
!     POUR INTERP = 'NON'
!        RECUPERATION DES DISCRETISATIONS ET DES NUMEROS DE RANGEMENT
!        ASSOCIES DANS LA STRUCTURE DE DONNEES "RESU_"
!     POUR INTERP = 'LIN', 'LOG', ...
!        RECUPERATION DES INSTANTS UTILISATEURS
!     ------------------------------------------------------------------
! IN  : INTERP : TYPE D'INTERPOLATION
! IN  : RESU_   : NOM DE LA STRUCTURE DE DONNEES
! IN  : MOTCLE : MOT CLE FACTEUR
! IN  : IOCC   : NUMERO D'OCCURENCE
! IN  : KDISC_  : NOM JEVEUX POUR STOCKER LES INSTANTS
! IN  : KRANG_  : NOM JEVEUX POUR STOCKER LES NUMEROS DE RANGEMENT
! OUT : NBDISC : NOMBRE D'INSTANTS/FREQUENCES TROUVES
! OUT : IER    : CODE RETOUR, = 0    : OK
!                             = 100  : PLUSIEURS CHAMPS TROUVES
!                             = 110  : AUCUN CHAMP TROUVE
!                             SINON  : NOOK
!     ------------------------------------------------------------------
    integer :: vali
    real(kind=8) :: valr
    character(len=4) :: type
    character(len=8) :: k8b, crit
    character(len=16) :: nomcmd
    character(len=19) :: listr
    character(len=8) :: kval
    complex(kind=8) :: cval
    character(len=19) :: resu_
    character(len=24) :: typres, kdisc_, krang_
!------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer :: i, ier, ier1, iocc, iord, iret
    integer :: ival, jdisc, jordr, jrang, l, laccr
    integer :: ldisc, lli, lt, n, nbi, nbi2, nbdisc
    integer :: nbtrou, nno, nto, nutrou(1)
    real(kind=8) :: epsi, tusr
    integer, pointer :: nume(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    resu_ = resu
    kdisc_ = kdisc
    krang_ = krang
!
    ier = 0
    nbdisc = 0
    type = 'R8  '
    call getres(k8b, k8b, nomcmd)
    call getvr8(motcle, 'PRECISION', iocc=iocc, scal=epsi, nbret=n)
    call getvtx(motcle, 'CRITERE', iocc=iocc, scal=crit, nbret=n)
!
    call jeexin(resu_//'.DISC', ier1)
    if (ier1 .gt. 0) then
!     --- CAS D'UNE SD DYNA_GENE (HARM_GENE OU TRAN_GENE)
        call jeveuo(resu_//'.DISC', 'L', ldisc)
        call jelira(resu_//'.DISC', 'LONMAX', nbi)
    else
!     --- CAS D'UNE SD RESU_LTAT
        call rslipa(resu_, 'INST', '&&RSTRAN.LIINST', ldisc, nbi)
    endif
!
    call jeexin(resu_//'.ORDR', iret)
    if (iret .eq. 0) then
        call utmess('F', 'ALGORITH17_26')
!        CALL WKVECT('&&RSTRAN.ORDR','V V I',NBI,JORDR)
!        DO 10 I = 1,NBI
!          ZI(JORDR+I-1) = I
!   10   CONTINUE
    else
        call jeveuo(resu_//'.ORDR', 'L', jordr)
        call jelira(resu_//'.ORDR', 'LONUTI', nbi2)
        if (nbi .ne. nbi2) then
            call utmess('F', 'ALGORITH17_27')
        endif
    endif
!
!     --- RECHERCHE A PARTIR D'UN NUMERO D'ORDRE ---
!
    call getvis(motcle, 'NUME_ORDRE', iocc=iocc, nbval=0, nbret=nno)
    if (nno .ne. 0) then
        nbdisc = -nno
        call wkvect(krang_, 'V V I', nbdisc, jrang)
        call wkvect(kdisc_, 'V V R8', nbdisc, jdisc)
        AS_ALLOCATE(vi=nume, size=nbdisc)
        call getvis(motcle, 'NUME_ORDRE', iocc=iocc, nbval=nbdisc, vect=nume,&
                    nbret=nno)
        do i = 0, nbdisc - 1
            do iord = 0, nbi - 1
                if (nume(1+i) .eq. zi(jordr+iord)) goto 30
            end do
            ier = ier + 110
            vali = nume(1+i)
            call utmess('A', 'UTILITAI8_17', si=vali)
            goto 40
 30         continue
            zi(jrang+i) = iord + 1
            zr(jdisc+i) = zr(ldisc+iord)
 40         continue
        end do
        goto 100
    endif
!
!     --- RECHERCHE A PARTIR D'UN INSTANT ---
!
    call getvr8(motcle, 'INST', iocc=iocc, nbval=0, nbret=lt)
    if (lt .eq. 0) then
        call getvid(motcle, 'LIST_INST', iocc=iocc, scal=listr, nbret=lli)
        if (lli .ne. 0) then
            call jeveuo(listr//'.VALE', 'L', laccr)
            call jelira(listr//'.VALE', 'LONMAX', nbdisc)
        else
            goto 80
        endif
    else
        nbdisc = -lt
        call wkvect('&&RSTRAN.INSTANTS', 'V V R', nbdisc, laccr)
        call getvr8(motcle, 'INST', iocc=iocc, nbval=nbdisc, vect=zr(laccr),&
                    nbret=l)
    endif
    call wkvect(krang_, 'V V I', nbdisc, jrang)
    call wkvect(kdisc_, 'V V R8', nbdisc, jdisc)
    do i = 0, nbdisc - 1
        tusr = zr(laccr+i)
        if (interp(1:3) .ne. 'NON') then
            zi(jrang+i) = i + 1
            zr(jdisc+i) = tusr
            goto 70
        endif
        call rsindi(type, ldisc, 1, jordr, ival,&
                    tusr, kval, cval, epsi, crit,&
                    nbi, nbtrou, nutrou, 1)
        if (nbtrou .eq. 0) then
            ier = ier + 110
            valr = tusr
            call utmess('A', 'UTILITAI8_18', sr=valr)
            goto 70
        else if (nbtrou.ne.1) then
            ier = ier + 100
            valr = tusr
            vali = -nbtrou
            call utmess('F', 'UTILITAI8_19', si=vali, sr=valr)
            goto 70
        endif
        do iord = 0, nbi - 1
            if (nutrou(1) .eq. zi(jordr+iord)) goto 60
        end do
 60     continue
        zi(jrang+i) = iord + 1
        zr(jdisc+i) = zr(ldisc+iord)
 70     continue
    end do
    goto 100
!
 80 continue
!
!     --- RECHERCHE A PARTIR D'UNE FREQUENCE ---
!
    call gettco(resu_(1:8), typres)
    if (typres(1:9) .eq. 'HARM_GENE') then
        call getvr8(motcle, 'FREQ', iocc=iocc, nbval=0, nbret=lt)
        if (lt .eq. 0) then
            call getvid(motcle, 'LIST_FREQ', iocc=iocc, scal=listr, nbret=lli)
            if (lli .ne. 0) then
                call jeveuo(listr//'.VALE', 'L', laccr)
                call jelira(listr//'.VALE', 'LONMAX', nbdisc)
            else
                goto 81
            endif
        else
            nbdisc = -lt
            call wkvect('&&RSTRAN.FREQUENCES', 'V V R', nbdisc, laccr)
            call getvr8(motcle, 'FREQ', iocc=iocc, nbval=nbdisc, vect=zr(laccr),&
                        nbret=l)
        endif
        call wkvect(krang_, 'V V I', nbdisc, jrang)
        call wkvect(kdisc_, 'V V R8', nbdisc, jdisc)
        do i = 0, nbdisc - 1
            tusr = zr(laccr+i)
            call rsindi(type, ldisc, 1, jordr, ival,&
                        tusr, kval, cval, epsi, crit,&
                        nbi, nbtrou, nutrou, 1)
            if (nbtrou .eq. 0) then
                ier = ier + 110
                valr = tusr
                call utmess('A', 'UTILITAI8_18', sr=valr)
                goto 71
            else if (nbtrou.ne.1) then
                ier = ier + 100
                valr = tusr
                vali = -nbtrou
                call utmess('F', 'UTILITAI8_19', si=vali, sr=valr)
                goto 71
            endif
            do iord = 0, nbi - 1
                if (nutrou(1) .eq. zi(jordr+iord)) goto 61
            end do
 61         continue
            zi(jrang+i) = iord + 1
            zr(jdisc+i) = zr(ldisc+iord)
 71         continue
        end do
        goto 100
    endif
!
!  --- PAR DEFAUT, TOUT ORDRE
!
 81 continue
!
    call getvtx(motcle, 'TOUT_INST', iocc=iocc, scal=k8b, nbret=nto)
    call getvtx(motcle, 'TOUT_ORDRE', iocc=iocc, scal=k8b, nbret=nto)
    nbdisc = nbi
    call wkvect(krang_, 'V V I', nbdisc, jrang)
    call wkvect(kdisc_, 'V V R8', nbdisc, jdisc)
    do iord = 0, nbdisc - 1
        zi(jrang+iord) = iord + 1
        zr(jdisc+iord) = zr(ldisc+iord)
    end do
!
100 continue
    call jedetr('&&RSTRAN.ORDR')
    AS_DEALLOCATE(vi=nume)
    call jedetr('&&RSTRAN.INSTANTS')
    call jedetr('&&RSTRAN.FREQUENCES')
    call jedetr('&&RSTRAN.LIINST')
!
    call jedema()
end subroutine
