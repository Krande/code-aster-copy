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
subroutine ritz99(nomres)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/moco99.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: nomres
!
! --------------------------------------------------------------------------------------------------
!
! CREATION D'UNE BASE MODALE DE TYPE RITZ (C A D QUELCONQUE)
!
! --------------------------------------------------------------------------------------------------
!
! NOMRES /I/ : NOM K8 DU RESULTAT
!
! --------------------------------------------------------------------------------------------------
!
    complex(kind=8) :: cbid
    integer(kind=8) :: vali(3)
    character(len=8) :: resu_basemo, resu_modeintf, k8b, intf, listam
    character(len=19) :: numref, nume1
    character(len=24) :: trang1, trang2, tempor, tempi, tempi2, rigi1, mass1, amor1
    character(len=24) :: valk(3), concep(3)
    aster_logical :: seul, has_numref
    integer(kind=8) :: iocc_ritz, iam, iamog, iamor, nbocc_modeintf, ibi4, ibi5
    integer(kind=8) :: ibi6, ibid, nbocc_basemo, idgl, idiff, idor, ier
    integer(kind=8) :: ii, inord, ioci, jamo2, jamog
    integer(kind=8) :: lnbm, lrang1, lrang2, n, nam
    integer(kind=8) :: nbamor, nbdef, nbg, nbgl, nbi, nbid, nbli, nbocc_ritz
    integer(kind=8) :: nbmod1, nbmod2, nbmoda, nbmodb, nbold(1), nbtot
    integer(kind=8) :: nbmm, nbbm, nbmi, iocc_modeintf, iocc_basemo
    integer(kind=8) :: nnum
    real(kind=8) :: bid, ebid
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    seul = ASTER_FALSE
    nbmod1 = 0
    nbmod2 = 0
    tempor = '&&RITZ99.GLOBAL'
    trang1 = '&&RITZ99.NUME.RANG1'
    trang2 = '&&RITZ99.NUME.RANG2'
    tempi = '&&RITZ99.LISTE'
    tempi2 = '&&RITZ99.LISTE2'
    concep(1) = ' '
!
! - RECUPERATION NUMEROTATION DE REFERENCE
!
    ! Test wether NUME_REF has been filled by the user
    call getvid('  ', 'NUME_REF', nbval=0, iocc=1, nbret=nnum)
    has_numref = (nnum /= 0)
    if (has_numref) then
        call getvid('  ', 'NUME_REF', iocc=1, scal=numref)
    else
        call dismoi('NUME_DDL', nomres, 'RESU_DYNA', repk=numref, arret='F', ier=ier)
    end if
!
! - DETERMINATION DU NOMBRE D'OCCURRENCE DE RITZ
!
    call getfac('RITZ', nbocc_ritz)
    ASSERT(nbocc_ritz .le. 2)
    iocc_basemo = 1
    iocc_modeintf = 1
    nbocc_modeintf = 0
    do iocc_ritz = 1, nbocc_ritz
        call getvid('RITZ', 'MODE_MECA', iocc=iocc_ritz, nbval=0, nbret=nbmm)
        call getvid('RITZ', 'BASE_MODALE', iocc=iocc_ritz, nbval=0, nbret=nbbm)
        call getvid('RITZ', 'MODE_INTF', iocc=iocc_ritz, nbval=0, nbret=nbmi)
        nbocc_modeintf = nbocc_modeintf-nbmi
        if (nbmi .ne. 0) then
            iocc_modeintf = iocc_ritz
        end if
        if ((nbmm .ne. 0) .or. (nbbm .ne. 0)) then
            iocc_basemo = iocc_ritz
        end if
    end do
!
    if (nbocc_ritz .eq. 2) then
        if (nbocc_modeintf .ne. 1) then
            call utmess('F', 'DEFIBASEMODALE1_51')
        end if
    else
        if (nbmm .eq. 0) then
            call utmess('F', 'DEFIBASEMODALE1_1')
        end if
    end if
!
! - DETERMINATION DU NOMBRE DE CONCEPT(S) MODE_* de BASE_MODALE
!
    nbocc_basemo = 0
    call getvid('RITZ', 'BASE_MODALE', iocc=iocc_basemo, scal=resu_basemo, nbret=nbocc_basemo)
!
! - DETERMINATION DU NOMBRE DE CONCEPT(S) MODE_* de MODE_INTF
!
    nbocc_modeintf = 0
    call getvid('RITZ', 'MODE_INTF', iocc=iocc_modeintf, scal=resu_modeintf, nbret=nbocc_modeintf)
!
! - DEBUT DE LA BOUCLE DE TRAITEMENT DE "BASE_MODALE"
!
    ! "base_modale" refers to a mode_meca concept previously produced by DEFI_BASE_MODALE
    !      so it may contain both dynamical and statical modes
    ! If the first modes to put into the result are defined by a previous "base_modale"
    if (nbocc_basemo .ne. 0) then
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_modeintf, scal=nbmod2, nbret=ibi4)
        call rsorac(resu_modeintf, 'LONUTI', ibid, bid, k8b, &
                    cbid, ebid, 'ABSOLU', nbold, 1, &
                    nbid)
        !
        ! nbmodb: number of modes to be ADDED
        if (ibi4 .eq. 0) then
            nbmodb = nbold(1)
        else
            nbmodb = min(nbmod2, nbold(1))
        end if
        !
        ! nbmod1: number of modes in initial base_modale
        ! nbdef: number of statical modes in initial base_modale
        call dismoi('NB_MODES_TOT', resu_basemo, 'RESULTAT', repi=nbmod1)
        call dismoi('NB_MODES_STA', resu_basemo, 'RESULTAT', repi=nbdef)
! ----- TOTAL NUMBER OF MODES IN OUTPUT CONCEPT
        nbtot = nbmod1+nbmodb
        if (nbtot .le. 0) then
            call utmess('F', 'DEFIBASEMODALE1_50')
        end if
! ----- ALLOCATION DE LA STRUCTURE DE DONNEES BASE_MODALE
        if (nomres .ne. resu_basemo) then
            call rscrsd('G', nomres, 'MODE_MECA', nbtot)
        else
            ! Nombre de numeros d'ordre dans la "base_modale" d'entrée
            call rsorac(resu_basemo, 'LONUTI', ibid, bid, k8b, &
                        cbid, ebid, 'ABSOLU', nbold, 1, &
                        nbid)
            if (nbtot .gt. nbold(1)) call rsagsd(nomres, nbtot)
            if (.not. has_numref) then
                call utmess('E', 'DEFIBASEMODALE1_9')
            end if
            ! Retrieve NUME_EQUA from NUME_DDL
            numref(15:19) = '.NUME'
            !
            intf = ' '
            call getvid('  ', 'INTERF_DYNA', iocc=1, nbval=0, nbret=ioci)
            if (ioci .lt. 0) then
                call getvid('  ', 'INTERF_DYNA', iocc=1, scal=intf, nbret=ioci)
            end if
        end if
        !
        ! All present modes in initial "base_modale" are copied and set to the
        !    right numbering using moco99
        if (nbmod1 .gt. 0) then
            call wkvect(trang1, 'V V I', nbmod1, lrang1)
            do ii = 1, nbmod1
                zi(lrang1+ii-1) = ii
            end do
            inord = 1
            call moco99(nomres, resu_basemo, nbmod1, zi(lrang1), inord, &
                        .true._1)
            call jedetr(trang1)
        end if
        !
        ! All modes to be added are copied and set to the right numbering
        !    using moco99
        if (nbmodb .gt. 0) then
            call wkvect(trang2, 'V V I', nbmodb, lrang2)
            do ii = 1, nbmodb
                zi(lrang2+ii-1) = ii
            end do
            call moco99(nomres, resu_modeintf, nbmodb, zi(lrang2), inord, &
                        .false._1)
            call jedetr(trang2)
        end if
        ! Number of dynamical modes in initial "base_modale"
        nbmoda = nbmod1-nbdef
        ! = total number of dynamical modes
        ! (only statical modes may be added)
        ! Total number of statical modes
        nbmodb = nbmodb+nbdef
        ! All work has been done
        goto 40
    end if
!
! - DETERMINATION DU NOMBRE DE CONCEPT(S) MODE_MECA
!
    ! nbgl: number of mode_meca concepts
    call getvid('RITZ', 'MODE_MECA', iocc=iocc_basemo, nbval=0, nbret=nbgl)
    nbgl = -nbgl
    if (nbgl .eq. 1) then
        call getvid('RITZ', 'MODE_MECA', iocc=iocc_basemo, scal=resu_basemo)
    end if
    if (nbgl .gt. 1) then
        call wkvect(tempor, 'V V K8', nbgl, idgl)
        call wkvect(tempi, 'V V I', nbgl, idor)
!  ---ON RECUPERE ICI LE NB DE VAL DE LA LISTE NMAX_MODE
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_basemo, nbval=0, nbret=nbli)
        nbli = -nbli
        if ((nbli .ne. 0) .and. (nbli .ne. nbgl)) then
            vali(1) = nbgl
            vali(2) = nbli
            call utmess('F', 'DEFIBASEMODALE1_31', ni=2, vali=vali)
        end if
        call getvid('RITZ', 'MODE_MECA', iocc=iocc_basemo, nbval=nbgl, vect=zk8(idgl), nbret=nbg)
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_basemo, nbval=nbli, vect=zi(idor), nbret=nbi)
    end if
!
!
! --- DETERMINATION NOMBRE ET NUMERO ORDRE MODE
!
    if (nbocc_modeintf .eq. 0) then
! ----ON N'A QUE MODE_MECA ET PAS DE MODES D'INTERFACE------
        seul = .true.
    end if
!
    ! nbmoda: number of dynamical modes in output
    if (nbgl .eq. 1) then
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_basemo, scal=nbmod1, nbret=ibi5)
        nbmoda = nbmod1
        call rsorac(resu_basemo, 'LONUTI', ibid, bid, k8b, &
                    cbid, ebid, 'ABSOLU', nbold(1), 1, &
                    nbid)
        if (ibi5 .eq. 0) then
            nbmoda = nbold(1)
        else
            nbmoda = min(nbmod1, nbold(1))
        end if
    else if (nbgl .gt. 1) then
        nbmoda = 0
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_basemo, scal=nbmod1, nbret=ibi5)
        call wkvect(tempi2, 'V V I', nbgl, lnbm)
        do iocc_ritz = 1, nbgl
            call rsorac(zk8(idgl+iocc_ritz-1), 'LONUTI', ibid, bid, k8b, &
                        cbid, ebid, 'ABSOLU', nbold, 1, &
                        nbid)
            if (ibi5 .eq. 0) then
                nbmoda = nbmoda+nbold(1)
                zi(lnbm+iocc_ritz-1) = nbold(1)
            else
                nbmoda = nbmoda+min(zi(idor+iocc_ritz-1), nbold(1))
                zi(lnbm+iocc_ritz-1) = min(zi(idor+iocc_ritz-1), nbold(1))
            end if
        end do
    end if
!
    if (nbmoda .gt. 0) then
        call wkvect(trang1, 'V V I', nbmoda, lrang1)
        do ii = 1, nbmoda
            zi(lrang1+ii-1) = ii
        end do
    end if
!
    ! if there are statical modes
    ! nbmodb: number of statical modes in output
    if (.not. seul) then
        call getvis('RITZ', 'NMAX_MODE', iocc=iocc_modeintf, scal=nbmod2, nbret=ibi6)
        call rsorac(resu_modeintf, 'LONUTI', ibid, bid, k8b, &
                    cbid, ebid, 'ABSOLU', nbold, 1, &
                    nbid)
        if (ibi6 .eq. 0) then
            nbmodb = nbold(1)
        else
            nbmodb = min(nbmod2, nbold(1))
        end if
        if (nbmodb .gt. 0) then
            call wkvect(trang2, 'V V I', nbmodb, lrang2)
            do ii = 1, nbmodb
                zi(lrang2+ii-1) = ii
            end do
        end if
    else
        nbmodb = 0
    end if
!
!
! --- DETERMINATION NOMBRE TOTAL
!
    nbtot = nbmoda+nbmodb
!
! --- ON AJOUTE LA LIST_AMOR--------------------------------------
!
    call getvid(' ', 'LIST_AMOR', nbval=0, nbret=nam)
    if (nam .ne. 0) then
        call getvid(' ', 'LIST_AMOR', scal=listam, nbret=n)
        call jelira(listam//'           .VALE', 'LONMAX', nbamor)
        if (nbamor .gt. nbmoda) then
            vali(1) = nbmoda
            vali(2) = nbamor
            vali(3) = nbmoda
            valk(1) = 'PREMIERS COEFFICIENTS'
            call utmess('A', 'ALGORITH16_18', sk=valk(1), ni=3, vali=vali)
            call wkvect('&&RITZ99.AMORTI', 'V V R8', nbmoda, jamog)
!
            call jeveuo(listam//'           .VALE', 'L', iamog)
            do iam = 1, nbmoda
                zr(jamog+iam-1) = zr(iamog+iam-1)
            end do
        else if (nbamor .lt. nbmoda) then
            call wkvect('&&RITZ99.AMORTI', 'V V R8', nbamor, jamog)
            call jeveuo(listam//'           .VALE', 'L', iamog)
            do iam = 1, nbamor
                zr(jamog+iam-1) = zr(iamog+iam-1)
            end do
            idiff = nbmoda-nbamor
            vali(1) = idiff
            vali(2) = nbmoda
            vali(3) = idiff
            call utmess('I', 'ALGORITH16_19', ni=3, vali=vali)
            call wkvect('&&RITZ99.AMORTI2', 'V V R8', nbmoda, jamo2)
            do iam = 1, nbamor
                zr(jamo2+iam-1) = zr(jamog+iam-1)
            end do
            do iam = nbamor+1, nbmoda
                zr(jamo2+iam-1) = zr(jamog+nbamor-1)
            end do
            jamog = jamo2
        else if (nbamor .eq. nbmoda) then
            call wkvect('&&RITZ99.AMORTI', 'V V R8', nbamor, jamog)
            call jeveuo(listam//'           .VALE', 'L', iamog)
            do iam = 1, nbamor
                zr(jamog+iam-1) = zr(iamog+iam-1)
            end do
        end if
!   ----ON AJOUTE LA LIST_AMOR COMME VALEURS DU PARAM 'AMOR_REDUIT'
!       DU RESULT1 (SI UN SEUL MODE_MECA)
        if (nbgl .eq. 1) then
            do iam = 1, nbmoda
                call rsadpa(resu_basemo, 'E', 1, 'AMOR_REDUIT', iam, &
                            0, sjv=iamor)
                zr(iamor) = zr(jamog+iam-1)
            end do
        end if
    end if
!
!
! --- ALLOCATION DE LA STRUCTURE DE DONNEES BASE_MODALE
!
    if (nbtot .gt. 0) then
        call rscrsd('G', nomres, 'MODE_MECA', nbtot)
    else
        call utmess('F', 'DEFIBASEMODALE1_50')
    end if
!
! --- COPIE DES MODES DYNAMIQUES EN CORRIGEANT LA NUMEROTATION (MOCO99)
!
    inord = 1
    if (nbmoda .gt. 0) then
        if (nbgl .eq. 1) then
            call moco99(nomres, resu_basemo, nbmoda, zi(lrang1), inord, &
                        ASTER_TRUE)
!
            call dismoi('NUME_DDL', resu_basemo, 'RESU_DYNA', repk=nume1)
            call dismoi('REF_RIGI_PREM', resu_basemo, 'RESU_DYNA', repk=rigi1, arret='C', ier=ier)
            call dismoi('REF_MASS_PREM', resu_basemo, 'RESU_DYNA', repk=mass1, arret='C', ier=ier)
            call dismoi('REF_AMOR_PREM', resu_basemo, 'RESU_DYNA', repk=amor1, arret='C', ier=ier)
        else if (nbgl .gt. 1) then
            do iocc_ritz = 1, nbgl
                call moco99(nomres, zk8(idgl+iocc_ritz-1), &
                            zi(lnbm+iocc_ritz-1), zi(lrang1), inord, &
                            ASTER_TRUE)
                resu_basemo = zk8(idgl+iocc_ritz-1)
            end do
            inord = inord+nbmoda
        end if
!
        call jedetr(trang1)
    end if
!
! --- COPIE DES MODES STATIQUES EN CORRIGEANT LA NUMEROTATION (MOCO99)
!
    if (.not. seul) then
        if (nbmodb .gt. 0) then
            call moco99(nomres, resu_modeintf, nbmodb, zi(lrang2), inord, &
                        ASTER_FALSE)
            call jedetr(trang2)
        end if
    end if
!
40  continue
!
    call jedetr(tempor)
    call jedetr(trang1)
    call jedetr(trang2)
    call jedetr(tempi)
    call jedetr(tempi2)
    call jedetr('&&RITZ99.AMORTI')
!
    call jedema()
end subroutine
