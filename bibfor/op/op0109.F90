! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine op0109()
    implicit none
!
!     COMMANDE : COMB_SISM_MODAL
!
!     ------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/ascalc.h"
#include "asterfort/asenap.h"
#include "asterfort/asexci.h"
#include "asterfort/asimpr.h"
#include "asterfort/asmsup.h"
#include "asterfort/assert.h"
#include "asterfort/asveri.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbliva.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/vprecu.h"
#include "asterfort/wkvect.h"
#include "asterfort/jexnum.h"
!
    integer :: vali(2)
!-----------------------------------------------------------------------
    integer :: iam, ibid, id, ifu, ii
    integer :: im, info, iret, j, jamo2
    integer :: jamor, jarm, jasy, jcsu, jdep, jdir
    integer :: jkno, jnsu, jopt, jordr, jrea, jspe
    integer :: lmod, lval, na, na1, nbamor
    integer :: nbfac, nbmode, nbopt, nbordr, nbpara, nbpari, nbpark
    integer :: nbparr, nbsup, nc, ncd, ncm, ncmt, ndepl
    integer :: neq, nf1, nf2, nimpr, nla, nmm, nmult
    integer :: nna, np, npm, nret1, nret2, ns, nt
    integer :: nty2, jomeg, jmasg
    real(kind=8) :: rundef, xcumul, xfm
!-----------------------------------------------------------------------
    parameter(nbpara=11)
    real(kind=8) :: r8b, temps, prec, xmastr, masuni, f1gup, f2gup
    real(kind=8) :: freq, facpar, masmod, mastot(3), zero, cumul(3)
    character(len=1) :: dir(3)
    character(len=3) :: corf
    character(len=4) :: ctyp
    character(len=8) :: k8b, resu, meca, psmo, stat, masse, typcmo, typcdi
    character(len=8) :: crit, tmas, noma
    character(len=8) :: nature, typcma, paraki(2), valeki(2)
    character(len=9) :: niveau
    character(len=16) :: nomcmd, concep, nomsy
    character(len=14) :: nume
    character(len=19) :: kvec, kval, kspect, kasysp, knoeu, knume
    character(len=19) :: liar
    character(len=24) :: desc, refd, nopara(nbpara)
    character(len=24) :: valk(3)
    aster_logical :: tronc, monoap, muapde, comdir, corfre, calmas
    complex(kind=8) :: c16b
!     ------------------------------------------------------------------
    data dir/'X', 'Y', 'Z'/
    data desc/'                   .SCDE'/
    data refd/'                   .REFD'/
    data kvec/'&&OP0109.VAL_PROPRE'/
    data kval/'&&OP0109.GRAN_MODAL'/
    data kspect/'&&OP0109.SPECTRE   '/
    data kasysp/'&&OP0109.GAMMA0'/
    data knoeu/'&&OP0109.NOM_SUPPOR'/
    data knume/'&&OP0109.NUME_ORDRE'/
    data nopara/&
     &  'FREQ', 'OMEGA2',&
     &  'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
     &  'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
     &  'MASS_EFFE_UN_DX', 'MASS_EFFE_UN_DY', 'MASS_EFFE_UN_DZ'/
!     ------------------------------------------------------------------
    call jemarq()
    r8b = 0.d0
    zero = 0.d0
    tronc = .false.
    comdir = .false.
    rundef = r8vide()
!
    call getres(resu, concep, nomcmd)
!
!     --- LECTURE MOT-CLE FACTEUR IMPRESSION ---
!
    call getvtx('IMPRESSION', 'NIVEAU', iocc=1, scal=niveau, nbret=nimpr)
    if (nimpr .eq. 0) niveau = 'TOUT     '
!
!     ----- RECUPERATION DES OPTIONS DE CALCUL -----
!
    call getvtx(' ', 'OPTION', nbval=0, nbret=ns)
    nbopt = -ns
    call wkvect('&&OP0109.OPTION', 'V V K16', nbopt, jopt)
    call getvtx(' ', 'OPTION', nbval=nbopt, vect=zk16(jopt), nbret=ns)
!
!     ----- RECUPERATION DES MODES -----
!
    call getvid(' ', 'MODE_MECA', scal=meca, nbret=nmm)
    call getvid(' ', 'MODE_CORR', scal=psmo, nbret=npm)
    if (npm .ne. 0) tronc = .true.
!
    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(meca, ' ', 0, knume, nbordr, &
                prec, crit, iret)
    if (iret .ne. 0) goto 999
    call jeveuo(knume, 'L', jordr)
    call dismoi('REF_MASS_PREM', meca, 'RESU_DYNA', repk=masse, arret='C')
    nomsy = 'DEPL'
    call vprecu(meca, nomsy, nbordr, zi(jordr), kvec, &
                nbpara, nopara(1), k8b, kval, k8b, &
                neq, nbmode, ctyp, nbpari, nbparr, &
                nbpark)
!
    call jeveuo(kvec, 'L', lmod)
    call jeveuo(kval, 'L', lval)
!
!    --- ON TESTE SI LES PARAMATRES REELS SONT BIEN PRESENTS
!        LE TEST CONSISTE A VERIFIER QUE MASS_EFFE_DX DU 1ER MODE
!        A UNE VALEUR REELE DIFFERENTE DE R8MAEM
!
    if (zr(lval+nbmode*2) .eq. r8vide()) then
!
        valk(1) = meca
        call utmess('F', 'SEISME_27', sk=valk(1))
!
    end if
!
!     ----- RECUPERATION DES AMORTISSEMENTS -----
! Traitement du cas AMOR_REDUIT
    call getvr8(' ', 'AMOR_REDUIT', nbval=0, nbret=na1)
    na = na1
    if (na .ne. 0) then
        nbamor = -na
        call wkvect('&&OP0109.AMORTISSEMENT', 'V V R', nbamor, jamor)
        if (na1 .ne. 0) then
            call getvr8(' ', 'AMOR_REDUIT', nbval=nbamor, vect=zr(jamor), nbret=na)
        end if
        if (nbamor .gt. nbmode) then
            vali(1) = nbamor
            vali(2) = nbmode
            call utmess('F', 'SEISME_11', ni=2, vali=vali)
        end if
        if (nbamor .lt. nbmode) then
            call wkvect('&&OP0109.AMORTISSEMEN2', 'V V R', nbmode, jamo2)
            do iam = 1, nbamor
                zr(jamo2+iam-1) = zr(jamor+iam-1)
            end do
            do iam = nbamor, nbmode
                zr(jamo2+iam-1) = zr(jamor+nbamor-1)
            end do
            nbamor = nbmode
            jamor = jamo2
!
        end if
    else
! Traitement du cas LIST_AMOR
        call getvid(' ', 'LIST_AMOR', scal=liar, nbret=nla)
        if (nla .ne. 0) then
            call jelira(liar//'.VALE', 'LONUTI', nbamor)
            if (nbamor .gt. nbmode) then
                vali(1) = nbamor
                vali(2) = nbmode
                call utmess('F', 'SEISME_11', ni=2, vali=vali)
            end if
            call jeveuo(liar//'.VALE', 'L', jarm)
            call wkvect('&&OP0109.AMORTISSEMENT', 'V V R', nbmode, jamor)
            do iam = 1, nbamor
                zr(jamor+iam-1) = zr(jarm+iam-1)
            end do
            if (nbamor .lt. nbmode) then
                do iam = nbamor, nbmode
                    zr(jamor+iam-1) = zr(jarm+nbamor-1)
                end do
            end if
            nbamor = nbmode
        else
!Traitement du cas AMOR_GENE
            call getvid(' ', 'AMOR_GENE', scal=liar, nbret=nla)
            if (nla .ne. 0) then
                call jeveuo(liar//'.DESC', 'L', jamor)
!               On verifie que AMOR_GENE est une matrice (+ diagonale)
                if (zi(jamor-1+1) .ne. 2) then
                    ASSERT(.false.)
                end if
                if (zi(jamor-1+3) .ne. 1) then
                    ASSERT(.false.)
                end if
                call jelira(liar//'.VALM', 'LONO', nbamor)
                if (nbamor .gt. nbmode) then
                    vali(1) = nbamor
                    vali(2) = nbmode
                    call utmess('F', 'SEISME_11', ni=2, vali=vali)
                end if
                call jeveuo(jexnum(liar//'.VALM', 1), 'L', jarm)
                call wkvect('&&OP0109.AMORTISSEMENT', 'V V R', nbmode, jamor)
                do iam = 1, nbamor
                    zr(jamor-1+iam) = zr(jarm-1+iam)
                end do
                if (nbamor .lt. nbmode) then
                    do iam = nbamor, nbmode
                        zr(jamor+iam-1) = zr(jarm+nbamor-1)
                    end do
                end if
                nbamor = nbmode
                do iam = 1, nbmode
                    call rsadpa(meca, 'L', 1, 'MASS_GENE', iam, &
                                0, sjv=jmasg)
                    call rsadpa(meca, 'L', 1, 'OMEGA2', iam, &
                                0, sjv=jomeg)
                    zr(jamor-1+iam) = zr(jamor-1+iam)/(2*sqrt(zr(jomeg))*zr(jmasg))
                end do
            else
                ASSERT(.false.)
            end if
        end if
    end if
    if (nbamor .ne. nbmode) then
        vali(1) = nbamor
        vali(2) = nbmode
        call utmess('F', 'SEISME_13', ni=2, vali=vali)
    end if
!     ----- DIVERS RECOMBINAISON -----
    call getvtx('COMB_MODE', 'TYPE', iocc=1, scal=typcmo, nbret=ncm)
    call getvr8('COMB_MODE', 'DUREE', iocc=1, scal=temps, nbret=ncmt)
    call getvr8('COMB_MODE', 'FREQ_1', iocc=1, scal=f1gup, nbret=nf1)
    call getvr8('COMB_MODE', 'FREQ_2', iocc=1, scal=f2gup, nbret=nf2)
!
    call getvtx('COMB_DIRECTION', 'TYPE', iocc=1, scal=typcdi, nbret=ncd)
    if (ncd .ne. 0) comdir = .true.
    call getvtx('EXCIT', 'NATURE', iocc=1, scal=nature, nbret=nna)
!
    call infmaj()
    call infniv(ifu, info)
!
    corfre = .false.
    call getvtx(' ', 'CORR_FREQ', scal=corf, nbret=nc)
    if (corf .eq. 'OUI') corfre = .true.
!
    if (info .eq. 1 .or. info .eq. 2) then
        valk(1) = meca
        valk(2) = typcmo
        valk(3) = zk16(jopt)
        vali(1) = nbmode
        call utmess('I', 'SEISME_15', nk=3, valk=valk, si=vali(1))
        do j = 2, nbopt
            call utmess('I', 'SEISME_16', sk=zk16(jopt+j-1))
        end do
        if (nna .ne. 0) then
            call utmess('I', 'SEISME_17', sk=nature)
        end if
        if (ncd .ne. 0) then
            call utmess('I', 'SEISME_18', sk=typcdi)
        end if
    end if
!     ----- RECUPERATION DES EXCITATIONS -----
    call utmess('I', 'SEISME_65')
    call wkvect('&&OP0109.DIRECTION', 'V V I', 3, jdir)
    call wkvect('&&OP0109.NB_SUPPOR', 'V V I', 3, jnsu)
    call asexci(masse, meca, zr(jamor), nbmode, corfre, &
                info, zi(jdir), monoap, muapde, kspect, &
                kasysp, nbsup, zi(jnsu), knoeu, nopara, &
                zi(jordr))
    call jeveuo(kasysp, 'E', jasy)
    call jeveuo(kspect, 'E', jspe)
    if (.not. monoap) then
        call jeveuo(knoeu, 'E', jkno)
    else
        jkno = 1
    end if
!     ----- VERIFICATION DE LA COHERENCE DES REQUETES    -----
!     -----  SUR LES COMPOSANTES DANS LE CAS CORRELE     -----
    typcma = ' '
    call getfac('COMB_DEPL_APPUI', ndepl)
    if (info .eq. 1 .or. info .eq. 2) then
        if ((.not. monoap) .and. (.not. muapde)) then
            call getvtx('COMB_MULT_APPUI', 'TYPE_COMBI', iocc=1, scal=typcma, nbret=nty2)
            call getfac('COMB_MULT_APPUI', nmult)
            if (ndepl .ne. 0 .and. nmult .eq. 0) then
                call utmess('F', 'SEISME_14')
            end if
            if (nty2 .ne. 0) then
                call utmess('I', 'SEISME_19', sk=typcma)
            end if
        else if ((.not. monoap) .and. (muapde)) then
            call getfac('COMB_MULT_APPUI', nmult)
            if (nmult .ne. 0) then
                call utmess('A', 'SEISME_28')
            end if
        end if
    end if
!
!
!     ----- MASSE DE LA STRUCTURE ---
    calmas = .false.
    xmastr = 1.d0
    call getvid(' ', 'MASS_INER', scal=tmas, nbret=nt)
    if (nt .ne. 0) then
!        VERIFICATION DES PARAMETRES DE LA TABLE 'TMAS'
        call dismoi('NOM_MAILLA', masse, 'MATR_ASSE', repk=noma)
        call tbexp2(tmas, 'LIEU')
        call tbexp2(tmas, 'MASSE')
        call tbliva(tmas, 1, 'LIEU', [ibid], [r8b], &
                    [c16b], noma, k8b, [r8b], 'MASSE', &
                    k8b, ibid, xmastr, c16b, k8b, &
                    iret)
        if (iret .eq. 2) then
            call utmess('F', 'SEISME_20', sk=tmas)
        else if (iret .eq. 3) then
            call tbexp2(tmas, 'ENTITE')
            paraki(1) = 'LIEU'
            paraki(2) = 'ENTITE'
            valeki(1) = noma
            valeki(2) = 'TOUT'
            call tbliva(tmas, 2, paraki, [ibid], [r8b], &
                        [c16b], valeki, k8b, [r8b], 'MASSE', &
                        k8b, ibid, xmastr, c16b, k8b, &
                        iret)
            if (iret .ne. 0) then
                call utmess('F', 'SEISME_20', sk=tmas)
            end if
        end if
        calmas = .true.
    else
        calmas = .true.
        xmastr = zero
        xcumul = zero
        do id = 1, 3
            if (zi(jdir+id-1) .eq. 1) then
                do im = 1, nbmode
                    masmod = zr(lval+nbmode*(1+id)+im-1)
                    masuni = zr(lval+nbmode*(7+id)+im-1)
                    if (masuni .ne. rundef) then
                        xmastr = xmastr+masmod
                        xcumul = xcumul+masuni
                    else
                        calmas = .false.
                        xmastr = 1.d0
                        goto 24
                    end if
                end do
                xmastr = xmastr/xcumul
                goto 24
            end if
        end do
24      continue
    end if
!
    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'MASS_EFFE') then
        if (calmas) then
            call utmess('I', 'SEISME_46')
        else
            call utmess('I', 'SEISME_48')
        end if
        mastot(1) = zero
        mastot(2) = zero
        mastot(3) = zero
        cumul(1) = zero
        cumul(2) = zero
        cumul(3) = zero
        do im = 1, nbmode
            ii = 0
            freq = zr(lval+im-1)
            do id = 1, 3
                if (zi(jdir+id-1) .eq. 1) then
                    facpar = zr(lval+nbmode*(4+id)+im-1)
                    masmod = zr(lval+nbmode*(1+id)+im-1)
                    masuni = zr(lval+nbmode*(7+id)+im-1)
                    mastot(id) = mastot(id)+masmod
                    if (masuni .ne. rundef) then
                        xfm = masuni
                    else
                        xfm = masmod/xmastr
                    end if
                    cumul(id) = cumul(id)+xfm
                    if (ii .eq. 0) then
                        ii = 1
                        if (calmas) then
                            call utmess('I', 'SEISME_47', si=im, sk=dir(id), &
                                        nr=5, valr=[freq, facpar, masmod, xfm, cumul(id)])
                        else
                            call utmess('I', 'SEISME_49', si=im, sk=dir(id), &
                                        nr=3, valr=[freq, facpar, masmod])
                        end if
                    else
                        if (calmas) then
                            call utmess('I', 'SEISME_58', sk=dir(id), &
                                        nr=4, valr=[facpar, masmod, xfm, cumul(id)])
                        else
                            call utmess('I', 'SEISME_76', sk=dir(id), &
                                        nr=2, valr=[facpar, masmod])
                        end if
                    end if
                end if
            end do
        end do
        if (calmas) then
            call utmess('I', 'SEISME_50', sr=xmastr)
        end if
        call utmess('I', 'SEISME_51', sr=xmastr)
        do id = 1, 3
            xfm = mastot(id)/xmastr
            if (calmas) then
                if (zi(jdir+id-1) .eq. 1) then
                    call utmess('I', 'SEISME_52', sk=dir(id), &
                                nr=2, valr=[mastot(id), 100.d0*xfm])
                end if
            else
                if (zi(jdir+id-1) .eq. 1) then
                    call utmess('I', 'SEISME_75', sk=dir(id), sr=mastot(id))
                end if
            end if
        end do
    end if
!     --- RECUPERATION DES MODES STATIQUES ---
    call getvtx(' ', 'MULTI_APPUI', nbval=0, nbret=nret1)
    call getfac('DEPL_MULT_APPUI', nret2)
    if ((nret1 .ne. 0) .and. (nret2 .eq. 0)) then
        call utmess('A', 'SEISME_31')
    end if
    call getvid('DEPL_MULT_APPUI', 'MODE_STAT', iocc=1, scal=stat, nbret=ns)
!     --- VERIFICATION - SI GUPTA -> PAS DE MULTI_APPUI ---
    if ((typcmo .eq. 'GUPTA') .and. (nret1 .ne. 0)) then
        call utmess('F', 'SEISME_32')
    end if
!     --- VERIFICATION - SI GUPTA -> F1 < F2 ---
    if ((typcmo .eq. 'GUPTA') .and. (f1gup .ge. f2gup)) then
        call utmess('F', 'SEISME_33')
    end if
!     --- VERIFICATION DES MODES ---
    if (masse .ne. ' ') then
! dans les cas non standards (sous-structuration) on passe sans verif
        call asveri(zk16(jopt), nbopt, meca, psmo, stat, &
                    tronc, monoap, nbsup, zi(jnsu), zk8(jkno), &
                    zi(jdir), zi(jordr), nbmode)
    end if
!     ----- CAS DU MULTI-SUPPORT -----
    if (.not. monoap) then
        call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=nume)
        if (nume .eq. ' ') then
            call utmess('F', 'SEISME_40')
        end if
        call wkvect('&&OP0109.REAC_SUP', 'V V R', nbsup*nbmode*3, jrea)
        call wkvect('&&OP0109.DEPL_SUP', 'V V R', nbsup*3, jdep)
        call wkvect('&&OP0109.TYPE_COM', 'V V I', nbsup*3, jcsu)
        call asmsup(masse, meca, nbmode, neq, nbsup, &
                    zi(jnsu), zk8(jkno), zi(jdir), zr(jrea), zi(jcsu), &
                    nume, zi(jordr))
        call getfac('COMB_DEPL_APPUI', nbfac)
        if (nbfac .ne. 0) call asenap(masse)
    else
        jrea = 1
        jdep = 1
        jcsu = 1
    end if
!
!     --- CALCUL DES REPONSES ---
!
    call ascalc(resu, masse, meca, psmo, stat, &
                nbmode, neq, zi(jordr), zk16(jopt), nbopt, &
                zi(jdir), monoap, muapde, nbsup, zi(jnsu), &
                typcmo, temps, comdir, typcdi, tronc, &
                zr(jamor), zr(jspe), zr(jasy), zk8(jkno), zr(jrea), &
                zr(jdep), zi(jcsu), corfre, f1gup, f2gup)
    if ((.not. monoap) .and. comdir) then
        call utmess('I', 'SEISME_74', sk=typcdi)
    end if
    if (ndepl .ne. 0) call asimpr(nbsup, zi(jcsu), zk8(jkno))
!
!
999 continue
    call titre()
!
!
!     -- CREATION DE L'OBJET .REFD SI NECESSAIRE:
!     -------------------------------------------
!    call refdaj(' ', resu, -1, ' ', 'INIT',&
!                ' ', iret)
!
!
!      DEPLACEMENT: (QN/MN)*DNM, FORCE: (QN/MN*W2)*DNM.
    call jedema()
end subroutine
