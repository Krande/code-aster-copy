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
! aslint: disable=W1504
!
subroutine dladap(result, tinit, lcrea, lamort, neq, &
                  imat, masse, rigid, amort, dep0, &
                  vit0, acc0, fexte, famor, fliai, &
                  nchar, nveca, liad, lifo, modele, &
                  mate, mateco, carele, charge, infoch, &
                  fomult, numedd, nume, numrep, ds_energy, &
                  sd_obsv, mesh)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/etausr.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/dismoi.h"
#include "asterfort/dlarch.h"
#include "asterfort/dlfdyn.h"
#include "asterfort/dlfext.h"
#include "asterfort/dltcrr.h"
#include "asterfort/enerca.h"
#include "asterfort/extdia.h"
#include "asterfort/frqapp.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmarpc.h"
#include "asterfort/recpar.h"
#include "asterfort/sigusr.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "asterfort/nmobse.h"
#include "asterfort/lobs.h"
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL MECANIQUE TRANSITOIRE PAR INTEGRATION DIRECTE
!     AVEC  METHODE EXPLICITE :  DIFFERENCES CENTREES AVEC PAS
!     ADAPTATIF
!
! --------------------------------------------------------------------------------------------------
!
!  IN  : TINIT     : INSTANT DE CALCUL INITIAL
!  IN  : LCREA     : LOGIQUE INDIQUANT SI IL Y A REPRISE
!  IN  : LAMORT    : LOGIQUE INDIQUANT SI IL Y A AMORTISSEMENT
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : IMAT      : TABLEAU D'ADRESSES POUR LES MATRICES
!  IN  : MASSE     : MATRICE DE MASSE
!  IN  : RIGID     : MATRICE DE RIGIDITE
!  IN  : AMORT     : MATRICE D'AMORTISSEMENT
!  IN  : NCHAR     : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!  IN  : NVECA     : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!  IN  : LIAD      : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!  IN  : LIFO      : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!  IN  : MODELE    : NOM DU MODELE
!  IN  : MATE      : NOM DU CHAMP DE MATERIAU
!  IN  : CARELE    : CARACTERISTIQUES DES POUTRES ET COQUES
!  IN  : CHARGE    : LISTE DES CHARGES
!  IN  : INFOCH    : INFO SUR LES CHARGES
!  IN  : FOMULT    : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!  IN  : NUMEDD    : NUME_DDL DE LA MATR_ASSE RIGID
!  IN  : NUME      : NUMERO D'ORDRE DE REPRISE
!  VAR : DEP0      : TABLEAU DES DEPLACEMENTS A L'INSTANT N
!  VAR : VIT0      : TABLEAU DES VITESSES A L'INSTANT N
!  VAR : ACC0      : TABLEAU DES ACCELERATIONS A L'INSTANT N
! IN  NUMREP : NUMERO DE REUSE POUR LA TABLE PARA_CALC
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: neq, imat(*), liad(*), nchar, nveca, nume, numrep
    character(len=8) :: masse, rigid, amort, result
    character(len=24) :: modele, carele, charge, fomult, mate, mateco, numedd
    character(len=24) :: infoch, lifo(*)
    real(kind=8) :: dep0(*), vit0(*), acc0(*), tinit
    real(kind=8) :: fexte(*), famor(*), fliai(*)
    aster_logical :: lamort, lcrea
    type(NL_DS_Energy), intent(inout) :: ds_energy
    integer(kind=8), parameter :: nbtyar = 6
    integer(kind=8) :: ifm, niv, alarm
    integer(kind=8) :: iv1, iv2, ieq, perc, last_prperc, freqpr
    integer(kind=8) :: jdepl
    integer(kind=8) :: jvite, jvit2
    integer(kind=8) :: jacce, jacc2
    integer(kind=8) :: jind1, jind2
    integer(kind=8) :: jvip1, jvip2
    integer(kind=8) :: jvmin, jvmin1, jvmin2
    integer(kind=8) :: istoc
    integer(kind=8) :: nddl
    integer(kind=8) :: iwk0
    integer(kind=8) :: vali(3)
    character(len=4) :: typ1(nbtyar)
    character(len=8) :: k8b
    character(len=8) :: vvar
    character(len=16) :: typres, k16bid, typear(nbtyar)
    character(len=19) :: masse1, rigid1, amort1, k19bid
    character(len=24) :: sop
    character(len=24) :: ndeeq
    real(kind=8) :: tps1(4), tfin
    real(kind=8) :: cmp, cdp, err, dti, dt1, dt2, dtmin, pas1
    real(kind=8) :: temps, temp2, dtmax, t_obs
    real(kind=8) :: dtarch, tarch, tarchi
    real(kind=8) :: epsi, r8val, freq, rtmp
    real(kind=8) :: tjob
    real(kind=8) :: pas2
    real(kind=8) :: valr(3)
    integer(kind=8) :: iwk1, iwk2, iforc1, iexcl
    integer(kind=8) :: nper, nrmax, nr, npas, ipas, iparch, iarchi
    integer(kind=8) :: nnc, nbexcl, nbipas, iveri, nbordr, nbiter
    integer(kind=8) :: nbpasc, ifnobi, ifcibi
    integer(kind=8) :: adeeq
    integer(kind=8) :: ibid
    aster_logical :: ener, l_obsv
    real(kind=8), pointer :: vale(:) => null()
    character(len=19), intent(inout) :: sd_obsv
    character(len=*), intent(in) :: mesh
    integer(kind=8) :: ipas_obs
    integer(kind=8) :: idepmoi, ivitmoi, iaccmoi
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
!====
! 1. LES DONNEES DU CALCUL
!====
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, niv)
!
! 1.2. ==> NOM DES STRUCTURES
!     --- RECUPERATION NOM DE LA COMMANDE ---
!
    call getres(k8b, typres, k16bid)
!
    ndeeq = numedd(1:8)//'      .NUME.DEEQ'
!
!     --- VECTEURS DE TRAVAIL SUR BASE VOLATILE ---
!
    call wkvect('&&DLADAP.FORCE1', 'V V R', neq, iforc1)
    call wkvect('&&DLADAP.F0', 'V V R', neq, iwk0)
    call wkvect('&&DLADAP.F1', 'V V R', neq, iwk1)
    call wkvect('&&DLADAP.F2', 'V V R', neq, iwk2)
    call wkvect('&&DLADAP.DEPL', 'V V R', neq, jdepl)
    call vtcreb('&&DLADAP.DEP2', 'V', 'R', nume_ddlz=numedd)
    call jeveuo('&&DLADAP.DEP2      '//'.VALE', 'E', vr=vale)
    call wkvect('&&DLADAP.VITE', 'V V R', neq, jvite)
    call wkvect('&&DLADAP.VIT2', 'V V R', neq, jvit2)
    call wkvect('&&DLADAP.VIP1', 'V V R', neq, jvip1)
    call wkvect('&&DLADAP.VIP2', 'V V R', neq, jvip2)
    call wkvect('&&DLADAP.ACCE', 'V V R', neq, jacce)
    call wkvect('&&DLADAP.ACC2', 'V V R', neq, jacc2)
    call wkvect('&&DLADAP.VMIN', 'V V R', neq, jvmin)
    call wkvect('&&DLADAP.VMIN1', 'V V R', neq, jvmin1)
    call wkvect('&&DLADAP.VMIN2', 'V V R', neq, jvmin2)
    call wkvect('&&DLADAP.IND1', 'V V I', neq, jind1)
    call wkvect('&&DLADAP.IND2', 'V V I', neq, jind2)
    call jeveuo(ndeeq, 'L', adeeq)
!
    call jeveuo('&&NMCH1P.DEPMOI    .VALE', 'E', idepmoi)
    call jeveuo('&&NMCH1P.VITMOI    .VALE', 'E', ivitmoi)
    call jeveuo('&&NMCH1P.ACCMOI    .VALE', 'E', iaccmoi)
!
!
    epsi = r8prem()
    npas = 0
    iarchi = nume
    ener = ds_energy%l_comp
!
! 1.4. ==> PARAMETRES D'INTEGRATION
!
    call getvr8('INCREMENT', 'INST_FIN', iocc=1, scal=tfin, nbret=ibid)
    call getvr8('INCREMENT', 'PAS', iocc=1, scal=dti, nbret=ibid)
    if (ibid .eq. 0) then
        call utmess('F', 'DYNALINE1_11')
    end if
    if (dti .eq. 0.d0) then
        call utmess('F', 'DYNALINE1_12')
    end if
    dtmax = 0.d0
    call recpar(neq, dti, dtmax, zr(jvmin), vvar, &
                cmp, cdp, dtmin, nper, nrmax)
    nbordr = int((tfin-tinit)/dti)+2
    nbipas = int((tfin-tinit)/dti)+1
!
! 1.5. ==> EXTRACTION DIAGONALE M ET CALCUL VITESSE INITIALE
!
    call dismoi('SUR_OPTION', masse, 'MATR_ASSE', repk=sop, arret='C', &
                ier=ibid)
    if (sop .eq. 'MASS_MECA_DIAG') then
        call extdia(masse, numedd, 2, zr(iwk1))
    else
        call utmess('F', 'DYNALINE1_13')
    end if
!
    do ieq = 1, neq
!
        if (zr(iwk1+ieq-1) .ne. 0.d0) then
!
            zr(iwk1+ieq-1) = 1.0d0/zr(iwk1+ieq-1)
            nddl = zi(adeeq+2*ieq-1)
            do iv1 = 1, neq
                iv2 = ieq+iv1
                if (iv2 .le. neq) then
                    if (zi(adeeq+2*iv2-1) .eq. nddl) then
                        zi(jind2+ieq-1) = iv2
                        goto 152
                    end if
                else
                    goto 152
                end if
            end do
!
152         continue
!
            do iv1 = 1, neq
                iv2 = ieq-iv1
                if (iv2 .gt. 0) then
                    if (zi(adeeq+2*iv2-1) .eq. nddl) then
                        zi(jind1+ieq-1) = iv2
                        goto 154
                    end if
                else
                    goto 154
                end if
            end do
!
154         continue
!
        end if
!
    end do
!
! 1.6. ==> AFFECTATION DES VECTEURS INITIAUX
!
    do ieq = 1, neq
        zr(jdepl+ieq-1) = dep0(ieq)
        zr(jvite+ieq-1) = vit0(ieq)-0.5d0*dti*acc0(ieq)
        zr(jvip1+ieq-1) = vit0(ieq)
        zr(jacce+ieq-1) = acc0(ieq)
        zr(jvmin1+ieq-1) = 1.d-15
        zr(jvmin2+ieq-1) = 1.d-15
    end do
!
! 1.7. ==> --- ARCHIVAGE ---
!
    tarchi = tinit
    nbexcl = 0
    typear(1) = 'DEPL'
    typear(2) = 'VITE'
    typear(3) = 'ACCE'
    if (ener) then
        typear(4) = 'FORC_EXTE'
        typear(5) = 'FORC_AMOR'
        typear(6) = 'FORC_LIAI'
    else
        typear(4) = '         '
        typear(5) = '         '
        typear(6) = '         '
    end if
    call getvis('ARCHIVAGE', 'PAS_ARCH', iocc=1, scal=iparch, nbret=ibid)
    if (ibid .eq. 0) iparch = 1
    dtarch = dti*iparch
    call getvtx('ARCHIVAGE', 'CHAM_EXCLU', iocc=1, nbval=0, nbret=nnc)
    if (nnc .ne. 0) then
        nbexcl = -nnc
        call getvtx('ARCHIVAGE', 'CHAM_EXCLU', iocc=1, nbval=nbexcl, vect=typ1, &
                    nbret=nnc)
    end if
!
    if (nbexcl .eq. nbtyar) then
        call utmess('F', 'ARCHIVAGE_14')
    end if
    do iexcl = 1, nbexcl
        if (typ1(iexcl) .eq. 'DEPL') then
            typear(1) = '    '
        else if (typ1(iexcl) .eq. 'VITE') then
            typear(2) = '    '
        else if (typ1(iexcl) .eq. 'ACCE') then
            typear(3) = '    '
        end if
    end do
!
!====
! 2. CREATION DES CONCEPTS RESULTAT
!====
!
    call dltcrr(result, neq, nbordr, iarchi, 'PREMIER(S)', &
                tinit, lcrea, typres, masse, rigid, &
                amort, dep0, vit0, acc0, fexte, &
                famor, fliai, numedd, nume, nbtyar, &
                typear)
!
!
    call titre()
!
!====
! 3. CALCUL : BOUCLE SUR LES PAS DE TEMPS
!====
!
    ipas = 0
    nbiter = 0
    iveri = 0
    tjob = 0.d0
    nbpasc = 0
    temps = tinit
    tarch = tinit+dtarch
    dt1 = 0.d0
    dt2 = dti
    call uttcpu('CPU.DLADAP', 'INIT', ' ')
!
30  continue
!
    if (ener) then
        do ieq = 1, neq
            fexte(ieq) = fexte(ieq+neq)
        end do
    end if
!
    freqpr = 5
    if (niv .eq. 2) freqpr = 1
    last_prperc = 0
!
    if (temps .lt. tfin) then
! - observation
        t_obs = temps
        ipas_obs = ipas
!
        istoc = 0
        err = 100.d0
        nr = 0
        if (iveri .eq. 0) then
            call uttcpu('CPU.DLADAP', 'DEBUT', ' ')
        else
            if (mod(iveri, nbpasc) .eq. 0) then
                call uttcpu('CPU.DLADAP', 'DEBUT', ' ')
            end if
        end if
!
!        --- DERNIER PAS DE TEMPS ? ---
        if (temps+dt2 .gt. tfin) dt2 = tfin-temps
101     continue
        if (err .gt. 1.d0 .and. nr .lt. nrmax) then
            nbiter = nbiter+1
            pas1 = (dt1+dt2)*0.5d0
            pas2 = dt2*0.5d0
            do ieq = 0, neq-1
!            --- VITESSES AUX INSTANTS INTERMEDIAIRES ------
                zr(jvit2+ieq) = zr(jvite+ieq)+pas1*zr(jacce+ieq)
!            --- DEPLACEMENTS AUX INSTANTS 'TEMPS+DT2' ---------
                vale(ieq+1) = zr(jdepl+ieq)+(dt2*zr(jvit2+ieq))
            end do
! ------------- CALCUL DU SECOND MEMBRE F*
            r8val = temps+dt2
            call dlfext(nveca, nchar, r8val, neq, liad, &
                        lifo, charge, infoch, fomult, modele, &
                        mate, mateco, carele, numedd, zr(iforc1))
!
            if (ener) then
                do ieq = 1, neq
                    fexte(ieq+neq) = zr(iforc1+ieq-1)
                end do
            end if
!
! ------------- FORCE DYNAMIQUE F* = F* - K DEP - C VIT
            call dlfdyn(imat(1), imat(3), lamort, neq, vale, &
                        zr(jvit2), zr(iforc1), zr(iwk0))
!
! ------------- RESOLUTION DE M . A = F ET CALCUL DE VITESSE STOCKEE
!           --- RESOLUTION AVEC FORCE1 COMME SECOND MEMBRE ---
            do ieq = 1, neq
                zr(jacc2+ieq-1) = zr(iwk1+ieq-1)*zr(iforc1+ieq-1)
!           --- VITESSE AUX INSTANTS 'TEMPS+DT2' ---
                zr(jvip2+ieq-1) = zr(jvit2+ieq-1)+pas2*zr(jacc2+ieq-1)
            end do
!
!        --- CALCUL DE VMIN ---
            if (vvar(1:4) .eq. 'MAXI') then
                do ieq = 0, neq-1
                    rtmp = abs(zr(jvite+ieq)*1.d-02)
                    zr(jvmin+ieq) = max(zr(jvmin+ieq), rtmp)
                end do
            else if (vvar(1:4) .eq. 'NORM') then
                do ieq = 0, neq-1
                    if (zr(iwk1+ieq) .ne. 0.d0) then
                        zr(jvmin1+ieq) = 1.d-02*zr(jvit2+(zi(jind1+ieq)-1))
                        zr(jvmin2+ieq) = 1.d-02*zr(jvit2+(zi(jind2+ieq)-1))
                    end if
                    rtmp = 1.d-15
                    zr(jvmin+ieq) = max(rtmp, zr(jvmin1+ieq), zr(jvmin2+ieq))
                end do
            end if
!
!        --- CALCUL DE FREQ. APPARENTE ET ERREUR ---
            call frqapp(dt2, neq, zr(jdepl), vale, zr(jacce), &
                        zr(jacc2), zr(jvmin), freq)
            err = nper*freq*dt2
!
!       --- REDUCTION DU PAS DE TEMPS ---
            if (err .gt. 1.d0) dt2 = dt2/cdp
            if (dt2 .le. dtmin .and. abs(tfin-(temps+dt2)) .gt. epsi) then
                call utmess('F', 'DYNALINE1_17')
            end if
            nr = nr+1
!       LES DEUX LIGNES SUIVANTES SIMULENT LE WHILE - CONTINUE
            goto 101
        else if (err .gt. 1.d0 .and. nr .ge. nrmax) then
            dt2 = dt2*cdp
            valr(1) = temps+dt2
            valr(2) = err
            valr(3) = dt2
            call utmess('A', 'DYNALINE1_16', si=nrmax, nr=3, valr=valr)
        end if
!
        dt1 = dt2
        temp2 = temps+dt2
!
!        --- AUGMENTATION DU PAS SI ERREUR TROP FAIBLE ---
        if (err .lt. 0.75d0) then
            if (npas .eq. 5) then
                dt2 = cmp*dt2
                dt2 = min(dt2, dti)
                npas = 4
            end if
            npas = npas+1
        else
            npas = 0
        end if
        ipas = ipas+1
!
!
!       --- ARCHIVAGE EVENTUEL DANS L'OBJET SOLUTION ---
        if ((temps .le. tarch .and. temp2 .ge. tarch) .or. (temp2 .eq. tfin)) then
            istoc = 0
            alarm = 1
            if ((temp2-tarch) .le. (tarch-temps)) then
                tarchi = temp2
                call dlarch(result, neq, istoc, iarchi, ' ', &
                            alarm, temp2, nbtyar, typear, masse, &
                            vale, zr(jvip2), zr(jacc2), fexte(1+neq), famor(1+neq), &
                            fliai(1+neq))
            else
                tarchi = temps
                call dlarch(result, neq, istoc, iarchi, ' ', &
                            alarm, temps, nbtyar, typear, masse, &
                            zr(jdepl), zr(jvip1), zr(jacce), fexte, famor, &
                            fliai)
            end if
            tarch = tarch+dtarch
        end if
!
        perc = int(100.d0*((temps-tinit)/(tfin-tinit)))
        if (perc .ne. last_prperc) then
            if (mod(perc, freqpr) .eq. 0) then
                call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, ipas], nr=2, &
                            valr=[temps, tarch-dtarch])
                last_prperc = perc
            end if
        end if
!
!
        if (ener) then
            masse1 = masse//'           '
            amort1 = amort//'           '
            rigid1 = rigid//'           '
            call wkvect('FNODABID', 'V V R', 2*neq, ifnobi)
            call wkvect('FCINEBID', 'V V R', 2*neq, ifcibi)
! ON CALCULE LA VITESSE A T N-1
            call enerca(k19bid, zr(jdepl), zr(jvip1), vale, zr(jvip2), &
                        masse1, amort1, rigid1, fexte, famor, &
                        fliai, zr(ifnobi), zr(ifcibi), lamort, .true._1, &
                        .false._1, ds_energy, '&&DLADAP')
            call jedetr('FNODABID')
            call jedetr('FCINEBID')
        end if
!
! ------------- ARCHIVAGE DES PARAMETRES
!
        call nmarpc(ds_energy, numrep, temps)
!
! ------------- TRANSFERT DES NOUVELLES VALEURS DANS LES ANCIENNES
        temps = temp2
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vale, b_incx, zr(jdepl), b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jvit2), b_incx, zr(jvite), b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jvip2), b_incx, zr(jvip1), b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jacc2), b_incx, zr(jacce), b_incy)
!
! - observation
        l_obsv = ASTER_FALSE
        call lobs(sd_obsv, ipas_obs, t_obs, l_obsv)
!
        if (l_obsv) then
! stock depl/vite/acce dans le champ pour Observation
!
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vale, b_incx, zr(idepmoi), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvit2), b_incx, zr(ivitmoi), b_incy)
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jacc2), b_incx, zr(iaccmoi), b_incy)
!
! make observation
            call nmobse(mesh, sd_obsv, t_obs)
        end if
!
! ------------- VERIFICATION DU TEMPS DE CALCUL RESTANT
        if (iveri .eq. 0) then
            call uttcpu('CPU.DLADAP', 'FIN', ' ')
            call uttcpr('CPU.DLADAP', 4, tps1)
            tjob = tps1(1)
            if (tps1(4) .eq. 0.d0) tps1(4) = 1.d-02
            nbpasc = int(1.d-02*(tps1(1)/tps1(4)))+1
        else
            if (mod(iveri, nbpasc) .eq. 0) then
                call uttcpu('CPU.DLADAP', 'FIN', ' ')
                call uttcpr('CPU.DLADAP', 4, tps1)
                if (tps1(1) .le. max(tjob/100.d0, 15.d0)) then
                    goto 999
                end if
                if (tps1(4) .eq. 0.d0) tps1(4) = 1.d-02
                nbpasc = int(1.d-02*(tjob/tps1(4)))+1
            end if
        end if
        iveri = iveri+1
!
        goto 30
    end if
!
999 continue
!
!====
! 4. ARCHIVAGE DU DERNIER INSTANT DE CALCUL POUR LES CHAMPS QUI ONT
!    ETE EXCLUS DE L'ARCHIVAGE AU FIL DES PAS DE TEMPS
!====
!
    if (nbexcl .ne. 0) then
!
        do iexcl = 1, nbexcl
            typear(iexcl) = typ1(iexcl)
        end do
!
        alarm = 0
        call dlarch(result, neq, istoc, iarchi, 'DERNIER(S)', &
                    alarm, temps, nbtyar, typear, masse, &
                    zr(jdepl), zr(jvip1), zr(jacce), fexte(neq+1), famor(neq+1), &
                    fliai(neq+1))
    end if
!
!====
! 5. LA FIN
!====
!
! --- VERIFICATION SI INTERRUPTION DEMANDEE PAR SIGNAL USR1
!
    if (etausr() .eq. 1) then
        call sigusr()
    end if
!
    if (tps1(1) .le. max(tjob/100.d0, 15.d0)) then
        vali(1) = ipas
        vali(2) = iarchi
        vali(3) = nbpasc
        valr(1) = tarchi
        valr(2) = nbpasc*tps1(4)
        valr(3) = tps1(1)
        call utmess('Z', 'DYNAMIQUE_11', ni=3, vali=vali, nr=3, &
                    valr=valr, num_except=ASTER_TIMELIMIT_ERROR)
    end if
!
    vali(1) = ipas
    vali(2) = nbiter
    call utmess('I', 'DYNALINE1_21', ni=2, vali=vali, nr=8, &
                valr=valr)
!
!     --- DESTRUCTION DES OBJETS DE TRAVAIL ---
!
    call jedetc('V', '.CODI', 20)
    call jedetc('V', '.MATE_CODE', 9)
!
    call jedema()
end subroutine
