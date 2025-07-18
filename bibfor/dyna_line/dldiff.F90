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
subroutine dldiff(result, force1, lcrea, lamort, neq, &
                  imat, masse, rigid, amort, dep0, &
                  vit0, acc0, fexte, famor, fliai, &
                  t0, nchar, nveca, liad, lifo, &
                  modele, mate, mateco, carele, charge, infoch, &
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
#include "asterc/r8depi.h"
#include "asterfort/dismoi.h"
#include "asterfort/dlarch.h"
#include "asterfort/dldif0.h"
#include "asterfort/dltcrr.h"
#include "asterfort/dltins.h"
#include "asterfort/dyarch.h"
#include "asterfort/extdia.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sigusr.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "asterfort/nmobse.h"
#include "asterfort/lobs.h"
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL MECANIQUE TRANSITOIRE PAR INTEGRATION DIRECTE
!     AVEC  METHODE EXPLICITE :  DIFFERENCES CENTREES
!
! --------------------------------------------------------------------------------------------------
!
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
!  IN  : T0        : INSTANT DE CALCUL INITIAL
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
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: neq, imat(*), liad(*), nchar, nveca, nume, numrep
    character(len=8) :: masse, rigid, amort
    character(len=24) :: modele, carele, charge, fomult, mate, mateco, numedd
    character(len=24) :: infoch, lifo(*)
    character(len=8) :: result
    character(len=19) :: force1
    real(kind=8) :: dep0(*), vit0(*), acc0(*), t0
    real(kind=8) :: fexte(*), famor(*), fliai(*)
    aster_logical :: lamort, lcrea
    type(NL_DS_Energy), intent(inout) :: ds_energy
    integer(kind=8), parameter :: nbtyar = 6
    integer(kind=8) :: iwk0, iwk1, iwk2
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ieq, iexcl, perc, freqpr, last_prperc
    integer(kind=8) :: ivite1, ivite2, iacce1, iarchi
    integer(kind=8) :: ibid
    integer(kind=8) :: alarm, archiv
    integer(kind=8) :: ipepa, igrpa
    integer(kind=8) :: ipas, istop, istoc, jstoc
    integer(kind=8) :: jnbpa, jbint, jlpas
    integer(kind=8) :: npatot, nbgrpa, nbptpa
    integer(kind=8) :: nbexcl, nbordr
    character(len=4) :: typ1(nbtyar)
    character(len=8) :: nomres
    character(len=16) :: typres, nomcmd, typear(nbtyar)
    character(len=19) :: lisarc
    character(len=24) :: lisins, lispas, libint, linbpa
    character(len=24) :: sop
    real(kind=8) :: tps1(4), tps2(4), lastarch
    real(kind=8) :: dt, dtm, dtmax, temps, dt1, tf
    real(kind=8) :: omeg, deuxpi
    real(kind=8) :: r8bid
    integer(kind=8) :: vali(2)
    real(kind=8) :: valr(2)
    aster_logical :: ener, l_obsv
    real(kind=8), pointer :: vale(:) => null()
    character(len=19), intent(inout) :: sd_obsv
    character(len=*), intent(in) :: mesh
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
    call getres(nomres, typres, nomcmd)
!
! 1.3. ==> VECTEURS DE TRAVAIL SUR BASE VOLATILE ---
!
    call wkvect('&&DLDIFF.F0', 'V V R', neq, iwk0)
    call wkvect('&&DLDIFF.F1', 'V V R', neq, iwk1)
    call wkvect('&&DLDIFF.F2', 'V V R', neq, iwk2)
    call vtcreb('&&DLDIFF.DEPL1', 'V', 'R', nume_ddlz=numedd)
    call jeveuo('&&DLDIFF.DEPL1     '//'.VALE', 'E', vr=vale)
    call wkvect('&&DLDIFF.VITE1', 'V V R', neq, ivite1)
    call wkvect('&&DLDIFF.VITE2', 'V V R', neq, ivite2)
    call wkvect('&&DLDIFF.ACCE1', 'V V R', neq, iacce1)
!
    deuxpi = r8depi()
    iarchi = nume
    ener = ds_energy%l_comp
!
! 1.4. ==> PARAMETRES D'INTEGRATION
!
    call dltins(nbgrpa, lispas, libint, linbpa, npatot, &
                t0, lisins)
    call jeveuo(lispas, 'L', jlpas)
    call jeveuo(libint, 'L', jbint)
    call jeveuo(linbpa, 'L', jnbpa)
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
    dt1 = zr(jlpas)
    r8bid = dt1/2.d0
    do ieq = 1, neq
        if (zr(iwk1+ieq-1) .ne. 0.d0) then
            zr(iwk1+ieq-1) = 1.0d0/zr(iwk1+ieq-1)
        end if
        vit0(1+ieq) = vit0(1+ieq)-r8bid*acc0(1+ieq)
    end do
!
! 1.6. ==> --- ARCHIVAGE ---
!
    lisarc = '&&DLDIFF.ARCHIVAGE'
    call dyarch(npatot, lisins, lisarc, nbordr, 1, &
                nbexcl, typ1)
    call jeveuo(lisarc, 'E', jstoc)
!
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
! 2. CREATION DU CONCEPT RESULTAT
!====
!
    t0 = zr(jbint)
    call dltcrr(result, neq, nbordr, iarchi, ' ', &
                t0, lcrea, typres, masse, rigid, &
                amort, dep0, vit0, acc0, fexte, &
                famor, fliai, numedd, nume, nbtyar, &
                typear)
!
!
    call titre()
!
!
!
!====
! 3. CALCUL
!====
!
! 3.1. ==> BOUCLE SUR LES GROUPES DE PAS DE TEMPS
    istop = 0
    ipas = 0
!
    call uttcpu('CPU.DLDIFF.1', 'INIT', ' ')
    call uttcpu('CPU.DLDIFF.2', 'INIT', ' ')
!
    do igrpa = 1, nbgrpa
!
! 3.1.1. ==> PREALABLES
!
        call uttcpu('CPU.DLDIFF.1', 'DEBUT', ' ')
        dt = zr(jlpas-1+igrpa)
        nbptpa = zi(jnbpa-1+igrpa)
        t0 = zr(jbint-1+igrpa)
        tf = zr(jbint+igrpa)
!
! 3.1.2. ==> VERIFICATION DU PAS DE TEMPS
!
        call extdia(rigid, numedd, 2, zr(iwk2))
        ibid = 0
        dtmax = dt
        do ieq = 1, neq
            if (zr(iwk1+ieq-1) .ne. 0.d0) then
                omeg = sqrt(zr(iwk2+ieq-1)*zr(iwk1+ieq-1))
                dtm = 5.d-02*deuxpi/omeg
                if (dtmax .gt. dtm) then
                    dtmax = dtm
                    ibid = 1
                end if
            end if
        end do
!
        if (ibid .eq. 1) then
            vali(1) = nint((tf-t0)/dtmax)
            vali(2) = igrpa
            valr(1) = dt
            valr(2) = dtmax
            call utmess('F', 'DYNAMIQUE_12', ni=2, vali=vali, nr=2, &
                        valr=valr)
        end if
! ==> FIN DE VERIFICATION
!
!
        freqpr = 5
        if (niv .eq. 2) freqpr = 1
        last_prperc = 0
! 3.1.3. ==> BOUCLE SUR LES NBPTPA "PETITS" PAS DE TEMPS
!
        do ipepa = 1, nbptpa
            ipas = ipas+1
            if (ipas .gt. npatot) goto 99
            istoc = 0
            temps = t0+dt*ipepa
            call uttcpu('CPU.DLDIFF.2', 'DEBUT', ' ')
            archiv = zi(jstoc+ipas-1)
!
            call dldif0(result, force1, neq, istoc, iarchi, &
                        lamort, imat, masse, rigid, amort, &
                        dep0, vit0, acc0, vale, zr(ivite1), &
                        zr(iacce1), zr(ivite2), fexte(1), famor(1), fliai(1), &
                        nchar, nveca, liad, lifo, modele, &
                        ener, mate, mateco, carele, charge, &
                        infoch, fomult, numedd, dt, temps, &
                        zr(iwk0), zr(iwk1), archiv, nbtyar, typear, &
                        numrep, ds_energy)
!
            if (archiv .eq. 1) lastarch = temps
            perc = int(100.d0*(real(ipas)/real(npatot)))
            if (perc .ne. last_prperc) then
                if (mod(perc, freqpr) .eq. 0) then
                    call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, ipas], nr=2, &
                                valr=[temps, lastarch])
                    last_prperc = perc
                end if
            end if

            ! - SI OBSERVATION

            l_obsv = ASTER_FALSE
            call lobs(sd_obsv, ipas, temps, l_obsv)
            if (l_obsv) then
                call nmobse(mesh, sd_obsv, temps)
            end if
!
!
! 3.5. ==> VERIFICATION DU TEMPS DE CALCUL RESTANT
!
            call uttcpu('CPU.DLDIFF.2', 'FIN', ' ')
            call uttcpr('CPU.DLDIFF.2', 4, tps2)
            if (tps2(1) .lt. 5.d0 .or. tps2(4) .gt. tps2(1)) then
                if (ipepa .ne. npatot) then
                    istop = 1
                    vali(1) = igrpa
                    vali(2) = ipepa
                    valr(1) = tps2(4)
                    valr(2) = tps2(1)
                    goto 99
                end if
            end if
!
! ---------- FIN DE LA BOUCLE SUR LES NBPTPA "PETITS" PAS DE TEMPS
        end do
        call uttcpu('CPU.DLDIFF.1', 'FIN', ' ')
        call uttcpr('CPU.DLDIFF.1', 4, tps1)
        if (tps1(1) .lt. 5.d0 .and. igrpa .ne. nbgrpa) then
            istop = 1
            vali(1) = igrpa
            vali(2) = ipepa
            valr(1) = tps1(4)
            valr(2) = tps1(1)
            goto 99
        end if
! ------- FIN BOUCLE SUR LES GROUPES DE PAS DE TEMPS
    end do
!
99  continue
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
!
        call dlarch(result, neq, istoc, iarchi, ' ', &
                    alarm, temps, nbexcl, typear, masse, &
                    vale, zr(ivite1), zr(iacce1), fexte(1+neq), famor(1+neq), &
                    fliai(1+neq))
!
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
    if (istop .eq. 1) then
        call utmess('Z', 'DYNAMIQUE_10', ni=2, vali=vali, nr=2, &
                    valr=valr, num_except=ASTER_TIMELIMIT_ERROR)
    end if
!
!     --- DESTRUCTION DES OBJETS DE TRAVAIL ---
!
    call jedetc('V', '.CODI', 20)
    call jedetc('V', '.MATE_CODE', 9)
!
    call jedema()
!
end subroutine
