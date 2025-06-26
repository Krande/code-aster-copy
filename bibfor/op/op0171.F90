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
subroutine op0171()
!
    use listLoad_module
    use loadTherCompute_module
!
    implicit none
!
#include "asterc/etausr.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterf_types.h"
#include "asterfort/addModelLigrel.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/cresol.h"
#include "asterfort/dismoi.h"
#include "asterfort/getMainPara.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnomsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/nttain.h"
#include "asterfort/nttcmv.h"
#include "asterfort/numero.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rssepa.h"
#include "asterfort/sigusr.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE_MO
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: phenom = "THER"
    aster_logical :: matcst, coecst, prem, reasmt, reasvt
    integer(kind=8) :: parcri(9), iifm, jlagp, jvPara
    integer(kind=8) :: k, iret
    integer(kind=8) :: itmaxl, iterl, ifm, niv, num
    integer(kind=8) :: iocc, n1, n2
    integer(kind=8) :: jtemp, jtempm, jtempp, j2nd, lonch, lglap
    real(kind=8) :: tpsthe(6), tpsnp1, testn, testr
    real(kind=8) :: tps1(7), tps2(4), tpex
    real(kind=8) :: parcrr(9), testi, epsr, epsl
    real(kind=8) :: r8aux(1)
    character(len=1) :: ci1, ci2, creas, ce1, ce2
    character(len=8) :: k8bid, answer
    character(len=16) :: k16bid, nomcvg
    character(len=19) :: solver, maprec
    character(len=24) :: listLoad, listLoadResu
    character(len=8) :: model, caraElem, materField
    character(len=24) :: mateco
    character(len=24) :: fieldInResult, vtemp, vtempm, vtempp, vec2nd
    character(len=24) :: result, modelLigrel, tempev, tempin
    character(len=24) :: timeMap, matass, numeDof, noojb
    character(len=24) :: cndirp, cnchci, cnchtp
    character(len=24) :: chlapm, chlapp, cnresi
    character(len=76) :: fmt
    integer(kind=8) :: vali(2)
    real(kind=8) :: valr(2)
    real(kind=8), pointer :: lagpm(:) => null()
    real(kind=8), pointer :: lagpp(:) => null()
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
!
    data maprec/'&&OP0171.MAPREC'/
    data cndirp, cnchtp/2*' '/
    data cnchci, cnresi/2*' '/
    data chlapm, chlapp/'&&OP0171.CLPM', '&&OP0171.CLPP'/
    data vtemp, vec2nd/'&&OP0171.TH', '&&OP0171.2ND'/
    data vtempm, vtempp/'&&OP0171.THM', '&&OP0171.THP'/
    data matass/'&&MTHASS'/
    data fmt/'(76(''*''))'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
    ce1 = ' '
    ce2 = ' '
    solver = '&&OP0171.SOLVER'

! - Get datastructure for results
    call getres(result, k16bid, k8bid)

! - Read main parameters
    call getMainPara(phenom, &
                     model, materField, mateco, caraElem, listLoad)

! - Save list of loads in results datastructure
    call nameListLoad(listLoadResu)
    call copisd('LISTE_CHARGES', 'G', listLoad, listLoadResu)

! - No structural elements for this command
    caraElem = ' '
    call dismoi('EXI_RDM', model, 'MODELE', repk=answer)
    if (answer .eq. 'OUI') then
        call utmess('F', 'THERNONLINE4_2')
    end if

! - Detect non-constant material parameters
    call dismoi('THER_F_INST', materField, 'CHAM_MATER', repk=answer)
    matcst = answer .eq. 'NON'

! - Detect non-constant loads
    call hasFuncLoad(listLoad, coecst)

! - Solver parameters
    call cresol(solver)
!
! --- RECUPERATION DU CRITERE DE CONVERGENCE
!
    nomcvg = 'CONVERGENCE'
    call getfac(nomcvg, iocc)
    if (iocc .eq. 1) then
        call getvr8(nomcvg, 'CRIT_TEMP_RELA', iocc=1, scal=parcrr(4), nbret=parcri(4))
        call getvr8(nomcvg, 'CRIT_ENTH_RELA', iocc=1, scal=parcrr(6), nbret=parcri(6))
!
        call getvis(nomcvg, 'ITER_GLOB_MAXI', iocc=1, scal=parcri(1), nbret=n1)
!
        call getvtx(nomcvg, 'ARRET', iocc=1, scal=answer, nbret=n1)
        parcri(9) = 0
        if (n1 .gt. 0) then
            if (answer .eq. 'NON') then
                parcri(9) = 1
            end if
        end if
    end if
    itmaxl = parcri(1)
    epsr = parcrr(4)
    epsl = parcrr(6)
!
! ======================================================================
!
    timeMap = result(1:8)//'.CHTPS'

! - Generate name of numbering object (numeDof)
    noojb = '12345678.00000.NUME.PRNO'
    call gnomsd(' ', noojb, 10, 14)
    numeDof = noojb(1:14)

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(model, nbLigr, listLigr)

! - Get list of LIGREL from loads
    call getListLoadLigrel(listLoad, nbLigr, listLigr)

! - Numbering
    call numero(numeDof, 'VG', &
                nbLigr, listLigr)
    AS_DEALLOCATE(vk24=listLigr)

!
    call vtcreb(vtemp, 'V', 'R', nume_ddlz=numeDof)
!
!
    call getvid('ETAT_INIT', 'EVOL_THER', iocc=1, scal=tempev, nbret=n1)
    if (n1 .gt. 0) then
        call getvis('ETAT_INIT', 'NUME_ORDRE', iocc=1, scal=num, nbret=n2)
        ASSERT(n2 .gt. 0)
        call rsexch('F', tempev, 'TEMP', num, tempin, iret)
        call vtcopy(tempin, vtemp, iret)
        if (iret .ne. 0) then
            call utmess("F", "FIELD0_6", sk='TEMP')
        end if
    end if
! ======================================================================
!
    call dismoi("NOM_LIGREL", model, "MODELE", repk=modelLigrel)
    r8aux(1) = 0.d0
    call mecact('V', chlapm, 'MODELE', modelLigrel, 'NEUT_R', &
                ncmp=1, nomcmp='X1', sr=r8aux(1))
!
    tpsnp1 = 0.d0
    prem = .true.
!
! ======================================================================
!
    call uttcpu('CPU.OP0171.1', 'INIT', ' ')
    call uttcpu('CPU.OP0171.1', 'DEBUT', ' ')
    call uttcpr('CPU.OP0171.1', 7, tps1)
    tpex = tps1(7)
    call uttcpu('CPU.OP0171.2', 'INIT', ' ')
!
    tpsthe(1) = tpsnp1
    tpsthe(2) = 0.d0
    tpsthe(3) = 0.d0
    tpsthe(4) = 0.d0
    tpsthe(5) = 0.d0
    tpsthe(6) = 0.d0
    write (ifm, fmt)
!
! --- DUPLICATION DES STRUCTURES DE DONNEES ET RECUPERATION D'ADRESSES
!
    call copisd('CHAMP_GD', 'V', vtemp(1:19), vtempm(1:19))
    call copisd('CHAMP_GD', 'V', vtemp(1:19), vtempp(1:19))
    call copisd('CHAMP_GD', 'V', vtemp(1:19), vec2nd(1:19))
    call jeveuo(vtemp(1:19)//'.VALE', 'E', jtemp)
    call jeveuo(vtempm(1:19)//'.VALE', 'E', jtempm)
    call jeveuo(vtempp(1:19)//'.VALE', 'E', jtempp)
    call jeveuo(vec2nd(1:19)//'.VALE', 'E', j2nd)
    call jelira(vec2nd(1:19)//'.VALE', 'LONMAX', lonch)
!
! --- COMPTEUR ET CRITERES D'ARRET
!
    iterl = 0
    testi = 1.d0
    testr = 1.d0
    reasvt = .true.
    reasmt = .true.
    write (ifm, fmt)
    write (ifm, 102)
102 format('*', 1x, 'ITERATION', 1x, '*', 1x, 'CRIT_TEMPER', 1x, '*', 1x, &
           'VALE_TEST_TEMPER', 1x, '*', 1x, 'CRIT_ENTHAL', 1x, '*', 1x, &
           'VALE_TEST_ENTHAL', 1x, '*')
101 format('*', 3x, i4, a1, 7x, 4(1pd11.3, a1, 3x), 3x, '*')
!
    write (ifm, fmt)
!
! ======================================================================
!        ITERATIONS DU PROBLEME DE TRANSPORT EN THERMIQUE N_LINEAIRE
! ======================================================================
!
200 continue
    call uttcpu('CPU.OP0171.2', 'DEBUT', ' ')
!
! --- ACTUALISATION EVENTUELLE DES VECTEURS ET DES MATRICES
!
    call nttcmv(model, mateco, caraElem, listLoad, numeDof, &
                solver, timeMap, tpsthe, tpsnp1, reasvt, &
                reasmt, creas, vtemp, vtempm, vec2nd, &
                matass, maprec, cndirp, cnchci, cnchtp)
    reasmt = .true.
    reasvt = .false.
!
! --- ARRET DES ITERATIONS
!
    if ((testi .gt. epsr .or. testr .gt. epsl) .and. iterl .lt. itmaxl) then
!
! *** ON CONTINUE...
!
        iterl = iterl+1
!
! - ITERATIONS INTERNES
!
        call nttain(model, mateco, caraElem, listLoad, numeDof, &
                    solver, timeMap, epsr, lonch, matass, &
                    maprec, cnchci, cnresi, vtemp, vtempm, &
                    vtempp, vec2nd, chlapm, chlapp, ci1, &
                    ci2, testi)
!
! - ACTUALISATION DU CHAMP ENTHALPIE
!
        if (prem) then
!
            call jelira(chlapp(1:19)//'.CELV', 'LONUTI', lglap)
            call jeveuo(chlapp(1:19)//'.CELV', 'L', jlagp)
            call copisd('CHAMP_GD', 'V', chlapp(1:19), chlapm(1:19))
            prem = .false.
!
        else
!
            call jeveuo(chlapm(1:19)//'.CELV', 'E', vr=lagpm)
            call jeveuo(chlapp(1:19)//'.CELV', 'L', vr=lagpp)
            testr = 0.d0
            testn = 0.d0
            do k = 1, lglap
                testr = testr+(lagpp(k)-lagpm(k))**2
                testn = testn+lagpp(k)**2
                lagpm(k) = lagpp(k)
            end do
            testr = sqrt(testr/testn)
!
        end if
!
! - EVALUATION DE LA CONVERGENCE ET AFFICHAGE
!
        iifm = iunifi('MESSAGE')
        write (iifm, 101) iterl, ce1, epsr, ce2, testi, ce1, epsl, ce2, testr
        call uttcpu('CPU.OP0171.2', 'FIN', ' ')
        call uttcpr('CPU.OP0171.2', 4, tps2)
!
! --- VERIFICATION SI INTERRUPTION DEMANDEE PAR SIGNAL USR1
!
        if (etausr() .eq. 1) then
            call sigusr()
        end if
!
! - Y A-T-IL ASSEZ DE TEMPS POUR REFAIRE UNE ITERATION ?
!
        if (tps2(4) .gt. 0.8d0*tps2(1)-tps2(4)) then
            vali(1) = iterl
            valr(1) = tps2(4)
            valr(2) = tps2(1)
            call utmess('Z', 'DISCRETISATION2_79', si=vali(1), nr=2, valr=valr, &
                        num_except=ASTER_TIMELIMIT_ERROR)
        end if
!
! - ON VA REFAIRE UNE ITERATION
!
        goto 200
!
! *** ON S'ARRETE... (CONVERGENCE OU NOMBRE MAX D'ITERATIONS ATTEINT)
!
    else
!
!
        if ((parcri(9) .eq. 0) .and. (iterl .ge. itmaxl)) then
            write (ifm, fmt)
            call utmess('Z', 'MECANONLINE9_7', num_except=ASTER_CONVERGENCE_ERROR)
        end if
!
    end if
!
! --- FIN DES ITERATIONS
!
    call uttcpu('CPU.OP0171.1', 'FIN', ' ')
    call uttcpr('CPU.OP0171.1', 7, tps1)
    write (ifm, fmt)
    write (ifm, '(A,21X,A,1PE10.2,21X,A)')&
     &                                 '*', 'DUREE:', tps1(7)-tpex, '*'
    write (ifm, fmt)
    write (ifm, '(/)')
!
! ======================================================================
!                   STOCKAGE DU RESULTAT
! ======================================================================
!
    call rscrsd('G', result, 'EVOL_THER', 1)
    call rsexch(' ', result, 'TEMP', 0, fieldInResult, iret)
    call rsadpa(result, 'E', 1, 'INST', 0, 0, sjv=jvPara)
    zr(jvPara) = 0.d0
    call copisd('CHAMP_GD', 'G', vtempp, fieldInResult)
    call rsnoch(result, 'TEMP', 0)
    call rssepa(result, 0, model, materField, caraElem, listLoadResu)
!
    call titre()
!
    call jedema()
end subroutine
