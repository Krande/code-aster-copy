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
! aslint: disable=W1501
!
subroutine ccfnrn(option, resuin, resultOut, lisord, nbordr, &
                  resultType)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "blas/dcopy.h"
#include "blas/daxpy.h"
#include "blas/zaxpy.h"
#include "asterc/r8vide.h"
#include "asterfort/asasve.h"
#include "asterfort/ascova.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/gettco.h"
#include "asterfort/infniv.h"
#include "asterfort/ischar.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mcmult.h"
#include "asterfort/memam2.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/rsGetMainPara.h"
#include "asterfort/numecn.h"
#include "asterfort/pcptcc.h"
#include "asterfort/pteddl.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vecinc.h"
#include "asterfort/vecgme.h"
#include "asterfort/vechme.h"
#include "asterfort/vefnme_cplx.h"
#include "asterfort/vefpme.h"
#include "asterfort/vrcins.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "asterfort/medome_once.h"
#include "asterfort/verif_bord.h"
#include "asterfort/ascomb.h"
#include "asterfort/dylach.h"
#include "asterfort/lislec.h"
    integer(kind=8) :: nbordr
    character(len=8) :: resuin, resultOut
    character(len=16) :: option, resultType
    character(len=19) :: lisord
!  CALC_CHAMP - CALCUL DES FORCES NODALES ET DES REACTIONS NODALES
!  -    -                  -      -              -         -
! ----------------------------------------------------------------------
    mpi_int :: mpicou, mpibid
    integer(kind=8) :: jordr, iret, iordr, i, ic, jref, ifm, niv, ibid
    integer(kind=8) :: nuord, nh, jnmo, nbddl, lmat, jvPara, ind, iordk
    integer(kind=8) :: neq, jfo, lonch, lonnew, jfr, jfi, rang, nbproc, nbpas, nbordi
    integer(kind=8) :: lonc2, ltrav, j, inume, jddl, jddr, lacce, p, irelat, jordi
    integer(kind=8) :: cret, jldist, iaux1, k, jcnoch, ideb, ifin, ipas, jvcham, iaux2
    character(len=1) :: stop, ktyp, kbid
    character(len=2) :: codret
    character(len=6) :: nompro
    character(len=8) :: k8bid, kiord, ctyp, nomcmp(3), para, mesh, blanc8
    character(len=16) :: typmo, optio2, motfac, typrep, typmat
    character(len=19) :: listLoad
    character(len=19) :: ligrel, vebid, k19bid, partsd, massgen
    character(len=24) :: numeDofRefe, loadFuncJv, loadNameJv, loadInfoJv, vechmp, vachmp, cnchmp
    character(len=24) :: vecgmp, vacgmp, cncgmp, vefpip, vafpip, cnfpip, vfono(2)
    character(len=24) :: caraElem, cnchmpc
    character(len=24) :: vafono, vreno, vareno, sigma, chdepl, valk(3), chdepk
    character(len=24) :: mateco, materField, vafonr, vafoni, k24b, basemo, numeDof
    character(len=24) :: numeEqua, numeEquaNew
    character(len=24) :: chvive, chacve, masse, chvarc, compor, k24bid, chamno, chamnk
    character(len=24) :: strx, vldist, vcnoch, vcham, lisori
    character(len=24) :: bidon, chacce, kstr
    character(len=8) :: model, modelNew
    aster_logical :: exitim, lstr, lstr2, ldist, dbg_ob, dbgv_ob, ltest, lsdpar, lcpu, lbid
    aster_logical :: lPilo1, lPilo2, l_pmesh, l_complex, noLoads
    real(kind=8) :: etan, time, partps(3), omega2, coef(3), raux
    real(kind=8) :: rctfin, rctdeb, rctfini, rctdebi, freq
    real(kind=8), pointer :: cgmp(:) => null()
    real(kind=8), pointer :: chmp(:) => null()
    real(kind=8), pointer :: fono(:) => null()
    real(kind=8), pointer :: fonor(:) => null()
    real(kind=8), pointer :: fonoi(:) => null()
    real(kind=8), pointer :: fpip(:) => null()
    real(kind=8), pointer :: noch(:) => null()
    real(kind=8), pointer :: reno(:) => null()
    real(kind=8), pointer :: nldepl(:) => null()
    complex(kind=8) :: ci, cun, cmun
    complex(kind=8), pointer :: nochc(:) => null()
    complex(kind=8), pointer :: chmpc(:) => null()
    integer(kind=8), pointer :: v_list_store(:) => null()
    real(kind=8), pointer :: prbid(:) => null()
    complex(kind=8), pointer :: pcbid(:) => null()
    character(len=24), pointer :: refa(:) => null()
    parameter(nompro='CCFNRN')
    data chvarc/'&&CCFNRN.CHVARC'/
    data k24bid/' '/
    data nomcmp/'DX', 'DY', 'DZ'/
    integer(kind=8) :: iret2
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
    call infniv(ifm, niv)
    blanc8 = '        '
! SI PARALLELISME EN TEMPS: INITIALISATION CONTEXTE
    call pcptcc(1, ldist, dbg_ob, dbgv_ob, lcpu, &
                ltest, rang, nbproc, mpicou, nbordr, &
                nbpas, vldist, vcham, lisori, nbordi, &
                lisord, k24b, k24bid, lbid, ibid, &
                ibid, ibid, ibid, ibid, k24b, &
                ibid, ibid, kbid, k24b, prbid, &
                pcbid)
    call jeveuo(vldist, 'L', jldist)
    if (lcpu) call cpu_time(rctdeb)
!
    lonch = -999
    ci = dcmplx(0.D0, 1.D0)
    cun = dcmplx(1.D0, 0.D0)
    cmun = dcmplx(-1.D0, 0.D0)
    l_pmesh = ASTER_FALSE
!
    bidon = '&&'//nompro//'.BIDON'
    listLoad = '&&CCFNRN.LISTLOAD'
!
    if ((option .eq. 'REAC_NODA') .and. &
        ((resultType .eq. 'DYNA_TRANS') .or. (resultType .eq. 'DYNA_HARMO'))) then
        call utmess('A', 'CALCCHAMP_4')
    end if
!
    call jeveuo(lisord, 'L', jordr)
!
! ----ON VERIFIE SI DERRIERE UN CONCEPT MODE_MECA SE TROUVE UN MODE_DYN
    if (resultType(1:9) .eq. 'MODE_MECA') then
        call rsadpa(resuin, 'L', 1, 'TYPE_MODE', 1, &
                    0, sjv=jvPara, styp=k8bid)
        typmo = zk16(jvPara)
    else
        typmo = ' '
    end if
!
! -- LES CHAMPS SONT-ILS DE NATURE COMPLEXE ?
    l_complex = .false.
    iordr = zi(jordr)
    call rsexch(' ', resuin, 'DEPL', iordr, chdepl, &
                iret)
    if (iret .eq. 0) then
        call jelira(chdepl(1:19)//'.VALE', 'TYPE', cval=ktyp)
        l_complex = (ktyp == 'C')
    end if
!
! - Only one list of loads for REAC_NODA
!
    if (option .eq. 'REAC_NODA' .and. &
        (resultType .eq. 'EVOL_ELAS' .or. resultType .eq. 'EVOL_NOLI')) then
        call jeveuo(lisord, 'L', vi=v_list_store)
        call medome_once(resuin, v_list_store, nbordr, list_load_=k19bid)
    end if
!
!
! TRI DES OPTIONS SUIVANT TYPESD
    lmat = 0
    exitim = .false.
    if (resultType .eq. 'EVOL_ELAS' .or. resultType .eq. 'EVOL_NOLI') then
        exitim = .true.
    else if (resultType .eq. 'MODE_MECA' .or. resultType .eq. 'DYNA_TRANS') then
        call jeexin(resuin//'           .REFD', iret)
        if (iret .ne. 0) then
            call dismoi('REF_MASS_PREM', resuin, 'RESU_DYNA', repk=masse, arret='C')
            if (masse .ne. ' ') then
                call mtdscr(masse)
                call jeveuo(masse(1:19)//'.&INT', 'L', lmat)
            end if
        end if
        if (resultType .eq. 'DYNA_TRANS') exitim = .true.
    else if (resultType .eq. 'DYNA_HARMO') then
        call jeexin(resuin//'           .REFD', iret)
        if (iret .ne. 0) then
            call dismoi('REF_MASS_PREM', resuin, 'RESU_DYNA', repk=masse, arret='C')
            if (masse .ne. ' ') then
                call mtdscr(masse)
                call jeveuo(masse(1:19)//'.&INT', 'L', lmat)
            end if
        end if
    end if
    numeDofRefe = " "
    if (resultType .eq. 'MODE_MECA' .or. resultType .eq. 'DYNA_TRANS') then
        call dismoi('NUME_DDL', resuin, 'RESU_DYNA', repk=numeDofRefe, arret='C')
        call jeexin(resuin//'           .REFD', iret)
        if (iret .ne. 0) then
            ! TRAITEMENT DE ELIM_LAGR (traité comme un cas particulier de projection modale)
            call dismoi('BASE_MODALE', resuin, 'RESU_DYNA', repk=basemo)
            call gettco(basemo, typrep)
            if (typrep(1:16) .eq. 'MAILLAGE_SDASTER') then
                call dismoi('REF_MASS_PREM', resuin, 'RESU_DYNA', repk=massgen, arret='C')
                if (massgen .ne. ' ') then
                    call gettco(massgen, typmat)
                    if (typmat .eq. 'MATR_ASSE_ELIM_R') then
                        ! nom de la matrice de masse d'origine
                        call jeveuo(massgen//'.REFA', 'L', vk24=refa)
                        masse = refa(20) (1:19)
                        call mtdscr(masse)
                        call jeveuo(masse(1:19)//'.&INT', 'L', lmat)
                        call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=numeDofRefe)
                    end if
                end if
            end if
        end if
    end if
    caraElem = ' '
    loadNameJv = ' '
    mateco = ' '
    materField = ' '
    nuord = zi(jordr)
    model = ' '
    if (resultType .eq. 'EVOL_THER') then
        call rsGetMainPara("THER", resultOut, nuord, &
                           listLoad, model, materField, mateco, caraElem, &
                           noLoads)
        if (option .eq. "REAC_NODA" .and. noLoads) then
            call utmess('I', 'CALCCHAMP_54')
        end if
    else
        call rsGetMainPara("MECA", resultOut, nuord, &
                           listLoad, model, materField, mateco, caraElem, &
                           noLoads)
        if (option .eq. "REAC_NODA" .and. noLoads) then
            call utmess('I', 'CALCCHAMP_54')
        end if
    end if

    if (model(1:2) .eq. '&&') call utmess('F', 'CALCULEL3_50')
!
! SI PARALLELISME EN TEMPS: ON DEBRANCHE L'EVENTUEL PARALLELISME EN ESPACE
    call pcptcc(2, ldist, dbg_ob, lbid, lbid, &
                lbid, rang, ibid, mpibid, ibid, &
                ibid, k24b, k24b, k24b, ibid, &
                k19bid, model, partsd, lsdpar, ibid, &
                ibid, ibid, ibid, ibid, k24b, &
                ibid, ibid, kbid, k24b, prbid, &
                pcbid)
    if (nbproc .eq. 1 .and. niv > 1) then
        call utmess('I', 'PREPOST_25', sk=option)
    else if (nbproc .gt. 1) then
        call dismoi('NOM_MAILLA', resuin, 'RESULTAT', repk=mesh)
!
        if (ldist) then
            ASSERT(.not. l_pmesh)
            call utmess('I', 'PREPOST_22', si=nbordr, sk=option)
        else if (.not. l_pmesh) then
            if (lsdpar) then
                call utmess('I', 'PREPOST_23', sk=option)
            else
                call utmess('I', 'PREPOST_24', sk=option)
            end if
        end if
    end if
!
    loadFuncJv = listLoad//'.FCHA'
    loadNameJv = listLoad//'.LCHA'
    loadInfoJv = listLoad//'.INFC'

    call exlima(' ', 0, 'V', model, ligrel)
! ON REGARDE S'IL Y A DES ELEMENTS DE STRUCTURE UTILISANT LE CHAMP STRX_ELGA
    strx = ' '
    call dismoi('EXI_STRX', model, 'MODELE', repk=kstr)
    lstr = (kstr(1:3) .eq. 'OUI')
! Y A-T-IL DES ELEMENTS SACHANT CALCULER L'OPTION STRX_ELGA
    call dismoi('EXI_STR2', model, 'MODELE', repk=kstr)
    lstr2 = (kstr(1:3) .eq. 'OUI')
!
!
    if (lcpu) then
        call cpu_time(rctfin)
        write (ifm, *) '< ', rang, 'ccfnrn> Preparation_CPU=', rctfin-rctdeb
    end if
    time = 0.d0
! SI PARALLELISME EN TEMPS: GESTION DE L'INDICE DE DECALAGE
    ipas = 1
    numeDof = ' '
    numeEqua = ' '
    numeEquaNew = ' '
    lonch = -999
    do i = 1, nbordr
        if (lcpu) call cpu_time(rctdebi)
!
! FILTRE POUR EVENTUEL PARALLELISME EN TEMPS
        if (((zi(jldist+i-1) .eq. rang) .and. (ldist)) .or. (.not. ldist)) then
            call jemarq()
!
! SI PARALLELISME EN TEMPS: RECHARGE ADRESSES JEVEUX A CAUSE DU JEMARQ/JEDEMA LOCAL
            if (ldist) then
                call jeveuo(vcham, 'E', jvcham)
                call jeveuo(lisori, 'L', jordi)
                if (ipas .gt. 1) call jeveuo(vcnoch, 'E', jcnoch)
            end if
! SI PARALLELISME EN TEMPS: CALCUL DES INDICES DE DECALAGE
            call pcptcc(4, ldist, dbg_ob, lbid, lbid, &
                        lbid, rang, nbproc, mpibid, ibid, &
                        nbpas, k24b, k24b, k24b, ibid, &
                        k19bid, k24b, k24bid, lbid, i, &
                        ipas, ideb, ifin, irelat, k24b, &
                        ibid, ibid, kbid, k24b, prbid, &
                        pcbid)
            if (lcpu) call cpu_time(rctdeb)
!
            iordr = zi(jordr+i-1)
!
            vechmp = ' '
            vachmp = ' '
            cnchmp = ' '
            vecgmp = ' '
            vacgmp = ' '
            cncgmp = ' '
            vefpip = ' '
            vafpip = ' '
            cnfpip = ' '
            etan = 0.d0
            vfono(1) = ' '
            vfono(2) = ' '
            vafono = ' '
            vafonr = ' '
            vafoni = ' '
            vreno = '&&'//nompro//'           .RELR'
            vareno = '&&'//nompro//'           .RELR'
!
            nh = 0
            if (resultType(1:8) .eq. 'FOURIER_') then
                call rsadpa(resuin, 'L', 1, 'NUME_MODE', iordr, 0, sjv=jnmo)
                nh = zi(jnmo)
            end if
!
            call rsexch(' ', resuin, 'SIEF_ELGA', iordr, sigma, iret)
!
            if (iret .ne. 0) then
                call rsexch(' ', resultOut, 'SIEF_ELGA', iordr, sigma, iret2)
                if (iret2 .ne. 0) then
                    optio2 = 'SIEF_ELGA'
                    if (ldist) then
                        call calcop(optio2, ' ', resuin, resultOut, lisori, &
                                    nbordi, resultType, cret, 'V')
                    else
                        call calcop(optio2, ' ', resuin, resultOut, lisord, &
                                    nbordr, resultType, cret, 'V')
                    end if
                    call rsexch(' ', resultOut, 'SIEF_ELGA', iordr, sigma, iret)
                end if
            end if
!
            if (lstr) then
                call rsexch(' ', resuin, 'STRX_ELGA', iordr, strx, iret)
                if (iret .ne. 0 .and. lstr2) then
                    optio2 = 'STRX_ELGA'
                    if (ldist) then
                        call calcop(optio2, ' ', resuin, resultOut, lisori, &
                                    nbordi, resultType, cret, 'V')
                    else
                        call calcop(optio2, ' ', resuin, resultOut, lisord, &
                                    nbordr, resultType, cret, 'V')
                    end if
                    call rsexch(' ', resultOut, 'STRX_ELGA', iordr, strx, iret)
                end if
            end if
!
            call rsexch(' ', resuin, 'DEPL', iordr, chdepl, iret)
            if (iret .ne. 0) then
                call codent(iordr, 'G', kiord)
                valk(1) = kiord
                valk(2) = option
                call utmess('A', 'PREPOST5_3', nk=2, valk=valk)
                goto 280
            end if
!
!       -- CALCUL D'UN NUME_DDL "MINIMUM" POUR ASASVE :
            if (resultType .eq. 'MODE_MECA' .or. resultType .eq. 'DYNA_TRANS') then
                numeEqua = numeDofRefe(1:14)//'.NUME'
            else
! NUME_DDL QUI PEUT CHANGER AVEC LE PAS DE TEMPS: DONC ON TESTE CET EVENTUEL CHANGEMENT
! SI PARALLELISME EN TEMPS ACTIVE ET PAS DE TEMPS PARALLELISES: NUME
                call numecn(model, chdepl, numeEqua)
                numeEquaNew = numeEqua
                if ((ldist) .and. (ideb .ne. ifin)) then
! SI PARALLELISME EN TEMPS ET NPAS NON ATTEINT: NBPROC CHAM_NOS SIMULTANES
                    do k = ideb, ifin
                        iordk = zi(jordr+k-1)
                        chdepk = ' '
                        call rsexch(' ', resuin, 'DEPL', iordk, chdepk, iret)
                        call numecn(model, chdepk, numeEqua)
                        if (dbg_ob) then
                            write (ifm, *) '< ', rang, &
                                'ccfnrn> numeddl_avant/numddl_apres=', numeEquaNew, numeEqua
                        end if
                        if (numeEquaNew .ne. numeEqua) then
                            call utmess('F', 'PREPOST_16')
                        end if
                    end do
                else
! SI PARALLELISME EN TEMPS et NPAS ATTEINT (RELIQUAT DE PAS DE TEMPS)
! ET SI NON PARALLELISME EN TEMPS
                    if (ldist) then
                        if (dbg_ob) then
                            write (ifm, *) '< ', rang, &
                                'ccfnrn> numeddl_avant/numddl_apres=', numeEqua, numeEquaNew
                        end if
                        if (numeEqua .ne. numeEquaNew) then
                            call utmess('F', 'PREPOST_16')
                        end if
                    end if
                end if
                numeEqua = numeEquaNew
            end if
            numeDof = numeEqua(1:14)
!
            call rsexch(' ', resuin, 'VITE', iordr, chvive, iret)
            if (iret .eq. 0) then
                chvive = '&&'//nompro//'.CHVIT_NUL'
                call copisd('CHAMP_GD', 'V', chdepl, chvive)
                call jelira(chvive(1:19)//'.VALE', 'LONMAX', nbddl)
                call jerazo(chvive(1:19)//'.VALE', nbddl, 1)
            end if
            call rsexch(' ', resuin, 'ACCE', iordr, chacve, iret)
            if (iret .eq. 0) then
                chacve = '&&'//nompro//'.CHACC_NUL'
                call copisd('CHAMP_GD', 'V', chdepl, chacve)
                call jelira(chacve(1:19)//'.VALE', 'LONMAX', nbddl)
                call jerazo(chacve(1:19)//'.VALE', nbddl, 1)
            end if
!
            if (exitim) then
                call rsadpa(resuin, 'L', 1, 'INST', iordr, 0, sjv=jvPara)
                time = zr(jvPara)
            end if
!
            call vrcins(model, materField, caraElem, time, chvarc(1:19), codret)
            call rsexch(' ', resuin, 'COMPORTEMENT', iordr, compor, iret)
!
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step1_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
!
! separation reel imag si dyna_harmo
            call vefnme_cplx(option, 'V', model, mateco, caraElem, &
                             compor, time, time, nh, ligrel, chvarc, sigma, sigma, &
                             strx, chdepl, vfono)
!       --- ASSEMBLAGE DES VECTEURS ELEMENTAIRES ---

            if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                call asasve(vfono(1), numeDof, 'R', vafono)
            else
! creation champ aux noeuds
                call vtcreb(vfono(1), 'V', 'R', nume_ddlz=numeDof)
                call asasve(vfono(1), numeDof, 'R', vafonr)
                call vtcreb(vfono(2), 'V', 'R', nume_ddlz=numeDof)
                call asasve(vfono(2), numeDof, 'R', vafoni)
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step2_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!       --- CREATION DE LA STRUCTURE CHAM_NO ---
            if ((ldist) .and. (ideb .ne. ifin)) then
! SI PARALLELISME EN TEMPS ET NPAS NON ATTEINT: NBPROC CHAM_NOS SIMULTANES
                p = 1
                do k = ideb, ifin
                    iordk = zi(jordr+k-1)
                    call rsexch(' ', resultOut, option, iordk, chamnk, iret)
! CAR LA VARIABLE CHAMNO DOIT ETRE CONNUE POUR L'IORDR COURANT
                    if (iordk .eq. iordr) chamno = chamnk
! ON DOIT PRENDRE LES MEMES DECISIONS QU'EN SEQUENTIEL: NETTOYAGE, MSG...
                    call jeexin(chamnk(1:19)//'.REFE', iret)
                    if (iret .ne. 0) then
                        call codent(iordk, 'G', kiord)
                        valk(1) = option
                        valk(2) = kiord
                        call utmess('I', 'PREPOST5_1', nk=2, valk=valk)
                        call detrsd('CHAM_NO', chamnk(1:19))
                    end if
                    zk24(jvcham+p-1) = ' '
                    zk24(jvcham+p-1) = chamnk
                    if (dbg_ob) write (ifm, *) '< ', rang, 'ccfnrn> p/k/chamnk=', p, k, chamnk
                    p = p+1
                end do
            else
! SINON, 1 SEUL A LA FOIS
! SI PARALLELISME EN TEMPS et NPAS ATTEINT (RELIQUAT DE PAS DE TEMPS)
! OU SI NON PARALLELISME EN TEMPS
                call rsexch(' ', resultOut, option, iordr, chamno, iret)
                call jeexin(chamno(1:19)//'.REFE', iret)
                if (iret .ne. 0) then
                    call codent(iordr, 'G', kiord)
                    valk(1) = option
                    valk(2) = kiord
                    call utmess('I', 'PREPOST5_1', nk=2, valk=valk)
                    call detrsd('CHAM_NO', chamno(1:19))
                end if
            end if
!
! CREATION DES SDS CHAM_NOS SIMPLE OU SIMULTANES
            if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                ktyp = 'R'
                if ((ldist) .and. (ideb .ne. ifin)) then
                    call vtcreb(chamno, 'G', 'R', nume_ddlz=numeDof, nb_equa_outz=neq, &
                                nbz=nbproc, vchamz=vcham)
                else
                    call vtcreb(chamno, 'G', 'R', nume_ddlz=numeDof, nb_equa_outz=neq)
                end if
                call jeveuo(chamno(1:19)//'.VALE', 'E', vr=noch)
            else
                ktyp = 'C'
                if ((ldist) .and. (ideb .ne. ifin)) then
                    call vtcreb(chamno, 'G', 'C', nume_ddlz=numeDof, nb_equa_outz=neq, &
                                nbz=nbproc, vchamz=vcham)
                else
                    call vtcreb(chamno, 'G', 'C', nume_ddlz=numeDof, nb_equa_outz=neq)
                end if
                call jeveuo(chamno(1:19)//'.VALE', 'E', vc=nochc)
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step3_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
!       --- REMPLISSAGE DE L'OBJET .VALE DU CHAM_NO ---
            call jelira(chamno(1:19)//'.VALE', 'LONMAX', lonnew)
! SI PARALLELISME EN TEMPS:
! POUR L'INSTANT, ON SUPPOSE QUE TOUS LES CHAM_NOS SONT DE LONGUEUR IDENTIQUE
! ON TESTE SI C'EST LE CAS SUR LES NBPROCS PAS DE TEMPS CONTIGUES ET SUR LE PAS PRECEDENT
            call pcptcc(6, ldist, dbg_ob, lbid, lbid, &
                        lbid, rang, ibid, mpibid, ibid, &
                        ibid, k24b, k24b, k24b, ibid, &
                        k19bid, k24b, k24bid, lbid, ibid, &
                        ipas, ibid, ibid, ibid, k24b, &
                        lonnew, lonch, kbid, k24b, prbid, &
                        pcbid)
            lonch = lonnew
!
            if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                call jeveuo(vafono, 'L', jfo)
                call jeveuo(zk24(jfo) (1:19)//'.VALE', 'L', vr=fono)
            else
                call jeveuo(vafonr, 'L', jfr)
                call jeveuo(zk24(jfr) (1:19)//'.VALE', 'L', vr=fonor)
                call jeveuo(vafoni, 'L', jfi)
                call jeveuo(zk24(jfi) (1:19)//'.VALE', 'L', vr=fonoi)
                do j = 0, lonch-1
                    nochc(1+j) = dcmplx(fonor(1+j), fonoi(1+j))
                end do
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step4_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
!       --- STOCKAGE DES FORCES NODALES ---
            if (option .eq. 'FORC_NODA') then
                if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, fono, b_incx, noch, b_incy)
                end if
                goto 270
            end if
!
!       --- CALCUL DES FORCES NODALES DE REACTION
!
            if (loadNameJv .ne. ' ') then
                partps(1) = time
!
! --- CHARGES NON PILOTEES (TYPE_CHARGE: 'FIXE_CSTE')
! --- SI LDIST, ON NE VERIFIE QU'AU PREMIER PAS
                if ((.not. ldist) .or. (ldist .and. (ipas .eq. 1))) then
                    if (ligrel(1:8) .ne. model) then
                        stop = 'C'
!               -- on verifie que le ligrel contient bien les mailles de bord
                        call verif_bord(model, ligrel)
                    else
                        stop = 'S'
                    end if
                end if
!
                if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                    call vechme(stop, &
                                model, caraElem, materField, mateco, &
                                loadNameJv, loadInfoJv, &
                                partps, &
                                vechmp, varcCurrZ_=chvarc, &
                                ligrelCalcZ_=ligrel, nharm_=nh)
                    call asasve(vechmp, numeDof, 'R', vachmp)
                    call ascova('D', vachmp, loadFuncJv, 'INST', time, &
                                'R', cnchmp)
!
! --- CHARGES SUIVEUSE (TYPE_CHARGE: 'SUIV')
                    call detrsd('CHAMP_GD', bidon)
                    call vtcreb(bidon, 'G', 'R', nume_ddlz=numeDof, nb_equa_outz=neq)
                    call vecgme('S', &
                                model, caraElem, materField, mateco, compor, &
                                loadNameJv, loadInfoJv, &
                                partps(1), partps(1), &
                                chdepl, bidon, &
                                chvive, chacve, strx, &
                                vecgmp, &
                                ligrelCalcZ_=ligrel)

                    call asasve(vecgmp, numeDof, 'R', vacgmp)
                    call ascova('D', vacgmp, loadFuncJv, 'INST', time, &
                                'R', cncgmp)
                else
                    call rsadpa(resuin, 'L', 1, 'FREQ', iordr, 0, sjv=jvPara)
                    freq = zr(jvPara)
                    if (ligrel(1:8) .ne. model) then
!pour les DYNA_HARMO
!pour l instant je ne fais le calcul de REAC_NODA que sur le model en entier
!(gestion FONC_MULT_C : fastidieuse)
                        call utmess('F', 'PREPOST3_96')
                    end if
                    motfac = 'EXCIT'
                    if ((.not. ldist .and. (i .eq. 1)) .or. (ldist .and. (ipas .eq. 1))) then
                        call lislec(motfac, 'MECANIQUE', 'V', listLoad)
                    else
                        call jedetr(cnchmpc(1:19)//'.REFE')
                        call jedetr(cnchmpc(1:19)//'.DESC')
                        call jedetr(cnchmpc(1:19)//'.VALE')
                    end if
                    vebid = '&&VEBIDON'
                    vechmp = '&&VECHMP'
                    call dylach(model, materField, mateco, caraElem, listLoad, &
                                numeDof, vebid, vechmp, vebid, vebid)
                    para = 'FREQ'
                    cnchmpc = '&&'//nompro//'.CHARGE'
                    call vtcreb(cnchmpc, 'V', 'C', nume_ddlz=numeDof, nb_equa_outz=neq)
                    call ascomb(listLoad, vechmp, 'C', para, freq, &
                                cnchmpc)
                end if
!
!
! --- POUR UN EVOL_NOLI, PRISE EN COMPTE DES FORCES PILOTEES
                if (resultType .eq. 'EVOL_NOLI') then
! - CHARGES PILOTEES (TYPE_CHARGE: 'FIXE_PILO')
                    call vefpme('S', &
                                model, caraElem, materField, mateco, &
                                loadNameJv, loadInfoJv, &
                                partps, &
                                chdepl, bidon, &
                                bidon, &
                                vefpip, ligrel)
                    call asasve(vefpip, numeDof, 'R', vafpip)
                    call ascova('D', vafpip, loadFuncJv, 'INST', time, &
                                'R', cnfpip)
!
! ------------- Loads with continuation method
                    lPilo1 = ischar(listLoad, 'DIRI', 'PILO')
                    lPilo2 = ischar(listLoad, 'NEUM', 'PILO')
                    if (lPilo1 .or. lPilo2) then
                        call rsadpa(resuin, 'L', 1, 'ETA_PILOTAGE', iordr, &
                                    0, sjv=jvPara, istop=0)
                        etan = zr(jvPara)
                        if (etan .eq. r8vide()) then
                            call utmess('F', 'CALCCHAMP_8')
                        end if
                    else
                        etan = 0.d0
                    end if
                end if
!
! --- CALCUL DU CHAMNO DE REACTION PAR DIFFERENCE DES FORCES NODALES
! --- ET DES FORCES EXTERIEURES MECANIQUES NON SUIVEUSES
                if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                    call jeveuo(cnchmp(1:19)//'.VALE', 'L', vr=chmp)
                    call jeveuo(cncgmp(1:19)//'.VALE', 'L', vr=cgmp)
                else
                    call jeveuo(cnchmpc(1:19)//'.VALE', 'L', vc=chmpc)
                end if
                if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                    do j = 0, lonch-1
                        noch(1+j) = fono(1+j)-chmp(1+j)-cgmp(1+j)
                    end do
                else
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zaxpy(b_n, cmun, chmpc, b_incx, nochc, &
                               b_incy)
                end if
                if (resultType .eq. 'EVOL_NOLI') then
                    call jeveuo(cnfpip(1:19)//'.VALE', 'L', vr=fpip)
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, -1.d0*etan, fpip, b_incx, noch, &
                               b_incy)
                end if
            else
!         --- CALCUL DU CHAMNO DE REACTION PAR RECOPIE DE FORC_NODA
                if (resultType .ne. 'DYNA_HARMO' .and. .not. l_complex) then
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, fono, b_incx, noch, b_incy)
                end if
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step5_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
!       --- TRAITEMENT DES MODE_MECA ---
            if (resultType .eq. 'MODE_MECA' .and. typmo(1:8) .eq. 'MODE_DYN' .and. &
                .not. l_complex) then
                call rsadpa(resuin, 'L', 1, 'OMEGA2', iordr, &
                            0, sjv=jvPara, styp=ctyp)
                omega2 = zr(jvPara)
                call jeveuo(chdepl(1:19)//'.VALE', 'L', vr=nldepl)
                call jelira(chdepl(1:19)//'.VALE', 'LONMAX', lonc2)
                call wkvect('&&'//nompro//'.TRAV', 'V V R', lonc2, ltrav)
                if (lmat .eq. 0) call utmess('F', 'PREPOST3_81', sk=option)
                call mrmult('ZERO', lmat, nldepl, zr(ltrav), 1, &
                            .true._1)
                b_n = to_blas_int(lonch)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, -1.d0*omega2, zr(ltrav), b_incx, noch, b_incy)
                call jedetr('&&'//nompro//'.TRAV')
!
!       --- TRAITEMENT DES MODE_STAT ---
            elseif (resultType .eq. 'MODE_MECA' .and. typmo(1:8) .eq. 'MODE_STA' &
                    .and. .not. l_complex) then
                call rsadpa(resuin, 'L', 1, 'TYPE_DEFO', iordr, 0, sjv=jvPara)
                if (zk16(jvPara) (1:9) .eq. 'FORC_IMPO') then
                    call rsadpa(resuin, 'L', 1, 'NUME_DDL', iordr, 0, sjv=jvPara)
                    inume = zi(jvPara)
                    noch(inume) = noch(inume)-1.d0
                else if (zk16(jvPara) (1:9) .eq. 'ACCE_IMPO') then
                    call jelira(chdepl(1:19)//'.VALE', 'LONMAX', lonc2)
                    call rsadpa(resuin, 'L', 1, 'COEF_X', iordr, 0, sjv=jvPara)
                    coef(1) = zr(jvPara)
                    call rsadpa(resuin, 'L', 1, 'COEF_Y', iordr, 0, sjv=jvPara)
                    coef(2) = zr(jvPara)
                    call rsadpa(resuin, 'L', 1, 'COEF_Z', iordr, 0, sjv=jvPara)
                    coef(3) = zr(jvPara)
                    call wkvect('&&'//nompro//'.POSI_DDL', 'V V I', 3*lonc2, jddl)
                    call pteddl('NUME_DDL', numeDof, 3, nomcmp, lonc2, &
                                tabl_equa=zi(jddl))
                    call wkvect('&&'//nompro//'.POSI_DDR', 'V V R', lonc2, jddr)
                    iaux1 = lonc2-1
                    iaux2 = jddl+ind
                    do ic = 1, 3
                        ind = lonc2*(ic-1)
                        raux = coef(ic)
                        do j = 0, iaux1
                            zr(jddr+j) = zr(jddr+j)+zi(jddl+ind+j)*raux
                        end do
                    end do
                    call wkvect('&&'//nompro//'.TRAV', 'V V R', lonc2, ltrav)
                    if (lmat .eq. 0) call utmess('F', 'PREPOST3_81', sk=option)
                    call mrmult('ZERO', lmat, zr(jddr), zr(ltrav), 1, &
                                .true._1)
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, -1.d0, zr(ltrav), b_incx, noch, b_incy)
                    call jedetr('&&'//nompro//'.POSI_DDR')
                    call jedetr('&&'//nompro//'.POSI_DDL')
                    call jedetr('&&'//nompro//'.TRAV')
                end if
!
!       --- TRAITEMENT DE DYNA_TRANS ---
            else if (resultType .eq. 'DYNA_TRANS') then
                call rsexch(' ', resuin, 'ACCE', iordr, chacce, &
                            iret)
                if (iret .eq. 0) then
                    call jeveuo(chacce(1:19)//'.VALE', 'L', lacce)
                    call wkvect('&&'//nompro//'.TRAV', 'V V R', lonch, ltrav)
                    if (lmat .eq. 0) call utmess('F', 'PREPOST3_81', sk=option)
                    call mrmult('ZERO', lmat, zr(lacce), zr(ltrav), 1, .true._1)
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, 1.d0, zr(ltrav), b_incx, noch, b_incy)
                    call jedetr('&&'//nompro//'.TRAV')
                else
                    call utmess('A', 'CALCULEL3_1')
                end if
!
!       --- TRAITEMENT DE DYNA_HARMO ---
            else if (resultType .eq. 'DYNA_HARMO' .or. &
                     (resultType .eq. 'MODE_MECA' .and. typmo(1:8) .eq. 'MODE_STA' &
                      .and. l_complex)) then
                call rsexch(' ', resuin, 'ACCE', iordr, chacce, iret)
                if (iret .eq. 0) then
                    call jeveuo(chacce(1:19)//'.VALE', 'L', lacce)
                    call wkvect('&&'//nompro//'.TRAV', 'V V C', lonch, ltrav)
                    if (lmat .eq. 0) call utmess('F', 'PREPOST3_81', sk=option)
                    call mcmult('ZERO', lmat, zc(lacce), zc(ltrav), 1, .true._1)
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zaxpy(b_n, cun, zc(ltrav), b_incx, nochc, b_incy)
                    call jedetr('&&'//nompro//'.TRAV')
                else
                    call utmess('A', 'CALCULEL3_1')
                end if
!
!       --- TRAITEMENT DE EVOL_NOLI ---
            else if (resultType .eq. 'EVOL_NOLI') then
                call rsexch(' ', resuin, 'ACCE', iordr, chacce, iret)
                if (iret .eq. 0) then
                    optio2 = 'M_GAMMA'
!
!           --- CALCUL DES MATRICES ELEMENTAIRES DE MASSE
                    call memam2(optio2, model, materField, mateco, caraElem, &
                                compor, time, chacce, vreno, 'V', &
                                ligrel)
!
!           --- ASSEMBLAGE DES VECTEURS ELEMENTAIRES ---
                    call asasve(vreno, numeDof, 'R', vareno)
                    call jeveuo(vareno, 'L', jref)
                    call jeveuo(zk24(jref) (1:19)//'.VALE', 'L', vr=reno)
                    b_n = to_blas_int(lonch)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, 1.d0, reno, b_incx, noch, b_incy)
                end if
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step6_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
270         continue
!
! SI PARALLELISME EN TEMPS: ACTIVATION TEST CANONIQUE POUR VERIFIER COM MPI
            if (ltest) then
                if (ktyp .eq. 'R') then
                    noch(:) = (iordr)*1.d0
                else if (ktyp .eq. 'C') then
                    call vecinc(lonch, (iordr)*ci, nochc)
                else
                    ASSERT(.False.)
                end if
            end if
! SI PARALLELISME EN TEMPS:  COM MPI CHAM_NOS.VALE DONT LES NOMS SONT STOCKES DANS VCHAM
            call pcptcc(7, ldist, dbg_ob, lbid, lbid, &
                        lbid, rang, nbproc, mpicou, ibid, &
                        ibid, k24b, vcham, k24b, ibid, &
                        k19bid, k24b, k24bid, lbid, ibid, &
                        ipas, ideb, ifin, irelat, k24b, &
                        ibid, lonch, ktyp, vcnoch, noch, &
                        nochc)
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step7_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
! POST-TRAITEMENTS
            if ((ldist) .and. (ideb .ne. ifin)) then
! SI PARALLELISME EN TEMPS ET NPAS NON ATTEINT: NBPROC CHAM_NOS SIMULTANES
                do k = ideb, ifin
                    iordk = zi(jordr+k-1)
                    call rsnoch(resultOut, option, iordk)
                    modelNew = ' '
                    if (resultType .eq. 'EVOL_THER') then
                        call rsGetMainPara("THER", resultOut, iordk, &
                                           listLoad, modelNew, materField, mateco, caraElem, &
                                           noLoads)
                        if (option .eq. "REAC_NODA" .and. noLoads) then
                            call utmess('I', 'CALCCHAMP_54')
                        end if
                    else
                        call rsGetMainPara("MECA", resultOut, iordk, &
                                           listLoad, modelNew, materField, mateco, caraElem, &
                                           noLoads)
                        if (option .eq. "REAC_NODA" .and. noLoads) then
                            call utmess('I', 'CALCCHAMP_54')
                        end if
                    end if
                    if (dbg_ob) then
                        write (ifm, *) '< ', rang, &
                            'ccfnrn> modele_avant/modele_apres=', model, modelNew
                    end if
                    if (model .ne. modelNew) then
                        call utmess('F', 'PREPOST_1')
                    else
                        model = modelNew
                    end if
                end do
            else
! EN SIMPLE
! SI PARALLELISME EN TEMPS et NPAS ATTEINT (RELIQUAT DE PAS DE TEMPS)
! ET SI NON PARALLELISME EN TEMPS
                call rsnoch(resultOut, option, iordr)
                modelNew = ' '
                if (resultType .eq. 'EVOL_THER') then
                    call rsGetMainPara("THER", resultOut, iordr, &
                                       listLoad, modelNew, materField, mateco, caraElem, &
                                       noLoads)
                    if (option .eq. "REAC_NODA" .and. noLoads) then
                        call utmess('I', 'CALCCHAMP_54')
                    end if
                else
                    call rsGetMainPara("MECA", resultOut, iordr, &
                                       listLoad, modelNew, materField, mateco, caraElem, &
                                       noLoads)
                    if (option .eq. "REAC_NODA" .and. noLoads) then
                        call utmess('I', 'CALCCHAMP_54')
                    end if
                end if
! CAS DE FIGURE DU RELIQUAT DE PAS DE TEMPS
                if (ldist) then
                    if (dbg_ob) then
                        write (ifm, *) '< ', rang, &
                            'ccfnrn> modele_avant/modele_apres=', model, modelNew
                    end if
                    if (model .ne. modelNew) call utmess('F', 'PREPOST_1')
                end if
                model = modelNew
            end if
            if (lcpu) then
                call cpu_time(rctfin)
                write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' step8_CPU=', rctfin-rctdeb
                call cpu_time(rctdeb)
            end if
!
! PARALLELISME EN TEMPS: TEST DE VERIFICATION
            call pcptcc(8, ldist, lbid, dbgv_ob, lbid, &
                        lbid, ibid, ibid, mpibid, ibid, &
                        ibid, k24b, vcham, k24b, ibid, &
                        k19bid, k24b, k24bid, lbid, ibid, &
                        ibid, ideb, ifin, ibid, chamno, &
                        ibid, ibid, kbid, k24b, prbid, &
                        pcbid)
!
            call detrsd('CHAMP_GD', '&&'//nompro//'.SIEF')
            call detrsd('VECT_ELEM', vfono(1) (1:8))
            call detrsd('VECT_ELEM', vfono(2) (1:8))
            call detrsd('VECT_ELEM', vreno(1:8))
            call detrsd('VECT_ELEM', vechmp(1:8))
            call detrsd('VECT_ELEM', vecgmp(1:8))
            call detrsd('VECT_ELEM', vefpip(1:8))
            call detrsd('CHAMP_GD', cnchmp(1:8)//'.ASCOVA')
            call detrsd('CHAMP_GD', cncgmp(1:8)//'.ASCOVA')
            call detrsd('CHAMP_GD', cnfpip(1:8)//'.ASCOVA')
            call jedetr(vachmp(1:8))
            call jedetr(vacgmp(1:8))
            call jedetr(vafpip(1:8))
            call jedetr(vachmp(1:6)//'00.BIDON')
            call jedetr(vacgmp(1:6)//'00.BIDON')
            call jedetr(vafpip(1:6)//'00.BIDON')
            call jedetr(vachmp(1:6)//'00.BIDON     .VALE')
            call jedetr(vacgmp(1:6)//'00.BIDON     .VALE')
            call jedetr(vafpip(1:6)//'00.BIDON     .VALE')
            call jedetr(vachmp(1:6)//'00.BIDON     .DESC')
            call jedetr(vacgmp(1:6)//'00.BIDON     .DESC')
            call jedetr(vafpip(1:6)//'00.BIDON     .DESC')
            call jedetr(vachmp(1:6)//'00.BIDON     .REFE')
            call jedetr(vacgmp(1:6)//'00.BIDON     .REFE')
            call jedetr(vafpip(1:6)//'00.BIDON     .REFE')
            call jedetr(vfono(1) (1:8)//'           .REFE')
            call jedetr(vfono(2) (1:8)//'           .REFE')
            call jedetr(vfono(1) (1:8)//'           .DESC')
            call jedetr(vfono(2) (1:8)//'           .VALE')
            call jedetr(vfono(1) (1:8)//'           .VALE')
            call jedetr(vfono(2) (1:8)//'           .DESC')
            call jedetr(vachmp(1:8)//'.ASCOVA')
            call jedetr(vacgmp(1:8)//'.ASCOVA')
            call jedetr(vafpip(1:8)//'.ASCOVA')
280         continue
! SI PARALLELISME EN TEMPS: GESTION DE L'INDICE DE DECALAGE
            if (ldist) ipas = ipas+1
!
            call jedema()
        end if
! FIN DU IF DISTRIBUTION POUR EVENTUEL PARALLELISME EN TEMPS
!
        if (lcpu) then
            call cpu_time(rctfini)
            write (ifm, *) '< ', rang, 'ccfnrn> Boucle i=', i, ' total_CPU=', rctfini-rctdebi
        end if
    end do
    call detrsd('CHAMP_GD', bidon)
!
! SI PARALLELISME EN TEMPS: NETTOYAGE DU CONTEXTE
    call pcptcc(3, ldist, dbg_ob, lbid, lbid, &
                lbid, rang, ibid, mpibid, ibid, &
                ibid, vldist, vcham, lisori, ibid, &
                k19bid, model, partsd, lsdpar, ibid, &
                ibid, ibid, ibid, ibid, k24b, &
                ibid, ibid, kbid, vcnoch, prbid, &
                pcbid)
!
    call jedema()
end subroutine
