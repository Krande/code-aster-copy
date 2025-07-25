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
subroutine peepot(resu, modele, mate, mateco, cara, &
                  nh, nbocc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterc/r8vide.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/celver.h"
#include "asterfort/chpve2.h"
#include "asterfort/compEnergyPotential.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim3.h"
#include "asterfort/gettco.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/meharm.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/peenca2.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/scalai.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vecint.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nh, nbocc
    character(len=*) :: resu, modele, mate, mateco, cara
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "ENER_POT"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nd, nr, ni, iret, np, nc, jord, jins, jad, nbordr, iord, numord, iainst
    integer(kind=8) :: ire1, ire2, nt, nm, ng, nbgrma, ig, jgr, nbma, nume, im
    integer(kind=8) :: iocc, jma, icheml, ier, nbMaiT, nbparr, nbpard, nbpaep, jnmo, ibid
    parameter(nbpaep=2, nbparr=6, nbpard=4)
    real(kind=8) :: prec, varpep(nbpaep), inst, valer(3), rundf
    character(len=1) :: base
    character(len=2) :: codret
    character(len=8) :: k8b, noma, resul, crit, nommai, typarr(nbparr), typard(nbpard), valk(2)
    character(len=8) :: nomgd
    character(len=16) :: typres, option, optio2, noparr(nbparr), nopard(nbpard)
    character(len=19) :: chelem, knum, kins, ligrel, tabtyp(3), chvarc, chvref, ligrel2
    character(len=19) :: field_node, field_elem
    character(len=24) :: chtime, typcha, chgeom, chcara(18), chtemp, chharm, chdisp
    character(len=24) :: compor, mlggma, nomgrm, valk2(2)
    aster_logical :: exitim, l_temp
    complex(kind=8) :: c16b
!
    mpi_int :: mpicow, mrang, mnbproc, mpicou
    aster_logical :: dbg_ob, lmonit
    integer(kind=8) :: rang, nbproc, k, ntsum, nmsum, nmmax, ngsum, ngmax
    integer(kind=8) :: decalig, decalim, jmntmg, jmigk, jmigi, jmim, niv, ifm, nbgr
    integer(kind=8) :: numpas, numloc, ietdeb, ietrat, ietfin, ietmax
    integer(kind=8) :: longt, icoef, mode, nel, idecgr, j, nbmasum, jnp, ind
    real(kind=8) :: retfin, ztot
    character(len=4) :: docu
    character(len=8) :: k8X, scal
    character(len=24) :: k24X
    character(len=24), pointer :: celk(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
    real(kind=8), pointer :: celv(:) => null()
!
    data noparr/'NUME_ORDRE', 'INST', 'LIEU', 'ENTITE', 'TOTALE',&
     &     'POUR_CENT'/
    data typarr/'I', 'R', 'K24', 'K8', 'R', 'R'/
    data nopard/'LIEU', 'ENTITE', 'TOTALE', 'POUR_CENT'/
    data typard/'K8', 'K8', 'R', 'R'/
    data tabtyp/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'ELEM#ENER_R'/
    data chvarc, chvref/'&&PEEPOT.VARC', '&&PEEPOT.VARC_REF'/
!
!     ------------------------------------------------------------------
    call jemarq()
    c16b = (0.d0, 0.d0)
!
    call infniv(ifm, niv)
! Afin de tracer le temps calcul de chaque etape (lmonit) et pour debugger (dbg_ob)
    lmonit = .false.
    dbg_ob = .false.
    if (lmonit) call system_clock(ietdeb, ietrat, ietmax)
    ntsum = 0
    nmsum = 0
    ngsum = 0
    nmmax = 0
    ngmax = 0
    k8X = 'XXXXXXXX'
    k24X = 'XXXXXXXXXXXXXXXXXXXXXXXX'
!
    base = 'V'
    rundf = r8vide()
    exitim = .false.
    inst = 0.d0
    chtemp = ' '
    chdisp = ' '
    typres = ' '
    call getvid(' ', 'CHAM_GD', scal=field_node, nbret=nd)
    if (nd .ne. 0) then
        call chpve2(field_node, 3, tabtyp, ier)
        call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
        call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
    end if
    call getvid(' ', 'RESULTAT', scal=resul, nbret=nr)
    call getvr8(' ', 'INST', scal=inst, nbret=ni)
    if (ni .ne. 0) exitim = .true.
    if (nr .ne. 0) then
        call gettco(resul, typres)
        if (typres(1:9) .eq. 'MODE_MECA') then
            noparr(2) = 'FREQ'
        else if (typres(1:9) .eq. 'EVOL_THER' .or. typres(1:9) .eq. 'EVOL_ELAS' .or. &
                 typres(1:9) .eq. 'MULT_ELAS' .or. typres(1:9) .eq. 'EVOL_NOLI' .or. &
                 typres(1:10) .eq. 'DYNA_TRANS') then
            noparr(2) = 'INST'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    option = 'ENER_POT'
    call mecham(option, modele, cara, nh, chgeom, &
                chcara, chharm, iret)
    if (iret .ne. 0) goto 90
    noma = chgeom(1:8)
    mlggma = noma//'.GROUPEMA'
!
    call exlim3('ENER_POT', 'V', modele, ligrel)
!
    knum = '&&PEEPOT.NUME_ORDRE'
    kins = '&&PEEPOT.INSTANT'
!
    if (nd .ne. 0) then
        nbordr = 1
        call wkvect(knum, 'V V I', nbordr, jord)
        zi(jord) = 1
        call wkvect(kins, 'V V R', nbordr, jins)
        zr(jins) = inst
        call tbcrsd(resu, 'G')
        call tbajpa(resu, nbpard, nopard, typard)
    else
        call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
        call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
        call rsutnu(resul, ' ', 0, knum, nbordr, &
                    prec, crit, iret)
        if (iret .ne. 0) goto 80
        call jeveuo(knum, 'L', jord)
!        --- ON RECUPERE LES INSTANTS ---
        call wkvect(kins, 'V V R', nbordr, jins)
        call jenonu(jexnom(resul//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            exitim = .true.
            do iord = 1, nbordr
                numord = zi(jord+iord-1)
                call rsadpa(resul, 'L', 1, 'INST', numord, &
                            0, sjv=iainst)
                zr(jins+iord-1) = zr(iainst)
            end do
        else
            call jenonu(jexnom(resul//'           .NOVA', 'FREQ'), iret)
            if (iret .ne. 0) then
                do iord = 1, nbordr
                    numord = zi(jord+iord-1)
                    call rsadpa(resul, 'L', 1, 'FREQ', numord, &
                                0, sjv=iainst)
                    zr(jins+iord-1) = zr(iainst)
                end do
            end if
        end if
        call tbcrsd(resu, 'G')
        call tbajpa(resu, nbparr, noparr, typarr)
    end if
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPELS GETVTX
! AFIN DE NE PAS LE REFAIRE POUR CHAQUE PAS DE TEMPS
!-----------------------------------------------------------------------------
    call wkvect('&&PEEPOT_jmntmg', 'V V I', 3*nbocc, jmntmg)
    if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> creation objet &&PEEPOT_jmntmg'
    do iocc = 1, nbocc
        call getvtx(option(1:9), 'TOUT', iocc=iocc, nbval=0, nbret=nt)
        call getvem(noma, 'MAILLE', option(1:9), 'MAILLE', iocc, &
                    0, k8b, nm)
        call getvem(noma, 'GROUP_MA', option(1:9), 'GROUP_MA', iocc, &
                    0, k8b, ng)
        zi(jmntmg+3*(iocc-1)) = nt
        ntsum = ntsum+abs(nt)
        zi(jmntmg+3*(iocc-1)+1) = nm
        nmmax = max(nmmax, abs(nm))
        nmsum = nmsum+abs(nm)
        zi(jmntmg+3*(iocc-1)+2) = ng
        ngmax = max(ngmax, abs(ng))
        ngsum = ngsum+abs(ng)
    end do
!
    ASSERT((ntsum .ge. 0) .and. (ngsum .ge. 0) .and. (nmsum .ge. 0))
    ASSERT((ngmax .ge. 0) .and. (ngmax .le. ngsum))
    ASSERT((nmmax .ge. 0) .and. (nmmax .le. nmsum))
    if (ngsum .gt. 0) then
        call wkvect('&&PEEPOT_jmigk', 'V V K24', ngsum, jmigk)
        call wkvect('&&PEEPOT_jmigi', 'V V I', ngsum, jmigi)
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> creation objets &&PEEPOT_jmigk/gi'
    end if
    if (nmsum .gt. 0) then
        call wkvect('&&PEEPOT_jmim', 'V V K8', nmsum, jmim)
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> creation objet &&PEEPOT_jmim'
    end if
    decalig = 0
    decalim = 0
    nbmasum = 0
    do iocc = 1, nbocc
        ng = zi(jmntmg+3*(iocc-1)+2)
        if (ng .ne. 0) then
            nbgrma = -ng
            call wkvect('&&PEEPOT_GROUPM', 'V V K24', nbgrma, jgr)
            call getvem(noma, 'GROUP_MA', option(1:9), 'GROUP_MA', iocc, &
                        nbgrma, zk24(jgr), ng)
            do ig = 1, nbgrma
                nomgrm = zk24(jgr+ig-1)
                zk24(jmigk-1+ig+decalig) = nomgrm
                call jeexin(jexnom(mlggma, nomgrm), iret)
                if (iret .eq. 0) then
                    call utmess('A', 'UTILITAI3_46', sk=nomgrm)
                    zk24(jmigk-1+ig+decalig) = k24X
                    zi(jmigi-1+ig+decalig) = -999
                    goto 140
                end if
                call jelira(jexnom(mlggma, nomgrm), 'LONUTI', nbma)
                if (nbma .eq. 0) then
                    call utmess('A', 'UTILITAI3_47', sk=nomgrm)
                    zk24(jmigk-1+ig+decalig) = k24X
                    zi(jmigi-1+ig+decalig) = -999
                    goto 140
                else
                    zi(jmigi-1+ig+decalig) = nbma
                    nbmasum = nbmasum+nbma
                end if
140             continue
            end do
            call jedetr('&&PEEPOT_GROUPM')
            decalig = decalig+nbgrma
! fin if sur nm (groupe de mailles)
        end if
        nm = zi(jmntmg+3*(iocc-1)+1)
        if (nm .ne. 0) then
            nbma = -nm
            call wkvect('&&PEEPOT_MAILLE', 'V V K8', nbma, jma)
            call getvem(noma, 'MAILLE', option(1:9), 'MAILLE', iocc, &
                        nbma, zk8(jma), nm)
            nbmasum = nbmasum+nbma
            call jelira(noma//'.TYPMAIL', 'LONMAX', nbMaiT)
            do im = 1, nbma
                nommai = zk8(jma+im-1)
                nume = char8_to_int(zk8(jma+im-1))
                if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                    call utmess('A', 'UTILITAI3_49', sk=nommai)
                    zk8(jmim-1+im+decalim) = k8X
                    goto 150
                else
                    zk8(jmim-1+im+decalim) = nommai
                end if
150             continue
            end do
            call jedetr('&&PEEPOT_MAILLE')
            decalim = decalim+nbma
! fin if sur nm (liste de mailles)
        end if
! fin boucle sur les iocc
    end do
!
!-----------------------------------------------------------------------------
! PREPARATION DE LA DISTRIBUTION DE TACHES MPI VIA
! FILTRE &PEECA2_vldist (POUR CELUI EN ESPACE-CONNECTIVITE INVERSE DE PEENCA2)
!-----------------------------------------------------------------------------
! Recuperation des donnees MPI pour le //isme en espace de peenca2 (actif par defaut)
    call asmpi_comm('GET_WORLD', mpicow)
    call asmpi_comm('GET', mpicou)
    if (mpicow .ne. mpicou) then
        ASSERT(.False.)
    end if
    call asmpi_info(mpicow, mrang, mnbproc)
    rang = to_aster_int(mrang)
    ASSERT(rang .ge. 0)
    nbproc = to_aster_int(mnbproc)
    ASSERT(nbproc .ge. 1)
!
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPEL PEENCA (step 1): OBJET POUR STOCKER LA CONNECTIVITE INVERSE
! POUR LA KIEME MAILLE DU GROUP_MA: (cf. PEENCA)
! &&PEEPOT_peenca(2*(k-1)+1)=numero du GREL
! &&PEEPOT_peenca(2*(k-1)+2)=indice de l element dans ce GREL
! AFIN DE NE PAS LE REFAIRE POUR CHAQUE MAILLE DE CALCUL DES GROUP_MA OU DES LISTE
! DE MAILLES ET POUR CHAQUE PAS DE TEMPS
!-----------------------------------------------------------------------------
    if (nbmasum .gt. 0) then
        call wkvect('&&PEEPOT_peenca', 'V V I', 2*nbmasum, jnp)
        call vecint(2*nbmasum, 0, zi(jnp))
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> creation objet &&PEEPOT_peenca', nbmasum
    end if
!
    if (lmonit) then
        call system_clock(ietfin)
        retfin = real(ietfin-ietdeb)/real(ietrat)
        write (ifm, *) '< ', rang, 'peepot> temps initialisation globale=', retfin
    end if
    numpas = 0
!
!-----------------------------------------------------------------------------
! BOUCLE PRINCIPALE: PAS DE TEMPS OU MODES OU...
!-----------------------------------------------------------------------------
!
    do iord = 1, nbordr
!
        if (lmonit) call system_clock(ietdeb, ietrat, ietmax)
        numpas = numpas+1
        numloc = iord-(numpas-1)*nbproc
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> iord/numpas/numloc=', iord, numpas, &
            numloc
        call jemarq()
        call jerecu('V')
        icheml = 0
        numord = zi(jord+iord-1)
        inst = zr(jins+iord-1)
        valer(1) = inst
        if (typres .eq. 'FOURIER_ELAS') then
            call rsadpa(resul, 'L', 1, 'NUME_MODE', numord, &
                        0, sjv=jnmo)
            call meharm(modele, zi(jnmo), chharm)
        end if
        chtime = ' '
        if (exitim) call mechti(noma, inst, rundf, rundf, chtime)
!
        if (nr .ne. 0) then
            call rsexch(' ', resul, 'EPOT_ELEM', numord, field_elem, &
                        iret)
            if (iret .gt. 0) then
                call rsexch(' ', resul, 'DEPL', numord, field_node, &
                            ire1)
                if (ire1 .gt. 0) then
                    call rsexch(' ', resul, 'TEMP', numord, field_node, &
                                ire2)
                    if (ire2 .gt. 0) goto 72
                    call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
                    call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
                else
                    call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
                    call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
                end if
            else
                call dismoi('TYPE_SUPERVIS', field_elem, 'CHAMP', repk=typcha)
                call dismoi('NOM_GD', field_elem, 'CHAMP', repk=nomgd)
            end if
        end if
!
        if (typcha(1:7) .eq. 'CHAM_NO') then
            call vrcins(modele, mate, cara, inst, chvarc, &
                        codret)
            call vrcref(modele(1:8), mate(1:8), cara(1:8), chvref(1:19))
            if (nomgd(1:4) .eq. 'DEPL') then
                optio2 = 'EPOT_ELEM'
                l_temp = ASTER_FALSE
            else if (nomgd(1:4) .eq. 'TEMP') then
                optio2 = 'ETHE_ELEM'
                l_temp = ASTER_TRUE
            else
                call utmess('F', 'UTILITAI3_73')
            end if
        else if (typcha(1:9) .eq. 'CHAM_ELEM') then
            if (nomgd(1:4) .eq. 'ENER') then
                chelem = field_elem
                goto 30
            else
                call utmess('F', 'UTILITAI3_73')
            end if
        else
            call utmess('F', 'UTILITAI3_73')
        end if
        icheml = 1
        chelem = '&&PEEPOT.CHAM_ELEM'
        compor = mate(1:8)//'.COMPOR'
        ibid = 0
        if (l_temp) then
            chtemp = field_node
            chdisp = ' '
        else
            chdisp = field_node
            chtemp = ' '
        end if
        if (lmonit) then
            call system_clock(ietfin)
            retfin = real(ietfin-ietdeb)/real(ietrat)
            write (ifm, *) '< ', rang, 'peepot> temps initialisation iord=', iord, retfin
            call system_clock(ietdeb, ietrat, ietmax)
        end if
        call compEnergyPotential(optio2, modele, ligrel, compor, l_temp, &
                                 chdisp, chtemp, chharm, chgeom, mateco, &
                                 chcara, chtime, chvarc, chvref, base, &
                                 chelem, iret)
30      continue
!
!-----------------------------------------------------------------------------
! MUTUALISATION POUR APPEL PEENCA (step 2) AFIN DE NE PAS LE REFAIRE POUR
! CHAQUE MAILLE DE CALCUL DES GROUP_MA OU DES LISTE DE MAILLES.
!
! VERIFICATIONS AU PREMIER PAS DE TEMPS, CALCUL #GREL (NBGR), NOM DU LIGREL (LIGREL2),
! TYPE DE CHAMPS (SCAL), ENERGIE TOTALE (ZTOT)
!-----------------------------------------------------------------------------
        if (iord .eq. 1) then
! on fait ces verifications qu'au premier pas de temps, cela suffit ici
            call celver(chelem, 'NBVARI_CST', 'STOP', ibid)
            call celver(chelem, 'NBSPT_1', 'STOP', ibid)
            call jelira(chelem//'.CELD', 'DOCU', cval=docu)
            if (docu .ne. 'CHML') then
                call utmess('F', 'CALCULEL3_52')
            end if
        end if
        call jeveuo(chelem//'.CELK', 'L', vk24=celk)
        call jeveuo(chelem//'.CELD', 'L', vi=celd)
        call jeveuo(chelem//'.CELV', 'L', vr=celv)
        ligrel2 = celk(1) (1:19)
        nbgr = nbgrel(ligrel2)
        if (iord .eq. 1) then
            scal = scalai(celd(1))
            if (scal(1:1) .ne. 'R') then
                call utmess('F', 'CALCULEL3_74', sk=scal)
            end if
        end if
        ztot = 0.d0
        do j = 1, nbgr
            mode = celd(celd(4+j)+2)
            if (mode .eq. 0) goto 34
            longt = digdel(mode)
            icoef = max(1, celd(4))
            longt = longt*icoef
            nel = nbelem(ligrel2, j)
            idecgr = celd(celd(4+j)+8)
            do k = 1, nel
                ztot = ztot+celv(idecgr+(k-1)*longt)
            end do
34          continue
        end do
        if (lmonit) then
            call system_clock(ietfin)
            retfin = real(ietfin-ietdeb)/real(ietrat)
            write (ifm, *) '< ', rang, 'peepot> temps compEnergyPotential iord=', iord, retfin
        end if
!
! CALCUL ENERGIE TOTALE DEJA DISPONIBLE
        varpep(1) = ztot
        varpep(2) = 100.d0
        decalig = 0
        decalim = 0
        if (numpas .eq. 1) then
            ind = -1
        else
            ind = 1
        end if
!
!-----------------------------------------------------------------------------
! BOUCLE SECONDAIRE: LISTE DE GROUP_MA OU DE MAILLES
!-----------------------------------------------------------------------------
!
        do iocc = 1, nbocc
            if (lmonit) call system_clock(ietdeb, ietrat, ietmax)
! Resultats getvtx deja lus une fois pour toute
            nt = zi(jmntmg+3*(iocc-1))
            nm = zi(jmntmg+3*(iocc-1)+1)
            ng = zi(jmntmg+3*(iocc-1)+2)
!
! Calcul sur 'TOUT'
            if (nt .ne. 0) then
                varpep(1) = ztot
                varpep(2) = 100.d0
                valk(1) = noma
                valk(2) = 'TOUT'
                if (nr .ne. 0) then
                    valer(2) = varpep(1)
                    valer(3) = varpep(2)
                    call tbajli(resu, nbparr, noparr, [numord], valer, &
                                [c16b], valk, 0)
                else
                    call tbajli(resu, nbpard, nopard, [numord], varpep, &
                                [c16b], valk, 0)
                end if
            end if
!
! Calcul sur GROUP_MA
            if (ng .ne. 0) then
                nbgrma = -ng
                valk2(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jmigk-1+ig+decalig)
                    if (nomgrm(1:24) .ne. k24X) then
                        nbma = zi(jmigi-1+ig+decalig)
                        ASSERT(nbma .ne. -999)
                        call jeveuo(jexnom(mlggma, nomgrm), 'L', jad)
                        call peenca2(chelem, nbpaep, varpep, nbma, zi(jad), &
                                     ligrel2, nbgr, ztot, ind, nbproc, &
                                     rang)
                        valk2(1) = nomgrm
                        if (nr .ne. 0) then
                            valer(2) = varpep(1)
                            valer(3) = varpep(2)
                            call tbajli(resu, nbparr, noparr, [numord], valer, &
                                        [c16b], valk2, 0)
                        else
                            call tbajli(resu, nbpard, nopard, [numord], varpep, &
                                        [c16b], valk2, 0)
                        end if
                    end if
                end do
                decalig = decalig+nbgrma
            end if
!
! Calcul sur liste de MA
            if (nm .ne. 0) then
                nbma = -nm
                valk(2) = 'MAILLE'
                do im = 1, nbma
                    nommai = zk8(jmim-1+im+decalim)
                    if (nommai .ne. k8X) then
                        nume = char8_to_int(nommai)
                        call peenca2(chelem, nbpaep, varpep, 1, [nume], &
                                     ligrel2, nbgr, ztot, ind, nbproc, &
                                     rang)
                        valk(1) = nommai
                        if (nr .ne. 0) then
                            valer(2) = varpep(1)
                            valer(3) = varpep(2)
                            call tbajli(resu, nbparr, noparr, [numord], valer, &
                                        [c16b], valk, 0)
                        else
                            call tbajli(resu, nbpard, nopard, [numord], varpep, &
                                        [c16b], valk, 0)
                        end if
                    end if
                end do
                decalim = decalim+nbma
            end if
!
            if (lmonit) then
                call system_clock(ietfin)
                retfin = real(ietfin-ietdeb)/real(ietrat)
                write (ifm, *) '< ', rang, 'peepot> temps peenca/tbajli iord/iocc=', &
                    iord, iocc, retfin
            end if
!
!-----------------------------------------------------------------------------
! FIN DE LA BOUCLE SECONDAIRE
!-----------------------------------------------------------------------------
        end do
!
        call jedetr('&&PEEPOT.PAR')
        if (icheml .ne. 0) call jedetr(chelem)
72      continue
        call jedema()
!
!-----------------------------------------------------------------------------
! FIN DE LA BOUCLE PRINCIPALE
!-----------------------------------------------------------------------------
    end do
!
! Nettoyage des objets de mutualisations et des buffers de com mpi
    call jedetr('&&PEEPOT_jmntmg')
    if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> destruction objet &&PEEPOT_jmntmg'
    if (ngsum .gt. 0) then
        call jedetr('&&PEEPOT_jmigk')
        call jedetr('&&PEEPOT_jmigi')
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> destruction objets &&PEEPOT_jmigk/gi'
    end if
    if (nmsum .gt. 0) then
        call jedetr('&&PEEPOT_jmim')
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> destruction objet &&PEEPOT_jmim'
    end if
    if (nbmasum .gt. 0) then
        call jedetr('&&PEEPOT_peenca')
        if (dbg_ob) write (ifm, *) '< ', rang, 'peepot> creation objet &&PEEPOT_peenca'
    end if
!
80  continue
    call jedetr(knum)
    call jedetr(kins)
!
90  continue
    call jedema()
end subroutine
