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
subroutine peecin(resu, modele, mate, mateco, cara, nh, nbocc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterc/r8depi.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/chpve2.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim3.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerecu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecact.h"
#include "asterfort/compEnergyKinetic.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/meharm.h"
#include "asterfort/peenca.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nh, nbocc
    character(len=*) :: resu, modele, mate, mateco, cara
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "ENER_CIN"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nd, nr, ni, iret, np, nc, jord, jins, jad, nbordr, iord, numord, iainst, jnmo, ibid
    integer(kind=8) :: ie, nt, nm, ng, nbgrma, ig, jgr, nbma, nume, im, lfreq, nbparr, nbpard
    integer(kind=8) :: nbpaep, iocc, jma, nf, inume, ier, nbMaiT
    parameter(nbpaep=2, nbparr=6, nbpard=4)
    real(kind=8) :: prec, xfreq, varpep(nbpaep), valer(3), inst
    real(kind=8) :: rundf
    character(len=1) :: base
    character(len=2) :: codret
    character(len=8) :: k8b, noma, resul, crit, nommai, nommas, typarr(nbparr), typard(nbpard)
    character(len=8) :: valk(2), nomgd
    character(len=16) :: typres, option, noparr(nbparr), nopard(nbpard), optmas, tabtyp(3)
    character(len=19) :: chelem, knum, kins, field_node, field_elem, ligrel, chvarc, chvref
    character(len=19) :: chdisp, chvite
    character(len=24) :: chmasd, chfreq, typcha, chtime, chgeom, chcara(18)
    character(len=24) :: chtemp, opt, mlggma, chharm, nomgrm, valk2(2)
    aster_logical :: exitim, l_modal
    complex(kind=8) :: c16b
!
    data noparr/'NUME_ORDRE', 'FREQ', 'LIEU', 'ENTITE', 'TOTALE',&
     &     'POUR_CENT'/
    data typarr/'I', 'R', 'K24', 'K8', 'R', 'R'/
    data nopard/'LIEU', 'ENTITE', 'TOTALE', 'POUR_CENT'/
    data typard/'K8', 'K8', 'R', 'R'/
    data tabtyp/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'ELEM#ENER_R'/
    data chvarc, chvref/'&&PEECIN.VARC', '&&PEECIN.VARC_REF'/
!     ------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
!
    base = 'V'
    rundf = r8vide()
    exitim = .false.
    inst = 0.d0
    chdisp = ' '
    chvite = ' '
    chtemp = ' '
    chfreq = ' '
    typres = ' '
    call getvid(' ', 'CHAM_GD', scal=field_node, nbret=nd)
    if (nd .ne. 0) then
        call chpve2(field_node, 3, tabtyp, ier)
        call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
        call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
    end if
    call getvr8(' ', 'FREQ', scal=xfreq, nbret=nf)
    call getvid(' ', 'RESULTAT', scal=resul, nbret=nr)
    call getvr8(' ', 'INST', scal=inst, nbret=ni)

    if (ni .ne. 0) exitim = .true.
    if (nr .ne. 0) then
        call gettco(resul, typres)
        if (typres(1:9) .eq. 'MODE_MECA') then
            noparr(2) = 'FREQ'
        else if (typres(1:9) .eq. 'EVOL_THER' .or. typres(1:9) .eq. 'EVOL_ELAS' .or. &
                 typres(1:9) .eq. 'EVOL_NOLI' .or. typres(1:10) .eq. 'DYNA_TRANS') then
            noparr(2) = 'INST'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    option = 'ENER_CIN'
    call mecham(option, modele, cara, nh, chgeom, &
                chcara, chharm, iret)
    if (iret .ne. 0) goto 90
    noma = chgeom(1:8)
    mlggma = noma//'.GROUPEMA'
!
    call exlim3('ENER_CIN', 'V', modele, ligrel)
!
    knum = '&&PEECIN.NUME_ORDRE'
    kins = '&&PEECIN.INSTANT'
    inume = 1
!
    if (nd .ne. 0) then
        if (nf .eq. 0) then
            xfreq = 1.d0
            call utmess('I', 'UTILITAI3_69')
        else
            call utmess('I', 'UTILITAI3_70')
            xfreq = (r8depi()*xfreq)**2
        end if
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
        call rsutnu(resul, ' ', 0, knum, nbordr, prec, crit, iret)
        if (iret .ne. 0) goto 80
        call jeveuo(knum, 'L', jord)
!        - DANS LE CAS OU CE N'EST PAS UN RESULTAT DE TYPE EVOL_NOLI -
!        --- ON RECUPERE L'OPTION DE CALCUL DE LA MATRICE DE MASSE ---
        if (typres(1:9) .ne. 'EVOL_NOLI') then
            call dismoi('REF_MASS_PREM', resul, 'RESU_DYNA', repk=nommas, arret='C')
            if (nommas .ne. ' ') then
                call dismoi('SUR_OPTION', nommas, 'MATR_ASSE', repk=opt, arret='C', ier=ie)
                if (ie .ne. 0) then
                    call utmess('A', 'UTILITAI3_71')
                else
                    if (opt(1:14) .eq. 'MASS_MECA_DIAG') inume = 0
                end if
            end if
        end if
!        --- ON VERIFIE SI L'UTILISATEUR A DEMANDE L'UTILISATION ---
!        --- D'UNE MATRICE DE MASSE DIAGONALE                    ---
!        --- DANS LA COMMANDE POST_ELEM                          ---
        call getvtx(option(1:9), 'OPTION', iocc=1, scal=optmas, nbret=nt)
        if (optmas(1:14) .eq. 'MASS_MECA_DIAG') then
            inume = 0
            call utmess('I', 'UTILITAI3_72')
        end if
!
        call wkvect(kins, 'V V R', nbordr, jins)
!            CAS D'UN CALCUL MODAL
!        --- ON RECUPERE LES FREQUENCES ---
        call jenonu(jexnom(resul//'           .NOVA', 'FREQ'), iret)
        if (iret .ne. 0) then
            do iord = 1, nbordr
                numord = zi(jord+iord-1)
                call rsadpa(resul, 'L', 1, 'FREQ', numord, &
                            0, sjv=iainst, styp=k8b, istop=0)
                zr(jins+iord-1) = zr(iainst)
            end do
        end if
!            CAS CALCUL TRANSITOIRE
!            RECUPERATION DES INSTANTS
        call jenonu(jexnom(resul//'           .NOVA', 'INST'), iret)
        if (iret .ne. 0) then
            exitim = .true.
            do iord = 1, nbordr
                numord = zi(jord+iord-1)
                call rsadpa(resul, 'L', 1, 'INST', numord, &
                            0, sjv=iainst, styp=k8b, istop=0)
                zr(jins+iord-1) = zr(iainst)
            end do
        end if
        call tbcrsd(resu, 'G')
        call tbajpa(resu, nbparr, noparr, typarr)
    end if
!
    chmasd = '&&PEECIN.MASD'
    call mecact('V', chmasd, 'MAILLA', noma, 'POSI', &
                ncmp=1, nomcmp='POS', si=inume)
!
    do iord = 1, nbordr
        call jemarq()
        call jerecu('V')
        l_modal = ASTER_FALSE
        numord = zi(jord+iord-1)
        inst = zr(jins+iord-1)
        ASSERT(inst .ne. rundf)
        valer(1) = inst
        if (typres .eq. 'FOURIER_ELAS') then
            call rsadpa(resul, 'L', 1, 'NUME_MODE', numord, &
                        0, sjv=jnmo, styp=k8b)
            call meharm(modele, zi(jnmo), chharm)
        end if
        chtime = ' '
        if (exitim) call mechti(noma, inst, rundf, rundf, chtime)
!
        if (nr .ne. 0) then
            call rsexch(' ', resul, 'ECIN_ELEM', numord, field_elem, iret)
            if (iret .gt. 0) then
                if (exitim) then
                    call rsexch(' ', resul, 'VITE', numord, field_node, iret)
                    if (iret .gt. 0) goto 72
                    call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
                    call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
                else
                    l_modal = ASTER_TRUE
                    call rsexch(' ', resul, 'DEPL', numord, field_node, iret)
                    if (iret .gt. 0) goto 72
                    call dismoi('NOM_GD', field_node, 'CHAMP', repk=nomgd)
                    call dismoi('TYPE_SUPERVIS', field_node, 'CHAMP', repk=typcha)
                end if
            else
                call dismoi('NOM_GD', field_elem, 'CHAMP', repk=nomgd)
                call dismoi('TYPE_SUPERVIS', field_elem, 'CHAMP', repk=typcha)
            end if
            if (exitim) then
                xfreq = 1.d0
            else
                call rsadpa(resul, 'L', 1, 'OMEGA2', numord, 0, sjv=lfreq)
                xfreq = zr(lfreq)
            end if
        end if
!
        chfreq = '&&PEECIN.OMEGA2'
        call mecact('V', chfreq, 'MAILLA', noma, 'OME2_R', &
                    ncmp=1, nomcmp='OMEG2', sr=xfreq)
!
        if (typcha(1:7) .eq. 'CHAM_NO') then
            if (nomgd(1:4) .eq. 'DEPL') then
                call vrcins(modele, mate, cara, inst, chvarc, codret)
                call vrcref(modele(1:8), mate(1:8), cara(1:8), chvref(1:19))
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
        chelem = '&&PEECIN.CHAM_ELEM'
        ibid = 0
        if (l_modal) then
            chvite = ' '
            chdisp = field_node
        else
            chvite = field_node
            chdisp = ' '
        end if

        call compEnergyKinetic(modele, ligrel, l_modal, &
                               chdisp, chvite, chfreq, chgeom, mateco, &
                               chcara, chmasd, chvarc, &
                               base, chelem, iret)
30      continue
!
!        --- ON CALCULE L'ENERGIE TOTALE ---
        call peenca(chelem, nbpaep, varpep, 0, [ibid])
!
        do iocc = 1, nbocc
            call getvtx(option(1:9), 'TOUT', iocc=iocc, nbval=0, nbret=nt)
            call getvem(noma, 'MAILLE', option(1:9), 'MAILLE', iocc, &
                        0, k8b, nm)
            call getvem(noma, 'GROUP_MA', option(1:9), 'GROUP_MA', iocc, &
                        0, k8b, ng)
            if (nt .ne. 0) then
                call peenca(chelem, nbpaep, varpep, 0, [ibid])
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
            if (ng .ne. 0) then
                nbgrma = -ng
                call wkvect('&&PEECIN_GROUPM', 'V V K24', nbgrma, jgr)
                call getvem(noma, 'GROUP_MA', option(1:9), 'GROUP_MA', iocc, &
                            nbgrma, zk24(jgr), ng)
                valk2(2) = 'GROUP_MA'
                do ig = 1, nbgrma
                    nomgrm = zk24(jgr+ig-1)
                    call jeexin(jexnom(mlggma, nomgrm), iret)
                    if (iret .eq. 0) then
                        call utmess('A', 'UTILITAI3_46', sk=nomgrm)
                        goto 40
                    end if
                    call jelira(jexnom(mlggma, nomgrm), 'LONUTI', nbma)
                    if (nbma .eq. 0) then
                        call utmess('A', 'UTILITAI3_47', sk=nomgrm)
                        goto 40
                    end if
                    call jeveuo(jexnom(mlggma, nomgrm), 'L', jad)
                    call peenca(chelem, nbpaep, varpep, nbma, zi(jad))
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
40                  continue
                end do
                call jedetr('&&PEECIN_GROUPM')
            end if
            if (nm .ne. 0) then
                nbma = -nm
                call wkvect('&&PEECIN_MAILLE', 'V V K8', nbma, jma)
                call getvem(noma, 'MAILLE', option(1:9), 'MAILLE', iocc, &
                            nbma, zk8(jma), nm)
                valk(2) = 'MAILLE'
                call jelira(noma//'.TYPMAIL', 'LONMAX', nbMaiT)
                do im = 1, nbma
                    nommai = zk8(jma+im-1)
                    nume = char8_to_int(nommai)
                    if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                        call utmess('A', 'UTILITAI3_49', sk=nommai)
                        goto 50
                    end if
                    call peenca(chelem, nbpaep, varpep, 1, [nume])
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
50                  continue
                end do
                call jedetr('&&PEECIN_MAILLE')
            end if
        end do
        call jedetr('&&PEECIN.PAR')
72      continue
        call jedema()
    end do
!
80  continue
    call jedetr(knum)
    call jedetr(kins)
!
90  continue
    call jedema()
end subroutine
