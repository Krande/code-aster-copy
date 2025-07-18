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
subroutine op0146()
    implicit none
!     OPERATEUR PROJ_SPEC_BASE
!     PROJECTION D UN OU PLUSIEURS SPECTRES DE TURBULENCE SUR UNE BASE
!     MODALE PERTURBEE PAR PRISE EN COMPTE DU COUPLAGE FLUIDE STRUCTURE
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/pasfre.h"
#include "asterfort/rebdfr.h"
#include "asterfort/sfifj.h"
#include "asterfort/specep.h"
#include "asterfort/specff.h"
#include "asterfort/spect1.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: lvmoy, n1, n2, nb, nbm
    integer(kind=8) :: dim, ij, nbabs
!-----------------------------------------------------------------------
    integer(kind=8) :: ideb, idisc, ifreq, ik, im1, im2
    integer(kind=8) :: inumo, ipf, ipv, is, ispect, iv, ivali
    integer(kind=8) :: ivite, iz, js, lbaref, lfsvi
    integer(kind=8) :: lfsvk, linds, lnoe, lnozo, lpasf, lspec
    integer(kind=8) :: nbpf, nbspec, nff, nfi, nmodf, nmodi, npoi
    integer(kind=8) :: npv, nspelo, nuzo, nzex
    real(kind=8) :: aire, alonto, freqf, freqi, pas, pui, pui2
    real(kind=8) :: pui2d, pui3d, vmoy, vmoyto, x1, x2, epsi
!-----------------------------------------------------------------------
    parameter(nb=1024)
    integer(kind=8) :: i3, ivitef, lfreq, lnumi, lnumj, lrefe
    real(kind=8) :: val, vitef
    aster_logical :: casint
    character(len=8) :: nomu, option, nomzon, nompro
    character(len=16) :: concep, cmd
    character(len=19) :: base, spectr, typflu
    character(len=24) :: vali, vite, freq, numo
    character(len=24) :: fsvi, fsvk, basref, pvite
    character(len=24) :: valk(3)
    character(len=24) :: chnumi, chnumj, chfreq, chvale
    integer(kind=8) :: mxval
    character(len=16), pointer :: vate(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(nomu, concep, cmd)
!
    call getvis(' ', 'NB_POIN  ', nbval=0, nbret=npoi)
    npoi = abs(npoi)
    nfi = 1
    nff = 1
    if (npoi .ne. 0) then
        call getvis(' ', 'NB_POIN  ', scal=nbpf)
        call getvr8(' ', 'FREQ_INIT', scal=freqi)
        call getvr8(' ', 'FREQ_FIN ', scal=freqf)
    else
        nfi = 0
        nff = 0
    end if
!
! --- 0.VERIFICATIONS AVANT EXECUTION ---
!
    if (npoi .ne. 0) then
        if (freqf .lt. freqi) then
            call utmess('F', 'MODELISA5_70')
        end if
        if (freqf .le. 0.d0 .or. freqi .le. 0.d0) then
            call utmess('F', 'MODELISA5_71')
        end if
        pui2 = log(dble(nbpf))/log(2.d0)
        pui = aint(pui2)
        pui2d = abs(pui2-pui)
        pui3d = abs(1.d0-pui2d)
        if (pui2d .ge. 1.d-06 .and. pui3d .ge. 1.d-06) then
            call utmess('F', 'MODELISA5_72')
        end if
    end if
!
! ----- FIN DES VERIFICATIONS AVANT EXECUTION -----
!
!
! --- 1.1.RECUPERATION DES SPECTRES ET VERIFICATIONS A L'EXECUTION ---
! ---     DE LA COMPATIBILITE DES SPECTRES SI COMBINAISON          ---
!
    call getvid(' ', 'SPEC_TURB', nbval=0, nbret=nbspec)
    nbspec = abs(nbspec)
    call wkvect('&&OP0146.TEMP.NOMS', 'V V K8', nbspec, lspec)
    call getvid(' ', 'SPEC_TURB', nbval=nbspec, vect=zk8(lspec))
!
    call wkvect('&&OP0146.TEMP.INDS', 'V V I', nbspec, linds)
    call wkvect('&&OP0146.TEMP.VMOY', 'V V R', nbspec, lvmoy)
!
    do is = 1, nbspec
        spectr = zk8(lspec+is-1)
        vali = spectr//'.VAIN'
        call jeveuo(vali, 'L', ivali)
        zi(linds+is-1) = zi(ivali)
    end do
!
    if (nbspec .gt. 1) then
        nspelo = 0
        do is = 1, nbspec
            if (zi(linds+is-1) .lt. 10) nspelo = nspelo+1
        end do
        if (nspelo .gt. 0 .and. nspelo .lt. nbspec) then
            call utmess('F', 'MODELISA5_73')
        end if
    end if
!
! --- 2.0 TRAITEMENT SPECIAL POUR SPEC-LONG-COR-5
!
    call jeveuo(zk8(lspec)//'           .VATE', 'L', vk16=vate)
    call wkvect(nomu//'.REFE', 'G V K16', 3, lrefe)
    if (vate(1) .eq. 'SPEC_CORR_CONV_3') then
        zk16(lrefe) = 'DEPL'
        zk16(lrefe+1) = 'TOUT'
    else
        call getvtx(' ', 'OPTION', scal=option)
        zk16(lrefe) = 'SPEC_GENE'
        zk16(lrefe+1) = option
    end if
    zk16(lrefe+2) = 'FREQ'
!
    if (vate(1) (1:14) .eq. 'SPEC_CORR_CONV') then
        call sfifj(nomu)
    else
!
! --- 2.RECUPERATION DES OBJETS DE LA BASE MODALE PERTURBEE ---
!
        call getvid(' ', 'BASE_ELAS_FLUI', scal=base)
!
        vite = base//'.VITE'
        freq = base//'.FREQ'
        numo = base//'.NUMO'
!
        call jeveuo(vite, 'L', ivite)
        call jelira(vite, 'LONUTI', npv)
        call getvr8(' ', 'VITE_FLUI', scal=vitef)
        call getvr8(' ', 'PRECISION', scal=epsi)
!
        ivitef = 0
        do i3 = 1, npv
            val = zr(ivite-1+i3)-vitef
            if (abs(val) .lt. epsi) then
                ivitef = i3
            end if
        end do
        if (ivitef .eq. 0) then
            call utmess('F', 'ALGELINE3_25', sr=vitef)
        end if
!
        call jeveuo(freq, 'L', ifreq)
        call jelira(freq, 'LONUTI', nbm)
        call jeveuo(numo, 'L', inumo)
!
        nbm = nbm/(2*npv)
!
!
! --- 2.1.RECUPERATION DES NOM DES PROFILS DE VITESSE ASSOCIES AUX ---
! ---     SPECTRES DANS LE CAS DES SPECTRES DE TYPE LONGUEUR DE    ---
! ---     TYPE LONGUEUR DE CORRELATION                             ---
!
        if (zi(linds) .lt. 10) then
            call wkvect('&&OP0146.TEMP.NOZ', 'V V K16', nbspec, lnozo)
            do is = 1, nbspec
                spectr = zk8(lspec+is-1)
                call jeveuo(spectr//'.VATE', 'L', vk16=vate)
                zk16(lnozo+is-1) = vate(3)
            end do
!
!
! --- 2.2.VERIFICATION DE L EXISTENCE DES ZONES ASSOCIEES DANS LE   ---
! ---     CONCEPT TYPE_FLUI_STRU ASSOCIE, POUR LES SPECTRES DE TYPE ---
! ---     LONGUEUR DE CORRELATION                                   ---
!
            basref = base//'.REMF'
            call jeveuo(basref, 'L', lbaref)
            typflu = zk8(lbaref)
            fsvi = typflu(1:19)//'.FSVI'
            fsvk = typflu(1:19)//'.FSVK'
            call jeveuo(fsvi, 'L', lfsvi)
            call jeveuo(fsvk, 'L', lfsvk)
            nzex = zi(lfsvi+1)
            pvite = zk8(lfsvk+4)
            pvite = pvite(1:19)//'.VALE'
            call jelira(pvite, 'LONUTI', lnoe)
            lnoe = lnoe/2
!
            do is = 1, nbspec
                nomzon = zk16(lnozo+is-1) (1:8)
                do iz = 1, nzex
                    if (nomzon .eq. zk8(lfsvk+3+iz)) goto 31
                end do
                valk(1) = zk8(lspec+is-1)
                valk(2) = nomzon
                valk(3) = typflu
                call utmess('F', 'MODELISA5_74', nk=3, valk=valk)
31              continue
                valk(1) = zk8(lspec+is-1)
                valk(2) = zk8(lfsvk+iz+3)
                call utmess('I', 'MODELISA5_75', nk=2, valk=valk)
            end do
!
! --- 2.2.ON VERIFIE QUE TOUS LES SPECTRES SONT ASSOCIES A DES ZONES ---
! ---     DIFFERENTES ET SONT DIFFERENTS                             ---
!
            do is = 1, nbspec-1
                do js = is+1, nbspec
                    if (zk8(lspec+is-1) .eq. zk8(lspec+js-1)) then
                        call utmess('F', 'MODELISA5_76')
                    end if
                    nompro = zk16(lnozo+is-1) (1:8)
                    nomzon = zk16(lnozo+js-1) (1:8)
                    if (nompro .eq. nomzon) then
                        valk(1) = zk8(lspec+is-1)
                        valk(2) = zk8(lspec+js-1)
                        valk(3) = nomzon
                        call utmess('F', 'MODELISA5_77', nk=3, valk=valk)
                    end if
                end do
            end do
!
! --- 2.3.CALCUL DES VITESSES MOYENNES DE CHAQUE ZONE D EXCITATION ---
! ---     ET DE LA VITESSE MOYENNE DE L ENSEMBLE DES ZONES         ---
!
!
            vmoyto = 0.d0
            alonto = 0.d0
! ---    BOUCLE SUR LES ZONES D EXCITATION DU FLUIDE
            do nuzo = 1, nzex
                pvite = zk8(lfsvk+3+nuzo)
                pvite = pvite(1:19)//'.VALE'
                call jeveuo(pvite, 'L', ipv)
! ---       RECHERCHE DES EXTREMITES DE LA ZONE 'NUZO'
                do ik = 1, lnoe
                    if (zr(ipv+lnoe+ik-1) .ne. 0.d0) then
                        n1 = ik
                        exit
                    end if
                end do
!
                do ik = lnoe, 1, -1
                    if (zr(ipv+lnoe+ik-1) .ne. 0.d0) then
                        n2 = ik
                        exit
                    end if
                end do
!
                aire = 0.d0
                x1 = zr(ipv+n1-1)
                x2 = zr(ipv+n2-1)
                do ik = n1+1, n2
                    aire = aire+( &
                           zr(ipv+lnoe+ik-1)+zr(ipv+lnoe+ik-2))*(zr(ipv+ik-1)-zr(ipv+i&
                           &k-2) &
                           )/2.d0
                end do
!
                vmoy = aire/(x2-x1)
                zr(lvmoy+nuzo-1) = vmoy
                vmoyto = vmoyto+aire
                alonto = alonto+(x2-x1)
!
! ---   FIN DE BOUCLE SUR LES ZONES D EXCITATION DU FLUIDE
            end do
!
            vmoyto = vmoyto/alonto
!
        end if
!
! --- 3.RECUPERATION DE L'OPTION DE CALCUL
!
        casint = .true.
        call getvtx(' ', 'OPTION', scal=option)
        if (option(1:4) .eq. 'DIAG') casint = .false.
!
! --- 4.DECOUPAGE DE LA BANDE DE FREQUENCE ---
!
!        ---- RECHERCHE DE LA FREQUENCE INITIALE  ET
!                       DE LA FREQUENCE FINALE
!        ---- RECHERCHE DES NUMEROS D ORDRE DES MODES PRIS EN COMPTE
!                       EN FONCTION DE FREQ_INIT ET FREQ_FIN.
!
        call rebdfr(zr(ifreq), nfi, nff, freqi, freqf, &
                    nmodi, nmodf, nbm, npv)
!
!
        dim = (nmodf-nmodi)+1
!
! --- 5.CREATION DE LA TABLE D'INTERSPECTRES ---
!
!
!
!
! --- CREATION D'UN VECTEUR DE TRAVAIL POUR STOCKER LA DISCRETISATION
! --- FREQUENTIELLE
!
        if (npoi .eq. 0) then
            call wkvect('&&OP0146.TEMP.PASF', 'V V R', dim*nb, lpasf)
            nbpf = dim*nb
        else
            call wkvect('&&OP0146.TEMP.PASF', 'V V R', nbpf, lpasf)
            pas = (freqf-freqi)/dble(nbpf-1)
            do ipf = 1, nbpf
                zr(lpasf+ipf-1) = freqi+dble(ipf-1)*pas
            end do
        end if
        call wkvect('&&OP0146.TEMP.DISC', 'V V R', 8*dim, idisc)
!
! --- CREATION DE CHAQUE INTERSPECTRE
!
        mxval = dim*(dim+1)/2
        chnumi = nomu//'.NUMI'
        call wkvect(chnumi, 'G V I', mxval, lnumi)
        chnumj = nomu//'.NUMJ'
        call wkvect(chnumj, 'G V I', mxval, lnumj)
        chvale = nomu//'.VALE'
        call jecrec(chvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                    mxval)
        chfreq = nomu//'.DISC'
        call wkvect(chfreq, 'G V R', nbpf, lfreq)
!
!
        iv = ivitef
!
        if (npoi .eq. 0) then
            call pasfre(zr(idisc), zr(ifreq), zr(lpasf), dim, nbm, &
                        iv, nmodi, freqi, freqf, nb)
        end if
!
        do ipf = 1, nbpf
            zr(lfreq-1+ipf) = zr(lpasf-1+ipf)
        end do
!
        ij = 0
        do im2 = nmodi, nmodf
!
            ideb = im2
            if (casint) ideb = nmodi
!
            do im1 = ideb, im2
                ij = ij+1
!
                zi(lnumi-1+ij) = zi(inumo+im2-1)
                zi(lnumj-1+ij) = zi(inumo+im1-1)
!
                call jecroc(jexnum(chvale, ij))
                if (zi(lnumi-1+ij) .eq. zi(lnumj-1+ij)) then
                    nbabs = nbpf
                else
                    nbabs = 2*nbpf
                end if
                call jeecra(jexnum(chvale, ij), 'LONMAX', nbabs)
                call jeecra(jexnum(chvale, ij), 'LONUTI', nbabs)
!
            end do
        end do
!
!
! --- 6.CALCUL DES INTERSPECTRES D'EXCITATIONS MODALES
! ---   BOUCLE SUR LE NOMBRE DE SPECTRES
!
        do is = 1, nbspec
            ispect = zi(linds+is-1)
            spectr = zk8(lspec+is-1)
            if (ispect .lt. 10) then
                nomzon = zk16(lnozo+is-1) (1:8)
                call spect1(casint, nomu, spectr, ispect, base, &
                            vitef, zi(inumo), nmodi, nmodf, nbm, &
                            nbpf, nomzon, zr(lvmoy+is-1), vmoyto)
            else if (ispect .eq. 11) then
                call specff(casint, nomu, spectr, base, zi(inumo), &
                            nmodi, nmodf, nbm, nbpf)
            else
                call specep(casint, nomu, spectr, base, vitef, &
                            zi(inumo), nmodi, nmodf, nbm, nbpf)
            end if
        end do
! FINSI ALTERNATIVE SPEC-LONG-COR-5
    end if
!
    call titre()
!
    call jedema()
end subroutine
