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

subroutine pofape()
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/anacri.h"
#include "asterfort/avgrno.h"
#include "asterfort/dtauno.h"
#include "asterfort/fgdoba.h"
#include "asterfort/fgdohs.h"
#include "asterfort/fgdowh.h"
#include "asterfort/fmcros.h"
#include "asterfort/fmpapa.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/rcvale.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/tbnuli.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!     COMMANDE POST_FATIGUE
!              CHARGEMENT PERIODIQUE
!     -----------------------------------------------------------------
!
    integer(kind=8) :: n1, n2, n3, n4, n5, n6, nbf, nbptot, nbpts, i, nbc, ibid, iordo
    integer(kind=8) :: ifonc1, ifonc, ilign, nbpar, nbpapf, j, nval, paract(35), nbeps
    integer(kind=8) :: ifonc2, ifonce, iordoe, ifonc3, ifoncp, iordop, nbepsp
    integer(kind=8) :: tdisp, nbnop, lisnoe(1), nbnot, nbordr, nnoini
    integer(kind=8) :: tspaq, k, jrwork, nbcmp, ordini
!
    real(kind=8) :: rbid, phmax, cissio, sphere, pcorr, val(2), vmax, vmin
    real(kind=8) :: domage(1), rcrit, vresu(24), resu(7), valpar(35)
    complex(kind=8) :: cbid
    aster_logical :: lhaigh, lke, post, fordef, plcicr
    aster_logical :: crsigm, crepst, crepse, crepsp, plcr2
    integer(kind=8) :: icodre(2), icodwo, icodba, icodhs
    character(len=8) :: k8b, nomten(6), nompar, nommat
    character(len=8) :: result, nomeps(6), nomepp(6)
    character(len=16) :: kdomm, nomcmd, pheno, criter, nomfor, typcha, forvie
    character(len=16) :: nomres(2), proaxe, nommet, forcri, cara
    character(len=19) :: k19b
    character(len=24) :: fvale(6), etvale(6), ptvale(6)
!     --- POST_FATI_MULT -----------------------------------------------
    parameter(nbpapf=50)
    character(len=3) :: typppf(nbpapf)
    character(len=16) :: nomppf(nbpapf)
    data nomppf/'CRITERE', 'VALE_CRITERE', 'PRES_HYDRO_MAX',&
     &               'AMPLI_CISSION', 'RAYON_SPHERE',&
     &               'VALE_MIN', 'VALE_MAX', 'DOMMAGE', 'NBRUP',&
     &               'DTAUMA', 'PHYDRM', 'NORMAX', 'NORMOY',&
     &               'EPNMAX', 'EPNMOY', 'DEPSPE', 'EPSPR1',&
     &               'SIGNM1', 'DENDIS', 'DENDIE', 'APHYDR',&
     &               'MPHYDR', 'DSIGEQ', 'SIGPR1', 'EPSNM1',&
     &               'INVA2S', 'DSITRE', 'DEPTRE', 'EPSPAC',&
     &               'RAYSPH', 'AMPCIS', 'DEPSEE', &
     &               'DTAUCR', 'DGAMCR', 'DSINCR', 'DEPNCR', &
     &               'MTAUCR', 'MGAMCR', 'MSINCR', 'MEPNCR', &
     &               'DGAMPC', 'DEPNPC', 'MGAMPC', 'MEPNPC', &
     &               'VNM1X', 'VNM1Y', 'VNM1Z',&
     &               'VNM2X', 'VNM2Y', 'VNM2Z'/
!
!
    data typppf/'K16', 'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R',&
     &                  'R', 'R', 'R', 'R', 'R', 'R', 'R'/
!
!     ---------------------------------------------------------------
!     ----------------------------------------------------------------
!
    call jemarq()
!
    cbid = (0.d0, 0.d0)
    rbid = 0.d0
    ibid = 0
!
    lke = .false.
    lhaigh = .false.
    nbc = 1
!
!
    call getres(result, k8b, nomcmd)
!
!     --- DETERMINATION DES CRITERES---
!
    criter = ' '
    call getvtx(' ', 'CRITERE', scal=criter, nbret=n1)
!
    typcha = ' '
    call getvtx(' ', 'TYPE_CHARGE', scal=typcha, nbret=n1)
!
    call getvid(' ', 'FORMULE_GRDEQ', scal=nomfor, nbret=nval)
    if (nval .eq. 0) then
        nomfor = '        '
    end if
!
    call getvid(' ', 'FORMULE_VIE', scal=forvie, nbret=nval)
    if (nval .eq. 0) then
        forvie = '        '
    end if
!
    call getvid(' ', 'FORMULE_CRITIQUE', scal=forcri, nbret=nval)
    if (nval .eq. 0) then
        forcri = '        '
    end if
!
!
    kdomm = ' '
    call getvtx(' ', 'DOMMAGE', scal=kdomm, nbret=n1)
!
! ---   NOM DE LA METHODE PERMETTANT DE DETERMINER LE CERCLE CIRCONSCRIT
    call getvtx(' ', 'METHODE', scal=nommet, nbret=nval)
    if (nval .eq. 0) then
        nommet = '        '
    end if
!
! ---   PROJECTION SUR UN AXE OU SUR DEUX AXES
!     (CHARGEMENT NON_PERIODIQUE UNIQUEMENT)
    call getvtx(' ', 'PROJECTION', scal=proaxe, nbret=nval)
    if (nval .eq. 0) then
        proaxe = '        '
    end if
!
!---    ANALYSER LE CRITERE
!  INITIALISER
    crsigm = .false.
    crepst = .false.
    crepse = .false.
    crepsp = .false.
    fordef = .false.
!
    call anacri(criter, nomfor, typcha, 'OUI', paract, &
                fordef, crsigm, crepst, crepse, crepsp)
!     --- RECUPERATION DE LA FONCTION CHARGEMENT ---
!
!CCCCCCCCC RECUPERER LA CONTRAINTE
    call getvid('HISTOIRE', 'SIGM_XX', iocc=1, scal=nomten(1), nbret=n1)
    call getvid('HISTOIRE', 'SIGM_YY', iocc=1, scal=nomten(2), nbret=n2)
    call getvid('HISTOIRE', 'SIGM_ZZ', iocc=1, scal=nomten(3), nbret=n3)
    call getvid('HISTOIRE', 'SIGM_XY', iocc=1, scal=nomten(4), nbret=n4)
    call getvid('HISTOIRE', 'SIGM_XZ', iocc=1, scal=nomten(5), nbret=n5)
    call getvid('HISTOIRE', 'SIGM_YZ', iocc=1, scal=nomten(6), nbret=n6)
    nbf = n1+n2+n3+n4+n5+n6
!
    if (nbf .ne. 0) then
        fvale(1) = nomten(1)//'           .VALE'
        call jelira(fvale(1), 'LONMAX', nbpts)
    end if
!
!CCCCCCCCC RECUPERER LA DEFORMATION TOTALE
    call getvid('HISTOIRE', 'EPS_XX', iocc=1, scal=nomeps(1), nbret=n1)
    call getvid('HISTOIRE', 'EPS_YY', iocc=1, scal=nomeps(2), nbret=n2)
    call getvid('HISTOIRE', 'EPS_ZZ', iocc=1, scal=nomeps(3), nbret=n3)
    call getvid('HISTOIRE', 'EPS_XY', iocc=1, scal=nomeps(4), nbret=n4)
    call getvid('HISTOIRE', 'EPS_XZ', iocc=1, scal=nomeps(5), nbret=n5)
    call getvid('HISTOIRE', 'EPS_YZ', iocc=1, scal=nomeps(6), nbret=n6)
    nbeps = n1+n2+n3+n4+n5+n6
!
    if (nbeps .ne. 0) then
        etvale(1) = nomeps(1)//'           .VALE'
        call jelira(etvale(1), 'LONMAX', nbpts)
    end if
!
!CCCCCCCCC RECUPERER LA DEFORMATION PLASTIQUE
    call getvid('HISTOIRE', 'EPSP_XX', iocc=1, scal=nomepp(1), nbret=n1)
    call getvid('HISTOIRE', 'EPSP_YY', iocc=1, scal=nomepp(2), nbret=n2)
    call getvid('HISTOIRE', 'EPSP_ZZ', iocc=1, scal=nomepp(3), nbret=n3)
    call getvid('HISTOIRE', 'EPSP_XY', iocc=1, scal=nomepp(4), nbret=n4)
    call getvid('HISTOIRE', 'EPSP_XZ', iocc=1, scal=nomepp(5), nbret=n5)
    call getvid('HISTOIRE', 'EPSP_YZ', iocc=1, scal=nomepp(6), nbret=n6)
    nbepsp = n1+n2+n3+n4+n5+n6
!
    if (nbepsp .ne. 0) then
        ptvale(1) = nomepp(1)//'           .VALE'
        call jelira(ptvale(1), 'LONMAX', nbpts)
    end if
!
!C  CONTRUIRE TABLEAU CONTRAINTE
    if (nbf .eq. 0) then
        if (crsigm) then
            call utmess('F', 'FATIGUE1_97')
        end if
        call wkvect('&&POFAPE.ORDO', 'V V R', nbpts/2*6, iordo)
    else
!
        nbptot = nbpts
        do i = 2, nbf
            fvale(i) = nomten(i)//'           .VALE'
            call jelira(fvale(i), 'LONMAX', nbpts)
            if (nbpts .ne. nbptot) then
                call utmess('F', 'FATIGUE1_21')
            end if
        end do
        call wkvect('&&POFAPE.ORDO', 'V V R', nbptot/2*nbf, iordo)
        call jeveuo(fvale(1), 'L', ifonc1)
        do i = 2, nbf
            call jeveuo(fvale(i), 'L', ifonc)
            do j = 1, nbptot/2
                if (zr(ifonc+j-1) .ne. zr(ifonc1+j-1)) then
                    call utmess('F', 'FATIGUE1_21')
                end if
                zr(iordo+(j-1)*nbf+i-1) = zr(ifonc+nbptot/2+j-1)
            end do
        end do
        nbptot = nbptot/2
        do j = 1, nbptot
            zr(iordo+(j-1)*nbf) = zr(ifonc1+nbptot+j-1)
!
        end do
!
    end if
!
!C  CONTRUIRE TABLEAU DEFORMATION TOTALE
    if (nbeps .eq. 0) then
        if (crepst) then
            call utmess('F', 'FATIGUE1_98')
        end if
        call wkvect('&&POFAPE.ORDOE', 'V V R', nbpts/2*6, iordoe)
    else
!
        nbptot = nbpts
        do i = 2, nbeps
            etvale(i) = nomeps(i)//'           .VALE'
            call jelira(etvale(i), 'LONMAX', nbpts)
            if (nbpts .ne. nbptot) then
                call utmess('F', 'FATIGUE1_21')
            end if
        end do
        call wkvect('&&POFAPE.ORDOE', 'V V R', nbptot*nbeps/2, iordoe)
        call jeveuo(etvale(1), 'L', ifonc2)
        do i = 2, nbeps
            call jeveuo(etvale(i), 'L', ifonce)
            do j = 1, nbptot/2
                if (zr(ifonce+j-1) .ne. zr(ifonc2+j-1)) then
                    call utmess('F', 'FATIGUE1_21')
                end if
                zr(iordoe+(j-1)*nbeps+i-1) = zr(ifonce+nbptot/2+j-1)
            end do
        end do
        nbptot = nbptot/2
        do j = 1, nbptot
            zr(iordoe+(j-1)*nbeps) = zr(ifonc2+nbptot+j-1)
        end do
    end if
!
!
!C  CONTRUIRE TABLEAU DEFORMATION PLASTIQUE
!
    if (nbepsp .eq. 0) then
        if (crepsp) then
            call utmess('F', 'FATIGUE1_99')
        end if
        call wkvect('&&POFAPE.ORDOP', 'V V R', nbpts/2*6, iordop)
    else
!
        nbptot = nbpts
        do i = 2, nbepsp
            ptvale(i) = nomepp(i)//'           .VALE'
            call jelira(ptvale(i), 'LONMAX', nbpts)
            if (nbpts .ne. nbptot) then
                call utmess('F', 'FATIGUE1_21')
            end if
        end do
        call wkvect('&&POFAPE.ORDOP', 'V V R', nbptot*nbepsp/2, iordop)
        call jeveuo(ptvale(1), 'L', ifonc3)
        do i = 2, nbepsp
            call jeveuo(ptvale(i), 'L', ifoncp)
            do j = 1, nbptot/2
                if (zr(ifoncp+j-1) .ne. zr(ifonc3+j-1)) then
                    call utmess('F', 'FATIGUE1_21')
                end if
                zr(iordop+(j-1)*nbeps+i-1) = zr(ifoncp+nbptot/2+j-1)
            end do
        end do
        nbptot = nbptot/2
        do j = 1, nbptot
            zr(iordop+(j-1)*nbeps) = zr(ifonc3+nbptot+j-1)
        end do
    end if
!
!CC  RECUPERER LE MATERIAU
!
    nommat = ' '
    call getvid(' ', 'MATER', scal=nommat, nbret=n1)
!
    if (crepse) then
        if ((nbeps+nbepsp) .eq. 0) then
            call utmess('F', 'FATIGUE1_95')
        end if
        if ((nbeps+nbepsp) .gt. 0) then
            call utmess('A', 'FATIGUE1_96')
        end if
    end if
!
!
!     --- CREATION DE LA TABLE ---
!
    call tbcrsd(result, 'G')
    call tbajpa(result, nbpapf, nomppf, typppf)
!
!
!
!CCCCCCCCCCCCCCCCCC
!
    call tbajli(result, 1, nomppf(1), [ibid], [rbid], &
                [cbid], criter, 0)
    call tbnuli(result, 1, nomppf(1), [ibid], [rbid], &
                [cbid], criter, [rbid], k8b, ilign)
    if (ilign .le. 0) ilign = 0
!
!
!
    do j = 1, 7
        resu(j) = 0.0d0
    end do
!
    if ((criter .eq. 'FORMULE_CRITERE') .or. (criter .eq. 'MATAKE_MODI_AV') .or. &
        (criter .eq. 'DANG_VAN_MODI_AV') .or. (criter .eq. 'FATESOCI_MODI_AV') .or. &
        (criter .eq. 'MATAKE_MODI_AC') .or. (criter .eq. 'DANG_VAN_MODI_AC')) then
!
! ! ANALYSER LE CRITERE
!         call anacri(criter, nomfor, typcha, 'OUI', paract,&
!                     fordef, lbid, lbid, lbid, lbid)
        post = .true.
! CONS TRUIRE UN VECTEUR WORK QUI CONTIENT CONTRAINE ET DEFORMATION
        nbcmp = 6
!
        call wkvect('&&POFAPE.ORDOCD', 'V V R', nbptot*nbcmp*3, jrwork)
!
        do j = 1, nbptot
            do k = 1, 6
                zr(jrwork+(j-1)*nbcmp*3+k-1) = zr(iordo+(j-1)*nbcmp+k-1)
                zr(jrwork+(j-1)*nbcmp*3+nbcmp+k-1) = zr(iordoe+(j-1)*nbcmp+k-1)
                zr(jrwork+(j-1)*nbcmp*3+nbcmp*2+k-1) = zr(iordop+(j-1)*nbcmp+k-1)
            end do
        end do
!
        tdisp = nbptot*nbcmp*3
        nbnot = 1
        lisnoe(1) = 1
        nbordr = nbpts/2
        nnoini = 1
        ordini = 1
        nbnop = 1
        tspaq = 18
        plcicr = .false.
        plcr2 = .false.
!
! POUR CHARGEMENT PERIODIQUE
        if (typcha .eq. 'PERIODIQUE') then
!
            call dtauno(jrwork, lisnoe, nbnot, nbordr, ordini, &
                        nnoini, nbnop, tspaq, nommet, criter, &
                        nomfor, kdomm, forvie, forcri, k8b, &
                        k19b, nommat, post, valpar, vresu)
!
            if ((paract(1) .eq. 1) .or. (paract(3) .eq. 1) .or. (paract(4) .eq. 1) .or. &
                (paract(5) .eq. 1) .or. (paract(6) .eq. 1)) then
!
                plcicr = .true.
            end if
!
! PLAN CRITIQUE DE TYPE CISSAILLEMENT DE DANG VAN-MATAKE
            if (plcicr) then
!
                call tbajli(result, 1, nomppf(10), [ibid], vresu(1), &
                            [cbid], k8b, ilign)
!
                do i = 1, 3
                    call tbajli(result, 1, nomppf(i+44), [ibid], vresu(i+1), &
                                [cbid], k8b, ilign)
!
                end do
!
                do i = 1, 4
                    call tbajli(result, 1, nomppf(i+11), [ibid], vresu(i+4), &
                                [cbid], k8b, ilign)
!
                end do
            end if
!
!POUR LES NOUVEUAX CRITERS DE PLAN CRITIQUE
!
            do i = 24, 35
                if (paract(i) .eq. 1) plcr2 = .true.
            end do
!
            if (plcr2) then
                do i = 1, 3
                    call tbajli(result, 1, nomppf(i+44), [ibid], vresu(i+1), &
                                [cbid], k8b, ilign)
!
                end do
!
                do i = 24, 35
                    if (paract(i) .eq. 1) then
                        call tbajli(result, 1, nomppf(i+9), [ibid], valpar(i), &
                                    [cbid], k8b, ilign)
                    end if
                end do
!
            end if
!
!
! POUR LES GRANDEURS HORS DES CRITERES A PLAN CRITIQUE
            if (paract(2) .eq. 1) then
                call tbajli(result, 1, nomppf(11), [ibid], valpar(2), &
                            [cbid], k8b, ilign)
            end if
            do i = 7, 23
                if (paract(i) .eq. 1) then
                    call tbajli(result, 1, nomppf(i+9), [ibid], valpar(i), &
                                [cbid], k8b, ilign)
                end if
            end do
!
!
            call tbajli(result, 1, nomppf(2), [ibid], vresu(9), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(9), [ibid], vresu(10), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(8), [ibid], vresu(11), &
                        [cbid], k8b, ilign)
!
!
! POUR CHARGEMENT NON-PERIODIQUE
        else if (typcha .eq. 'NON_PERIODIQUE') then
!
            call avgrno(zr(jrwork), tdisp, lisnoe, nbnot, nbordr, &
                        nnoini, nbnop, tspaq, criter, nomfor, &
                        kdomm, forvie, fordef, k8b, proaxe, &
                        nommat, k19b, post, resu)
!
            call tbajli(result, 1, nomppf(45), [ibid], resu(1), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(46), [ibid], resu(2), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(47), [ibid], resu(3), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(8), [ibid], resu(4), &
                        [cbid], k8b, ilign)
!
            call tbajli(result, 1, nomppf(48), [ibid], resu(5), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(49), [ibid], resu(6), &
                        [cbid], k8b, ilign)
            call tbajli(result, 1, nomppf(50), [ibid], resu(7), &
                        [cbid], k8b, ilign)
!
        end if
!
        goto 50
    end if
!
!
    nomres(1) = 'D0'
    nomres(2) = 'TAU0'
    nbpar = 0
    nompar = ' '
    call rcvale(nommat, 'FATIGUE ', nbpar, nompar, [rbid], &
                2, nomres, val, icodre, 2)
!
    if (criter .eq. 'CROSSLAND') then
!          -----------------------
        call fmcros(nbf, nbptot, zr(iordo), val(1), val(2), &
                    rcrit, phmax, cissio)
!
        call tbajli(result, 1, nomppf(2), [ibid], [rcrit], &
                    [cbid], k8b, ilign)
        call tbajli(result, 1, nomppf(3), [ibid], [phmax], &
                    [cbid], k8b, ilign)
        call tbajli(result, 1, nomppf(4), [ibid], [cissio], &
                    [cbid], k8b, ilign)
!
    else if (criter .eq. 'PAPADOPOULOS') then
!              --------------------------
        call fmpapa(nbf, nbptot, zr(iordo), val(1), val(2), &
                    rcrit, phmax, sphere)
!
        call tbajli(result, 1, nomppf(2), [ibid], [rcrit], &
                    [cbid], k8b, ilign)
        call tbajli(result, 1, nomppf(3), [ibid], [phmax], &
                    [cbid], k8b, ilign)
        call tbajli(result, 1, nomppf(5), [ibid], [sphere], &
                    [cbid], k8b, ilign)
!
    end if
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
    if (criter .ne. 'FORMULE_CRITERE') then
!     --- CORRECTION POUR CALCUL DU DOMMAGE ----
!
        call getvr8(' ', 'COEF_CORR', scal=pcorr, nbret=n1)
        if (n1 .ne. 0) then
            vmax = 2.d0*(rcrit+val(2))*pcorr
            vmin = 0.d0
        else
            vmax = 2.d0*(rcrit+val(2))*(val(1)/val(2))
            vmin = 0.d0
        end if
        call tbajli(result, 1, nomppf(6), [ibid], [vmin], &
                    [cbid], k8b, ilign)
        call tbajli(result, 1, nomppf(7), [ibid], [vmax], &
                    [cbid], k8b, ilign)
!
!         --- CALCUL DU DOMMAGE ELEMENTAIRE ---
!
!
!         --- CALCUL DU DOMMAGE ELEMENTAIRE DE WOHLER ---
!             ---------------------------------------
        if (kdomm .eq. 'WOHLER') then
            pheno = 'FATIGUE'
            call rccome(nommat, pheno, icodre(1))
            if (icodre(1) .eq. 1) then
                call utmess('F', 'FATIGUE1_24')
            end if
            cara = 'WOHLER'
            call rcpare(nommat, pheno, cara, icodwo)
            cara = 'A_BASQUIN'
            call rcpare(nommat, pheno, cara, icodba)
            cara = 'A0'
            call rcpare(nommat, pheno, cara, icodhs)
            if (icodwo .eq. 0) then
                call fgdowh(nommat, nbc, [vmin], [vmax], lke, &
                            [0.d0], lhaigh, [0.d0], domage)
            else if (icodba .eq. 0) then
                call fgdoba(nommat, nbc, [vmin], [vmax], lke, &
                            [0.d0], lhaigh, [0.d0], domage)
            else if (icodhs .eq. 0) then
                call fgdohs(nommat, nbc, [vmin], [vmax], lke, &
                            [0.d0], lhaigh, [0.d0], domage)
            end if
!
            call tbajli(result, 1, nomppf(8), [ibid], [domage], &
                        [cbid], k8b, ilign)
!
        else if (kdomm .eq. ' ') then
        else
            call utmess('F', 'FATIGUE1_20')
        end if
!
    end if
!
50  continue
!
    call jedetr('&&POFAPE.ORDO')
    call jedema()
!
end subroutine
