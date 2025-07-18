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
subroutine specep(casint, nomu, spectr, base, vite, &
                  nuor, imodi, imodf, nbm, nbpf)
    implicit none
!     PROJECTION D'UN SPECTRE D'EXCITATION TURBULENTE LOCALISEE (FORCES
!     ET MOMENTS PONCTUELS) SUR UNE BASE MODALE PERTURBEE PAR PRISE EN
!     COMPTE DU COUPLAGE FLUIDE STRUCTURE
!     APPELANT : OP0146 , OPERATEUR PROJ_SPEC_BASE
!-----------------------------------------------------------------------
! IN  : CASINT  : BOOLEEN, DONNE L'OPTION DE CALCUL
!       CASINT  = .TRUE.  => CALCUL DE TOUS LES INTERSPECTRES
!       CASINT  = .FALSE. => CALCUL DES AUTOSPECTRES UNIQUEMENT
! IN  : NOMU    : NOM UTILISATEUR
! IN  : SPECTR  : NOM DU CONCEPT SPECTRE
! IN  : BASE    : NOM DU CONCEPT MELASFLU
! IN  : VITE    : VITESSE ETUDIEE
! IN  : NUOR    : NUMEROS D'ORDRE DES MODES DU CONCEPT MELASFLU
! IN  : IMODI   : INDICE DU PREMIER MODE PRIS EN COMPTE
! IN  : IMODF   : INDICE DU DERNIER MODE PRIS EN COMPTE
! IN  : NBM     : NOMBRE DE MODES DU CONCEPT MELASFLU
! IN  : NBPF    : NOMBRE DE POINTS DE LA DISCRETISATION FREQUENTIELLE
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/axdipo.h"
#include "asterfort/deelpo.h"
#include "asterfort/dismoi.h"
#include "asterfort/exmano.h"
#include "asterfort/fointr.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/scalep.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    aster_logical :: casint
    integer(kind=8) :: imodi, imodf, nbm, nuor(nbm), nbpf, ij, nbval
    character(len=8) :: nomu
    character(len=19) :: spectr, base
    real(kind=8) :: vite
!
    integer(kind=8) :: ibid, dim, ival(2)
    real(kind=8) :: module
    real(kind=8) :: coefac(8), coefae(8), coefdc(6), coefde(6)
    aster_logical :: ltable, exiind
    character(len=8) :: caelem, modele, table, noma, nomno0
    character(len=16) :: config, nopart(2)
    character(len=19) :: typflu
    character(len=24) :: spvain, spvate, spvare, spnnoe
    character(len=24) :: chvale, chnumj, chtab
    character(len=24) :: remf, fsic, chrefe, mlgnma, chnumi
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iaxe, ideb, idec, iex, iex1
    integer(kind=8) :: iex2, ifsic, iinte, il, im1, im1b, im2
    integer(kind=8) :: im2b, imail, inat, iremf, iscal, ispin
    integer(kind=8) :: ispno, ispre, ispte, itypfl, iv, ivale, lwr
    integer(kind=8) :: nbexcp, nbma, nbmano, nbmr, numno0, lnumi, lnumj, i1, ind
    integer(kind=8) :: iprol, ire, iim, itab, nbfreq, isre, isim, ier2
    real(kind=8) :: beta, coedim, coef1, coef2, coefd, difphi, fr
    real(kind=8) :: frc, frref, phi1, phi2, phie
    real(kind=8) :: rhof, s0, scal11, scal12, scal21, scal22
    real(kind=8) :: sref, tolr, uabs
    real(kind=8), pointer :: freq(:) => null()
!-----------------------------------------------------------------------
    data nopart/'NUME_ORDRE_I', 'NUME_ORDRE_J'/
!
    data coefac/1.d-4, 1.9d-1, 7.d-2, 1.6d0,&
     &              2.7d-5, 1.9d-1, 7.d-2, 2.1d0/
!
    data coefae/1.d-5, 5.d-1, 2.d-2, 2.9d0,&
     &              3.3d-4, 1.d-1, 2.d-2, 2.8d0/
!
    data coefdc/4.d-5, 1.9d-1, 1.6d0,&
     &              2.7d-5, 1.9d-1, 2.1d0/
!
    data coefde/1.7d-5, 2.2d-1, 2.9d0,&
     &              1.7d-5, 1.9d-1, 2.8d0/
!
!-----------------------------------------------------------------------
    call jemarq()
!
!
! --- 1.TEST DE COMPATIBILITE TYPE DE SPECTRE/CONFIGURATION ETUDIEE ---
!
    remf = base//'.REMF'
    call jeveuo(remf, 'L', iremf)
    typflu = zk8(iremf)
    fsic = typflu//'.FSIC'
    call jeveuo(fsic, 'L', ifsic)
    itypfl = zi(ifsic)
    if (itypfl .ne. 2) then
        call utmess('F', 'MODELISA7_4')
    end if
!
!
! --- 2.RECUPERATION DU NOM DU CONCEPT MAILLAGE ---
!
    iv = 1
    write (chrefe, '(A8,A5,2I3.3,A5)') base(1:8), '.C01.', nuor(1), iv,&
     &                                 '.REFE'
    call dismoi("NOM_MAILLA", chrefe, 'CHAM_NO', repk=noma)
!
!
! --- 3.RECUPERATION DES INFORMATIONS CARACTERISTIQUES DU SPECTRE ---
!
    spvain = spectr//'.VAIN'
    spvate = spectr//'.VATE'
    spvare = spectr//'.VARE'
    spnnoe = spectr//'.NNOE'
!
    call jeveuo(spvain, 'L', ispin)
    ltable = .false.
    if (zi(ispin+1) .eq. 0) ltable = .true.
!
    call jeveuo(spvate, 'L', ispte)
    caelem = zk16(ispte+1) (1:8)
    modele = zk16(ispte+2) (1:8)
!
    if (ltable) then
!
        table = zk16(ispte+3) (1:8)
        chnumi = table//'.NUMI'
        chnumj = table//'.NUMJ'
        call jeveuo(chnumi, 'L', lnumi)
        call jeveuo(chnumj, 'L', lnumj)
        call jelira(chnumi, 'LONMAX', nbexcp)
!
    else
!
        nbexcp = 2
        config = zk16(ispte+4)
        call jeveuo(spvare, 'L', ispre)
        rhof = zr(ispre)
        call jeveuo(spnnoe, 'L', ispno)
        nomno0 = zk8(ispno)
!
!-------RECUPERATION DU DIAMETRE EXTERIEUR DE LA POUTRE, NECESSAIRE AU
!       DIMENSIONNEMENT DE L'EXCITATION GRAPPE2
!
        numno0 = char8_to_int(nomno0)
        mlgnma = noma//'.TYPMAIL'
        call jelira(mlgnma, 'LONMAX', nbma)
        call wkvect('&&SPECEP.TEMP.MAIL', 'V V I', nbma, imail)
        call exmano(noma, numno0, zi(imail), nbmano)
        if (nbmano .ne. 2) then
            call utmess('F', 'ALGELINE_70')
        end if
        call deelpo(caelem, noma, zi(imail), phi1)
        call deelpo(caelem, noma, zi(imail+1), phi2)
        tolr = r8prem()
        difphi = dble(abs(phi1-phi2))
        if (difphi .gt. phi1*tolr) then
            call utmess('F', 'ALGELINE_71')
        else
            phie = phi1
        end if
!
!-------CALCUL DE COEFFICIENTS DE DIMENSIONNEMENT
!
        coef1 = 282.d0/890.d0
        coef2 = 0.77d0*0.77d0/4.d0
        coefd = phie/8.9d-2
!
    end if
!
!
! --- 4.DETERMINATION DE L'AXE DIRECTEUR DE LA POUTRE ---
!
    call axdipo(noma, caelem, modele, iaxe)
!
!
! --- 5.CALCUL DES PRODUITS SCALAIRES PHII(XK).NK ET PHII'(XM).NM ---
! ---   XK POINTS D'APPLICATION DES FORCES PONCTUELLES            ---
! ---   XM POINTS D'APPLICATION DES MOMENTS PONCTUELS             ---
! ---   NK ET NM DIRECTIONS D'APPLICATION DES EXCITATIONS         ---
!
    nbmr = imodf-imodi+1
    call wkvect('&&SPECEP.TEMP.SCAL', 'V V R', nbexcp*nbmr, iscal)
    call scalep(spectr, noma, base, nuor, nbm, &
                imodi, nbmr, nbexcp, ltable, iaxe, &
                zr(iscal))
!
!
! --- 6.CALCUL DES INTERSPECTRES D'EXCITATIONS MODALES ---
! ---   BOUCLE SUR LE NOMBRE DE VITESSES               ---
!
    call wkvect('&&SPECEP.TEMP.SWR ', 'V V R', nbpf, lwr)
    dim = nbexcp*(nbexcp+1)/2
    dim = 2*nbpf*dim
    call wkvect('&&SPECEP.TEMP.INTE', 'V V R', dim, iinte)
!
! --- 6.1.RECUPERATION DE LA DISCRETISATION FREQUENTIELLE
    call jeveuo(nomu//'.DISC', 'L', lwr)
!
! --- 6.2.INTERPOLATION DES INTERSPECTRES A PROJETER
!
    if (ltable) then
!
        call wkvect('&&SPECEP.PROL', 'V V K24', 5, iprol)
        zk24(iprol) = 'FONCTION'
        zk24(iprol+1) = 'LIN LIN '
        zk24(iprol+2) = ' '
        zk24(iprol+3) = 'TOUTRESU'
        zk24(iprol+4) = 'CC      '
        call wkvect('&&SPECEP.IRE', 'V V R', nbpf, ire)
        call wkvect('&&SPECEP.IIM', 'V V R', nbpf, iim)
        chtab = table//'.VALE'
        call jeveuo(table//'.DISC', 'L', vr=freq)
!
        do iex2 = 1, nbexcp
            ival(2) = iex2
            do iex1 = 1, iex2
                iex = iex2*(iex2-1)/2+iex1
                ival(1) = iex1
                exiind = .false.
                do i1 = 1, nbexcp
                    if ((zi(lnumi-1+i1) .eq. ival(1)) .and. (zi(lnumj-1+i1) .eq. ival(2))) then
                        exiind = .true.
                        ind = i1
                    end if
                end do
                if (.not. exiind) then
                    call utmess('F', 'MODELISA2_89')
                end if
                call jeveuo(jexnum(chtab, ind), 'L', itab)
                call jelira(jexnum(chtab, ind), 'LONMAX', nbval)
                if (iex2 .eq. iex1) then
                    nbfreq = nbval
                    call jeexin('&&SPECEP.SRE', ibid)
                    if (ibid .eq. 0) then
                        call wkvect('&&SPECEP.SRE', 'V V R', nbfreq, isre)
                        call wkvect('&&SPECEP.SIM', 'V V R', nbfreq, isim)
                    end if
                    call fointr(' ', zk24(iprol), nbfreq, freq, zr(itab), &
                                nbpf, zr(lwr), zr(isre), ier2)
                    do i1 = 1, nbfreq
                        zr(isim-1+i1) = 0.d0
                    end do
                else
                    nbfreq = nbval/2
                    call jeexin('&&SPECEP.SRE', ibid)
                    if (ibid .eq. 0) then
                        call wkvect('&&SPECEP.SRE', 'V V R', nbfreq, isre)
                        call wkvect('&&SPECEP.SIM', 'V V R', nbfreq, isim)
                    end if
                    do i1 = 1, nbfreq
                        zr(ire-1+i1) = zr(itab+2*(i1-1))
                        zr(iim-1+i1) = zr(itab+2*(i1-1)+1)
                    end do
                    call fointr(' ', zk24(iprol), nbfreq, freq, zr(ire), &
                                nbpf, zr(lwr), zr(isre), ier2)
                    call fointr(' ', zk24(iprol), nbfreq, freq, zr(iim), &
                                nbpf, zr(lwr), zr(isim), ier2)
                end if
                do il = 1, nbpf
                    idec = 2*nbpf*(iex-1)+2*(il-1)
                    zr(iinte+idec) = zr(isre-1+il)
                    zr(iinte+idec+1) = zr(isim-1+il)
                end do
            end do
        end do
!
    else if (config(1:7) .eq. 'ASC_CEN') then
!
        uabs = dble(abs(vite))
!
        do iex2 = 1, nbexcp
            sref = coefac(4*(iex2-1)+1)
            frref = coefac(4*(iex2-1)+2)
            frc = coefac(4*(iex2-1)+3)
            beta = coefac(4*(iex2-1)+4)
            s0 = sref*(1.d0+(frref/frc)**(beta))
            inat = iex2-int(iex2/2)*2
            coedim = coef1*coefd*coefd
            if (inat .eq. 0) coedim = coedim*coefd*coefd*coef2
            iex = iex2*(iex2+1)/2
            do il = 1, nbpf
                idec = 2*nbpf*(iex-1)+2*(il-1)
                fr = zr(lwr+il-1)*phie/uabs
                module = 1.d0+(fr/frc)**(beta)
                module = s0/module
                zr(iinte+idec) = coedim*module
            end do
        end do
!
    else if (config(1:7) .eq. 'ASC_EXC') then
!
        uabs = dble(abs(vite))
!
        do iex2 = 1, nbexcp
            sref = coefae(4*(iex2-1)+1)
            frref = coefae(4*(iex2-1)+2)
            frc = coefae(4*(iex2-1)+3)
            beta = coefae(4*(iex2-1)+4)
            s0 = sref*(1.d0+(frref/frc)**(beta))
            inat = iex2-int(iex2/2)*2
            coedim = coef1*coefd*coefd
            if (inat .eq. 0) coedim = coedim*coefd*coefd*coef2
            iex = iex2*(iex2+1)/2
            do il = 1, nbpf
                idec = 2*nbpf*(iex-1)+2*(il-1)
                fr = zr(lwr+il-1)*phie/uabs
                module = 1.d0+(fr/frc)**(beta)
                module = s0/module
                zr(iinte+idec) = coedim*module
            end do
        end do
!
    else if (config(1:7) .eq. 'DES_CEN') then
!
        uabs = dble(abs(vite))
!
        do iex2 = 1, nbexcp
            s0 = coefdc(3*(iex2-1)+1)
            frc = coefdc(3*(iex2-1)+2)
            beta = coefdc(3*(iex2-1)+3)
            inat = iex2-int(iex2/2)*2
            coedim = coef1*coefd*coefd
            if (inat .eq. 0) coedim = coedim*coefd*coefd*coef2
            iex = iex2*(iex2+1)/2
            do il = 1, nbpf
                idec = 2*nbpf*(iex-1)+2*(il-1)
                fr = zr(lwr+il-1)*phie/uabs
                module = 1.d0+(fr/frc)**(beta)
                module = s0/module
                zr(iinte+idec) = coedim*module
            end do
        end do
!
    else if (config(1:7) .eq. 'DES_EXC') then
!
        uabs = dble(abs(vite))
!
        do iex2 = 1, nbexcp
            s0 = coefde(3*(iex2-1)+1)
            frc = coefde(3*(iex2-1)+2)
            beta = coefde(3*(iex2-1)+3)
            inat = iex2-int(iex2/2)*2
            coedim = coef1*coefd*coefd
            if (inat .eq. 0) coedim = coedim*coefd*coefd*coef2
            iex = iex2*(iex2+1)/2
            do il = 1, nbpf
                idec = 2*nbpf*(iex-1)+2*(il-1)
                fr = zr(lwr+il-1)*phie/uabs
                module = 1.d0+(fr/frc)**(beta)
                module = s0/module
                zr(iinte+idec) = coedim*module
            end do
        end do
!
    end if
!
! --- 6.3.PROJECTION DES INTERSPECTRES
!
    ij = 0
    chvale = nomu//'.VALE'
    do im2 = imodi, imodf
        ideb = im2
        if (casint) ideb = imodi
        do im1 = ideb, im2
            ij = ij+1
            call jeveuo(jexnum(chvale, ij), 'E', ivale)
            call jelira(jexnum(chvale, ij), 'LONMAX', nbval)
!
            im2b = im2-imodi+1
            im1b = im1-imodi+1
!
            if (ltable) then
!
                do il = 1, nbpf
!
                    do iex2 = 1, nbexcp
                        scal12 = zr(iscal+nbexcp*(im1b-1)+iex2-1)
                        scal22 = zr(iscal+nbexcp*(im2b-1)+iex2-1)
                        iex = iex2*(iex2+1)/2
                        idec = 2*nbpf*(iex-1)+2*(il-1)
                        if (nbval .eq. nbpf) then
                            zr(ivale+il-1) = zr(ivale+il-1)+scal12*scal22*zr(iinte+idec)
                        else
                            zr(ivale+2*(il-1)) = zr( &
                                                 ivale+2*(il-1))+scal12*scal22*zr(iinte+idec)
                        end if
                    end do
!
                    if (nbexcp .gt. 1) then
                        do iex2 = 2, nbexcp
                            scal12 = zr(iscal+nbexcp*(im1b-1)+iex2-1)
                            scal22 = zr(iscal+nbexcp*(im2b-1)+iex2-1)
                            do iex1 = 1, iex2-1
                                scal11 = zr(iscal+nbexcp*(im1b-1)+iex1-1)
                                scal21 = zr(iscal+nbexcp*(im2b-1)+iex1-1)
                                iex = iex2*(iex2-1)/2+iex1
                                idec = 2*nbpf*(iex-1)+2*(il-1)
                                if (nbval .eq. nbpf) then
                                    zr(ivale+il-1) = zr(ivale+il-1)+(scal11*scal22+scal12*sc&
                                                     &al21)*zr(iinte+idec)
                                else
                                    zr(ivale+2*(il-1)) = zr( &
                                                         ivale+2*(il-1))+(scal11*scal22+sca&
                                                         &l12*scal21)*zr(iinte+idec &
                                                         )
                                    zr(ivale+2*(il-1)+1) = zr( &
                                                           ivale+2*(il-1)+1)+(scal11*scal22-&
                                                           & scal12*scal21)*zr(iinte+idec+1 &
                                                           )
                                end if
                            end do
                        end do
                    end if
!
                end do
!
            else
!
                coedim = 0.25d0*coef1*coef1*rhof*rhof*uabs*uabs*uabs*phie*phie*phie
                do il = 1, nbpf
                    do iex2 = 1, nbexcp
                        scal12 = zr(iscal+nbexcp*(im1b-1)+iex2-1)
                        scal22 = zr(iscal+nbexcp*(im2b-1)+iex2-1)
                        iex = iex2*(iex2+1)/2
                        idec = 2*nbpf*(iex-1)+2*(il-1)
                        if (nbval .eq. nbpf) then
                            zr(ivale+il-1) = zr(ivale+il-1)+coedim*scal12*scal22*zr(iinte+i&
                                             &dec)
                        else
                            zr(ivale+2*(il-1)) = zr( &
                                                 ivale+2*(il-1) &
                                                 )+coedim*scal12*scal22*zr(iinte+idec &
                                                                           )
                        end if
                    end do
                end do
!
            end if
        end do
    end do
!
    call jedetr('&&SPECEP.TEMP.MAIL')
    call jedetr('&&SPECEP.TEMP.SCAL')
    call jedetr('&&SPECEP.TEMP.SWR ')
    call jedetr('&&SPECEP.TEMP.INTE')
!
    call jedema()
end subroutine
