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
subroutine op0131()
    implicit none
!     CALCUL DE REPONSE DYNAMIQUE SOUS FORME D INTERSPECTRE
!     EXCITATION ET REPONSES SONT DES INTERSPECTRES.
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8depi.h"
#include "asterfort/caldis.h"
#include "asterfort/calpim.h"
#include "asterfort/infmaj.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mlmatc.h"
#include "asterfort/preflx.h"
#include "asterfort/prekpr.h"
#include "asterfort/reciex.h"
#include "asterfort/recire.h"
#include "asterfort/recmod.h"
#include "asterfort/recmst.h"
#include "asterfort/rslipa.h"
#include "asterfort/titre.h"
#include "asterfort/transf.h"
#include "asterfort/vriale.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iderex, iderre, iret, ilnoex, ilcpex, i1, ilvaex, napexc, nindex
    integer(kind=8) :: nnoeex, ncmpex, nvasex, i2, nbptmd, iadfrq, ilamor, ij1, igim
    integer(kind=8) :: igre, ilmode, ilamsc, nmost1, iadsc3, n1, iadhii, imoddy
    integer(kind=8) :: ifreq2, ifreq1, nmost3, ndimre, ni, itail1, iexp, isign, iadj, j
    integer(kind=8) :: iadjs, iadg, iadjg, iadjgj, nj, ni1, ni2, i, nk, ij, nbmode
    integer(kind=8) :: nbddl, nbamor, npdsc3
    integer(kind=8) :: ispec, lfreq, lnumi, lnumj, lrefe, nbabs, num, ilfex2
    real(kind=8) :: bande(2), depi, fremin, fremax, pas, r8amor, r8freq, r8bid1
    real(kind=8) :: r8sign, r8omeg, r8omg2, pim, frefin
    complex(kind=8) :: xcj, xcgrep, xcg, xch
    character(len=4) :: typopt, excmod, frexci
    character(len=8) :: mtrmas, modmec, intexc, modsta, mailla, numer
    character(len=8) :: intrep, nomref
    character(len=8) :: chamat, celem, tymmec
    character(len=16) :: tyconc, nomcmd, graexc, grdmod, nocham
    character(len=24) :: lifex2, chvale, chfreq, chnumi, chnumj
    real(kind=8), pointer :: vpim(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    depi = r8depi()
!
!---0.---VERIFICATIONS
!
    call vriale()
!
    call infmaj()
!
!---1-----RECUPERATION DES ARGUMENTS DE LA COMMANDE
!---1.1---NOM DU RESULTAT, TYPE DU RESULTAT, ET NOM DE LA COMMANDE
!
    call getres(intrep, tyconc, nomcmd)
!
    nomref = intrep(1:8)
!
!---1.2---INTERSPECTRE EXCIT
!
    call reciex(intexc, iderex, nindex, nnoeex, ncmpex, &
                nvasex, graexc, excmod, napexc)
!
    ilnoex = 1
    call jeexin('&&OP0131.LISTENOEEXC', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LISTENOEEXC', 'L', ilnoex)
!
    ilcpex = 1
    call jeexin('&&OP0131.LISTECMPEXC', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LISTECMPEXC', 'L', ilcpex)
!
    ilvaex = 1
    call jeexin('&&OP0131.LVECTASSEXC', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LVECTASSEXC', 'L', ilvaex)
!
!---1.3---INTERSPECTRE REPONSE
!
    call recire(typopt, iderre, frexci, fremin, fremax, &
                pas, nbptmd)
!
!---1.4---MODES DYNAMIQUES
!
    call recmod(modmec, nbmode, nbamor, bande, tymmec, &
                grdmod)
    call rslipa(modmec, 'FREQ', '&&OP0131.LIFREQ', iadfrq, n1)
!
    ilamor = 1
    call jeexin('&&OP0131.LISTEAMOR', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LISTEAMOR', 'L', ilamor)
!
    ilmode = 1
    call jeexin('&&OP0131.LISTEMODES', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LISTEMODES', 'L', ilmode)
!
!---1.5---MODES STATIQUES
!
    call recmst(graexc, grdmod, nnoeex, ilnoex, ilcpex, &
                nmost1, modsta)
    call jeexin('&&OP0131.LISTADRMODSTAC', iret)
    if (iret .gt. 0) call jeveuo('&&OP0131.LISTADRMODSTAC', 'L', ilamsc)
!
!---1.6---CONSTITUTION DE LA LISTE COMPLETE DES ET DDLS EN CAS D EFFORT
!         ET DES ADRESSES DES MATRICES ELEMENTAIRES CORRESPONDANTES
!
    call prekpr(modmec, mtrmas, nbddl, numer, mailla, &
                chamat, celem)
!
!---2.----CALCUL DES COEFFICIENTS DE PARTICIPATION GENERALISE
!        IL Y EN A NAPEXC*NBMODE
!
    call calpim(graexc, excmod, napexc, nbmode, tymmec, &
                mtrmas, numer, nbddl, zk8(ilnoex), zk8(ilcpex), &
                nvasex, zk8(ilvaex))
    call jeveuo('&&OP0131.PIM', 'L', vr=vpim)
!
!---3-----MISE AU POINT DE LA DISCRETISATION FREQUENTIELLE
!   4-----ET RECALCUL DES DSP EXCITS DANS LA DISCRETISATION REPONSE
!
    call caldis(fremax, fremin, pas, frexci, nbptmd, &
                nbmode, zi(ilmode), zr(iadfrq), zr(ilamor), nindex, &
                npdsc3, frefin)
!
    call jeveuo('&&OP0131.DISCR3', 'E', iadsc3)
    lifex2 = '&&OP0131.LIADRFOE.FRQRE'
!
!
!---5.----PRISE EN COMPTE D UN COEFFICIENT MULTIPLICATEUR POUR LES
!         SOURCES DE PRESSION OU DE FORCE
!
!     ON MULTIPLIE LE DEBIT MASSE PAR RO ET OMEGA  POUR DES SOURCES
!     DE DEBIT
!     ON CONDITIONNE L INTERSPECTRE POUR LES SOURCES DE PRESSION
!     DE FACON A CE QU IL SOIT ENTIEREMENT CORRELE AVEC PRISE EN COMPTE
!     DE LA DISTANCE ENTRE LES NOEUDS D APPLICATION DE LA SOURCE
!
    if (graexc(1:5) .eq. 'SOUR_') then
        call preflx(graexc, mailla, chamat, celem, npdsc3, &
                    iadsc3, nindex, ilnoex, lifex2)
    end if
    call jeveuo(lifex2, 'L', ilfex2)
!
!---6-----CALCUL DE LA REPONSE EN DSP
!---6.1---CALCUL ET STOCKAGE DE HII
!
    call wkvect('&&OP0131.HII', 'V V C', nbmode*npdsc3, iadhii)
    do imoddy = 1, nbmode
        r8amor = zr(ilamor+imoddy-1)
        ifreq2 = iadfrq-1+zi(ilmode+imoddy-1)
        do ifreq1 = 1, npdsc3
            r8freq = zr(iadsc3+ifreq1-1)
            r8bid1 = r8freq/zr(ifreq2)
            call transf(r8bid1, r8amor, xch)
            xch = xch/(depi*zr(ifreq2))**2
            zc(iadhii-1+npdsc3*(imoddy-1)+ifreq1) = xch
        end do
    end do
!
!---6.2---CREATION DE LA MATRICE DE RETOUR DANS L ESPACE PHYSIQUE
!
    if (graexc .eq. 'DEPL_R') then
        nmost3 = nnoeex
    else
        nmost3 = 0
    end if
!
!------STOCKAGE DES NUMEROS DES DDLS REPONSE
!
    ndimre = nbmode+nmost3
!
!
    if (iderre .eq. 0) then
        nocham = 'DEPL_GENE'
    else if (iderre .eq. 1) then
        nocham = 'VITE_GENE'
    else
        nocham = 'ACCE_GENE'
    end if
!
    call wkvect(nomref//'.REFE', 'G V K16', 3, lrefe)
    zk16(lrefe) = nocham
    zk16(lrefe+1) = typopt
    zk16(lrefe+2) = 'FREQ'
!
    itail1 = ndimre*(1+ndimre)/2
    chnumi = nomref//'.NUMI'
    call wkvect(chnumi, 'G V I', itail1, lnumi)
    chnumj = nomref//'.NUMJ'
    call wkvect(chnumj, 'G V I', itail1, lnumj)
    chvale = nomref//'.VALE'
    call jecrec(chvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                itail1)
    chfreq = nomref//'.DISC'
    call wkvect(chfreq, 'G V R', npdsc3, lfreq)
!
    ij1 = 0
    do i1 = 1, ndimre
        do i2 = i1, ndimre
            ij1 = ij1+1
            zi(lnumi-1+ij1) = i1
            zi(lnumj-1+ij1) = i2
            if ((typopt .eq. 'TOUT') .or. (i1 .eq. i2)) then
                if (i1 .eq. i2) then
                    nbabs = npdsc3
                else
                    nbabs = 2*npdsc3
                end if
            else
                nbabs = 6
            end if
            call jecroc(jexnum(chvale, ij1))
            call jeecra(jexnum(chvale, ij1), 'LONMAX', nbabs)
            call jeecra(jexnum(chvale, ij1), 'LONUTI', nbabs)
        end do
    end do
!
!---6.5---COEFFICIENT EN PUISSANCE DE W : IEXP
!
    iexp = iderre-iderex
    if (graexc .ne. 'DEPL_R') iexp = iexp-2
!
!---6.6---SIGNE DEDUIT DE L ORDRE DE DERIVATION
!
    isign = mod(iderre-iderex+6, 2)
    if (isign .eq. 1) then
        r8sign = -1.d0
    else
        r8sign = 1.d0
    end if
!
!---6.7---CREATION DES MATRICES DE TRAVAIL
!
    n1 = nmost3+nbmode
    call wkvect('&&OP0131.J', 'V V C', n1*napexc, iadj)
    call wkvect('&&OP0131.JS', 'V V C', n1*napexc, iadjs)
    call wkvect('&&OP0131.G', 'V V C', napexc*napexc, iadg)
    call wkvect('&&OP0131.JG', 'V V C', n1*napexc, iadjg)
    call wkvect('&&OP0131.JGJ', 'V V C', n1*n1, iadjgj)
!
!---6.8---BOUCLE DE CALCUL SUR CHAQUE FREQUENCE REPONSE
!
    nj = napexc
    ni1 = nmost3
    ni2 = nbmode
    ni = ni1+ni2
!
    do ifreq1 = 1, npdsc3
!
        r8freq = zr(iadsc3+ifreq1-1)
        zr(lfreq-1+ifreq1) = r8freq
        r8omeg = depi*r8freq
        r8omg2 = r8omeg*r8omeg
        if (r8omeg .eq. 0.d0) then
            if (iexp .eq. 0) then
                r8bid1 = 0.d0
            else
                r8bid1 = 1.d0
            end if
        else
            r8bid1 = r8sign*(r8omeg**iexp)*r8omg2
        end if
!
!---------REMPLISSAGE DE LA MATRICE J ET DE JS
!
        do j = 1, nj
!  MATRICES DES MODES STATIQUES EST L IDENTITE
            if (j .le. ni1) then
                xcj = dcmplx(1.d0, 0.d0)
                zc(iadj-1+ni*(j-1)+j) = xcj
                zc(iadjs-1+nj*(j-1)+j) = dconjg(xcj)
            end if
            do i = 1, ni2
                pim = vpim(1+(j-1)*ni2+i-1)
                xcj = r8bid1*pim*zc(iadhii-1+npdsc3*(i-1)+ifreq1)
                zc(iadj-1+ni*(j-1)+ni1+i) = xcj
                zc(iadjs-1+nj*(i+ni1-1)+j) = dconjg(xcj)
            end do
        end do
!
!---------STOCKAGE GEXC A LA FREQUENCE IFREQ1
!
        do i = 1, napexc
            do j = 1, napexc
                if (j .ge. i) then
                    ij1 = (j*(j-1))/2+i
                    igre = zi(ilfex2-1+ij1)+npdsc3-1+2*(ifreq1-1)+1
                    igim = zi(ilfex2-1+ij1)+npdsc3-1+2*(ifreq1-1)+2
                    xcg = dcmplx(zr(igre), zr(igim))
                else
                    ij1 = (i*(i-1))/2+j
                    igre = zi(ilfex2-1+ij1)+npdsc3-1+2*(ifreq1-1)+1
                    igim = zi(ilfex2-1+ij1)+npdsc3-1+2*(ifreq1-1)+2
                    xcg = dcmplx(zr(igre), -zr(igim))
                end if
                zc(iadg-1+(j-1)*napexc+i) = xcg
            end do
        end do
!
!---------CALCUL JG A LA FREQUENCE IFREQ1
!
        call mlmatc(n1, napexc, napexc, zc(iadj), zc(iadg), &
                    zc(iadjg))
!
!---------CALCUL JGJ* A LA FREQUENCE IFREQ1
!
        call mlmatc(n1, napexc, n1, zc(iadjg), zc(iadjs), &
                    zc(iadjgj))
!
!----------SI REPONSE MODALE ON SORT ICI LES RESULTATS
!
        do i = 1, ndimre
            if (typopt .eq. 'DIAG ') then
                nk = i
            else
                nk = ndimre
            end if
            do j = i, nk
                ij = 0
                do num = 1, itail1
                    if ((i .eq. zi(lnumi-1+num)) .and. (j .eq. zi(lnumj-1+num))) then
                        ij = num
                        goto 650
                    end if
                end do
650             continue
                xcgrep = zc(iadjgj-1+(j-1)*ndimre+i)
                call jeveuo(jexnum(chvale, ij), 'E', ispec)
                if (i .eq. j) then
                    zr(ispec-1+ifreq1) = dble(xcgrep)
                else
                    zr(ispec-1+2*(ifreq1-1)+1) = dble(xcgrep)
                    zr(ispec-1+2*(ifreq1-1)+2) = dimag(xcgrep)
                end if
            end do
        end do
!
    end do
!
    if (typopt .eq. 'DIAG') then
        ij = 0
        do i = 1, ndimre
            do j = i+1, ndimre
                ij = ij+1
                call jeveuo(jexnum(chvale, ij), 'E', ispec)
                zr(ispec) = 0.d0
                zr(ispec+1) = frefin
                zr(ispec+2) = 0.d0
                zr(ispec+3) = 0.d0
                zr(ispec+4) = 0.d0
                zr(ispec+5) = 0.d0
            end do
        end do
    end if
!
    call titre()
!
!
    call jedetr('&&OP0131.LIFREQ')
    call jedema()
end subroutine
