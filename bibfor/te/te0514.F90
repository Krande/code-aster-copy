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
subroutine te0514(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/loncar.h"
#include "asterfort/ltequa.h"
#include "asterfort/ndcent.h"
#include "asterfort/padist.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/xdecou.h"
#include "asterfort/xdecov.h"
#include "asterfort/xdecqu.h"
#include "asterfort/xdecqv.h"
#include "asterfort/xdivte.h"
#include "asterfort/xelrex.h"
#include "asterfort/xnormv.h"
#include "asterfort/xxmmvd.h"
!
    character(len=16) :: option, nomte
!
!.......................................................................
!
!       CALCUL DES DONNEES TOPOLOGIQUES CONCERNANT CONCERNANT
!       LA DECOUPE DES ELEMENTS POUR L'INTEGRATION AVEC X-FEM
!                   (VOIR BOOK IV 27/10/04)
!
!  OPTION : 'TOPOSE' (X-FEM TOPOLOGIE DES SOUS-ELEMENTS)
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!......................................................................
!
!
    character(len=8) :: elp, noma, typma, enr, enr2, elrese(3), typsma
    integer(kind=8) :: igeom, jlsn, ifisc, nfisc, ptmaxi, pmmaxi
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jpmilt, zintmx, nfimax
    integer(kind=8) :: iadzi, iazk24, nno, nnn, jfisco, nsemax, jout6
    integer(kind=8) :: ninter, nit, nse, nnose, ise2, ncomph, ncompp, ncompc
    integer(kind=8) :: npts, cnse(6, 10), i, j, k, it, npi, ipt, ise, in, ni, cpt
    integer(kind=8) :: ndim, ibid, ndime, iad, jtab(7), jtab2(2), vali(2)
    integer(kind=8) :: nnc, npm, nmilie, nmfis, nbar, ar(12, 3)
    integer(kind=8) :: iret, nfiss, ifiss, ncomb, ninmax, nmmax, nbars, ars(12, 3)
    integer(kind=8) :: a1, a2, b1, b2, ncompm, iexit(2), joncno, jonact(8)
    parameter(ptmaxi=6, zintmx=5, pmmaxi=17, nsemax=6, nfimax=10)
    real(kind=8) :: nmil(3, 7), txlsn(28), ainter(ptmaxi*zintmx), rainter(4)
    real(kind=8) :: newpt(3), p(3), lonref, pinter(3*ptmaxi), lsn(3)
    real(kind=8) :: pmilie(3*pmmaxi), heav(nsemax*nfimax), u(3), v(3), normal(3)
    real(kind=8) :: xg(3), cridist, xref(81), ff(27), ptref(3), norme
    integer(kind=8) :: fisco(2*nfimax), fisc(2*nfimax), coupe(nfimax), zxain, ai, nnos
    parameter(ninmax=44, nmmax=264)
!
    integer(kind=8), dimension(:), allocatable :: ndoubl, ndoub2, ndoub3
!
    parameter(cridist=1.d-9)
    aster_logical :: deja, ajn, cut, lconnec_ok, pre1, jonc, condition_joncno
!
    data elrese/'SEG3', 'TRIA6', 'TETRA10'/
!......................................................................
!   LES TABLEAUX FISC, FISCO, PMILIE, PINTER ONT ETE
!   ALLOUE DE FACON STATIQUE POUR OPTIMISER LE CPU (CAR LES APPELS A
!   WKVECT DANS LES TE SONT COUTEUX).
!
!   On the other hand ndoubl, ndoub2, ndoub3 are too big to be allocated on the stack.
!   Hence the use of allocate. AS_ALLOCATE **is volontarily not used** given the performance
!   overhead and the fact it is unnecessary here (temporary arrays used only here)
    allocate (ndoubl(ninmax*(2**nfimax)))
    allocate (ndoub2(ninmax*(2**nfimax)))
    allocate (ndoub3(nmmax*(2**nfimax)))
!
    ASSERT(option .eq. 'TOPOSE')
    pmilie(:) = 0.d0
    ainter(:) = 0.d0
    rainter(:) = 0.d0
!
    call elref1(elp)
    call elrefe_info(fami='RIGI', ndim=ndime, nno=nno, nnos=nnos)
!
    call tecael(iadzi, iazk24, noms=0)
    noma = zk24(iazk24) (1:8)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
!     ATTENTION, NE PAS CONFONDRE NDIM ET NDIME  !!
!     NDIM EST LA DIMENSION DU MAILLAGE
!     NDIME EST DIMENSION DE L'ELEMENT FINI
!     PAR EXEMPLE, POUR LES ELEMENT DE BORDS D'UN MAILLAGE 3D :
!     NDIME = 2 ALORS QUE NDIM = 3
!
!     RECUPERATION DES ENTREES / SORTIE
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PLEVSET', 'L', jlsn)
    call jevech('PPINTTO', 'E', jpintt)
    call jevech('PCNSETO', 'E', jcnset)
    call jevech('PHEAVTO', 'E', jheavt)
    call jevech('PLONCHA', 'E', jlonch)
!
    call teattr('S', 'XFEM', enr, ibid)
!
    ncompm = 0
    if ((ibid .eq. 0) .and. ltequa(elp, enr)) then
        call jevech('PPMILTO', 'E', jpmilt)
        call tecach('OOO', 'PPMILTO', 'E', iret, nval=2, &
                    itab=jtab2)
        ncompm = jtab2(2)/ndim
    end if
!
    call tecach('OOO', 'PHEAVTO', 'E', iret, nval=7, &
                itab=jtab)
    ncomph = jtab(2)
    nfiss = jtab(7)
!
!   recuperation de PAINTTO pour les elements principaux
    jout6 = 0
    if (ndim .eq. ndime) then
        call jevech('PAINTTO', 'E', jout6)
    end if
!
!   Le modèle est-il HM-XFEM
    call teattr('C', 'HYDR1', enr2, iret)
    pre1 = (enr2 .eq. '1' .or. enr2 .eq. '2')
    condition_joncno = (pre1 .and. (enr(1:4) .eq. 'XH2C' .or. enr(1:4) .eq. 'XH3C'))
    joncno = 1
    if (condition_joncno) then
        call jevech('PJONCNO', 'E', joncno)
        do i = 1, 20
            zi(joncno-1+i) = 0
        end do
    end if
!
    do i = 1, nfimax
        coupe(i) = 0
        fisco(2*i-1) = 0
        fisco(2*i) = 0
        fisc(2*i-1) = 0
        fisc(2*i) = 0
    end do
    if (nfiss .gt. 1) then
        call jevech('PFISCO', 'L', jfisco)
        do i = 1, 2*nfiss
            fisco(i) = zi(jfisco+i-1)
        end do
    end if
!
    call tecach('OOO', 'PPINTTO', 'E', iret, nval=2, &
                itab=jtab2)
    ncompp = jtab2(2)/ndim
    call tecach('OOO', 'PCNSETO', 'E', iret, nval=2, &
                itab=jtab2)
    ncompc = jtab2(2)
!
! calcul du nombre de noeuds centraux
!
! initilisation : cas sans noeuds centraux
    nnc = 0
    nmil(:, :) = 0.d0
    txlsn(:) = 0.d0
!
! si l'element est incomplet, on calcule les informations concernant les les noeuds centraux
    if ((elp .eq. 'QU8') .or. (elp .eq. 'H20') .or. (elp .eq. 'P13') .or. (elp .eq. 'P15')) then
        call ndcent(igeom, ndim, zr(jlsn), nfiss, nmil, &
                    txlsn, nnc)
    end if
!
    iexit(1:2) = 0
73  continue
!
    npi = 0
    zxain = xxmmvd('ZXAIN')
    npm = 0
    ajn = .false.
    nmfis = 0
!
!     CALCUL D'UNE LONGUEUR CARACTERISTIQUE DE L'ELEMENT
    call loncar(ndim, typma, zr(igeom), lonref)
!
!     ON SUBDIVISE L'ELEMENT PARENT EN NIT SOUS-ELEMENTS
    call xdivte(elp, zi(jcnset), nit, nnose, iexit)
    cpt = nit
!     PROBLEME DE DIMENSSIONNEMENT DANS LES CATALOGUES
    ASSERT(ncompc .ge. nnose*ncomph)
!
    vali(1) = nfimax
    vali(2) = nfiss
    if (nfiss .gt. nfimax) then
        call utmess('F', 'XFEM2_6', ni=2, vali=vali)
    end if
!
!     BOUCLE SUR LES FISSURES
!
    do ifiss = 1, nfiss
!
        do i = 1, ifiss
            fisc(i) = 0
        end do
        ifisc = ifiss
        nfisc = 0
80      continue
        if (fisco(2*ifisc-1) .gt. 0) then
            nfisc = nfisc+1
            fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
            ifisc = fisco(2*ifisc-1)
            fisc(2*(nfisc-1)+1) = ifisc
            goto 80
        end if
!
!       BOUCLE SUR LES NIT TETRAS
        do it = 1, nit
!
!         DECOUPAGE EN NSE SOUS-ELEMENTS
            pinter(:) = 0.d0
            nmilie = 0
            ninter = 0
            npts = 0
            do i = 1, nsemax*nfimax
                heav(i) = 0.d0
            end do
!
            if (.not. iselli(elp)) then
                call xdecqu(nnose, it, ndim, zi(jcnset), jlsn, &
                            igeom, pinter, ninter, npts, ainter, &
                            pmilie, nmilie, nmfis, nmil, txlsn, &
                            zr(jpintt), zr(jpmilt), ifiss, nfiss, fisc, &
                            nfisc, cut, coupe, iexit, joncno, &
                            condition_joncno)
                if (iexit(1) .eq. 1) goto 73
!
                call xdecqv(nnose, it, zi(jcnset), zi(jheavt), zr(jlsn), &
                            igeom, ninter, npts, ndim, ainter, &
                            nse, cnse, heav, nsemax, pinter, &
                            pmilie, zr(jpintt), zr(jpmilt), cut, ncomph, &
                            nfisc, nfiss, ifiss, elp, fisc, &
                            lonref, txlsn, nmil)
            else
                call xdecou(ndim, elp, nno, nnose, it, &
                            zr(jpintt), zi(jcnset), zr(jlsn), fisc, igeom, &
                            nfiss, ifiss, pinter, ninter, npts, &
                            ainter, lonref, nfisc)
                call xdecov(ndim, elp, nno, nnose, it, &
                            zr(jpintt), zi(jcnset), zi(jheavt), ncomph, zr(jlsn), &
                            fisc, igeom, nfiss, ifiss, pinter, &
                            ninter, npts, ainter, nse, cnse, &
                            heav, nfisc, nsemax)
            end if
!
! ----- - BOUCLE SUR LES NINTER POINTS D'INTER : ARCHIVAGE DE PINTTO
            do ipt = 1, ninter
                do j = 1, ndim
                    newpt(j) = pinter(ndim*(ipt-1)+j)
                end do
                do i = 1, 4
                    rainter(i) = ainter(zxain*(ipt-1)+i)
                end do
!
!           VERIF SI EXISTE DEJA DANS PINTTO
                deja = .false.
                do i = 1, npi
                    do j = 1, ndim
                        p(j) = zr(jpintt-1+ndim*(i-1)+j)
                    end do
                    if (padist(ndim, p, newpt) .lt. (lonref*cridist)) then
                        deja = .true.
                        ni = i
                    end if
                end do
                ai = nint(ainter(zxain*(ipt-1)+1))
!             SI LE POINT D'INTER COINCIDE AVEC UN NOEUD
                if (ai .eq. 0) cycle
                if (.not. deja) then
                    npi = npi+1
                    ni = npi
!             NOMBRE TOTAL DE PT D'INTER LIMITE A LA TAILLE DE LA CARTE
                    if (npi .gt. ncompp) then
                        call utmess('F', 'XFEM_55')
                    end if
!             ARCHIVAGE DE PINTTO
                    do j = 1, ndim
                        zr(jpintt-1+ndim*(npi-1)+j) = newpt(j)
                    end do
!
!             ARCHIVAGE DE PAINTTO POUR LES ELEMENTS PRINCIPAUX
                    if (ndim .eq. ndime) then
                        do k = 1, zxain
                            zr(jout6-1+zxain*(npi-1)+k) = 0.d0
                        end do
                        if (ifiss .gt. 1) then
!             MARQUAGE POINT DE JONCTION DE FISSURES
                            if (abs(rainter(4)+1.d0) .le. r8prem()) then
                                zr(jout6-1+zxain*(npi-1)+4) = -1.d0
                            end if
                            call conare(typma, ar, nbar)
                            call xelrex(elp, nno, xref)
                            call reeref(elp, nno, zr(igeom), newpt, ndim, &
                                        ptref, ff)
                            u(:) = 0.d0
                            v(:) = 0.d0
                            normal(:) = 0.d0
                            do j = 1, nbar
                                do k = 1, ndim
                                    u(k) = xref((ar(j, 1)-1)*ndim+k)-xref((ar(j, 2)-1)*ndim+k)
                                    v(k) = ptref(k)-xref((ar(j, 2)-1)*ndim+k)
                                end do
                                call provec(u, v, normal)
                                call xnormv(3, normal, norme)
                                if (norme .lt. cridist) then
                                    zr(jout6-1+zxain*(npi-1)+1) = j
                                    zr(jout6-1+zxain*(npi-1)+3) = rainter(3)
                                    zr(jout6-1+zxain*(npi-1)+4) = rainter(4)
                                end if
                            end do
                        else
                            call conare(typma, ar, nbar)
                            typsma = elrese(ndim)
                            call conare(typsma, ars, nbars)
                            b1 = zi(jcnset+nnose*(it-1)+ars(nint(rainter(1)), 1)-1)
                            b2 = zi(jcnset+nnose*(it-1)+ars(nint(rainter(1)), 2)-1)
!
                            do j = 1, nbar
                                a1 = ar(j, 1)
                                a2 = ar(j, 2)
                                if (a1 .eq. b1 .and. a2 .eq. b2) then
                                    zr(jout6-1+zxain*(npi-1)+1) = j
                                    zr(jout6-1+zxain*(npi-1)+3) = rainter(3)
                                    zr(jout6-1+zxain*(npi-1)+4) = rainter(4)
                                else if (a1 .eq. b2 .and. a2 .eq. b1) then
                                    zr(jout6-1+zxain*(npi-1)+1) = j
                                    zr(jout6-1+zxain*(npi-1)+3) = rainter(3)
                                    zr(jout6-1+zxain*(npi-1)+4) = 1.d0-rainter(4)
                                end if
                            end do
                        end if
                    end if
!
!             MISE A JOUR DU CNSE (TRANSFORMATION DES 100 EN 1000...)
                    do ise = 1, nse
                        do in = 1, ndime+1
                            if (cnse(ise, in) .eq. 100+ipt) cnse(ise, in) = 1000+ni
                        end do
                    end do
                else
                    do ise = 1, nse
                        do in = 1, ndime+1
                            if (cnse(ise, in) .eq. 100+ipt) cnse(ise, in) = 1000+ni
                        end do
                    end do
                end if
            end do
!
!
! ------- BOUCLE SUR LES NMILIE POINTS MILIEUX : ARCHIVAGE DE PMILTO
            if (.not. iselli(elp)) then
                do ipt = 1, nmilie
                    do j = 1, ndim
                        newpt(j) = pmilie(ndim*(ipt-1)+j)
                    end do
!             VERIF SI EXISTE DEJA DANS PMILTO
                    deja = .false.
                    do i = 1, npm
                        do j = 1, ndim
                            p(j) = zr(jpmilt-1+ndim*(i-1)+j)
                        end do
                        if (padist(ndim, p, newpt) .lt. (lonref*cridist)) then
                            deja = .true.
                            ni = i
                        end if
                    end do
                    if (.not. deja) then
                        npm = npm+1
!               NOMBRE TOTAL DE POINTS MILIEUX LIMITE A PMMAX
                        if (npm .gt. ncompm) then
                            call utmess('F', 'XFEM_55')
                        end if
!               ARCHIVAGE DE PMILTO
                        do j = 1, ndim
                            zr(jpmilt-1+ndim*(npm-1)+j) = pmilie(ndim*(ipt-1)+j)
                        end do
!
!               MISE A JOUR DU CNSE (TRANSFORMATION DES 200 EN 2000...)
                        do ise = 1, nse
                            do in = 1, nnose
                                if (cnse(ise, in) .eq. 200+ipt) cnse(ise, in) = 2000+npm
                            end do
                        end do
                    else
                        do ise = 1, nse
                            do in = 1, nnose
                                if (cnse(ise, in) .eq. 200+ipt) cnse(ise, in) = 2000+ni
                            end do
                        end do
                    end if
!
                end do
            end if
! ------- BOUCLE SUR LES NSE SOUS-ELE : ARCHIVAGE DE PCNSETO, PHEAVTO
            do ise = 1, nse
                if (ise .eq. 1) then
                    ise2 = it
                else
                    cpt = cpt+1
                    ise2 = cpt
                end if
!         NOMBRE TOTAL DE SOUS-ELEMENTS LIMITE A LA TAILLE DE LA CARTE
!           ARCHIVAGE DE PHEAVTO
                do i = 1, ifiss
                    zi(jheavt-1+ncomph*(i-1)+ise2) = nint(heav(ifiss*( &
                                                               ise-1)+i))
                end do
!           ARCHIVAGE DE PCNSETO
                lconnec_ok = .true.
                do in = 1, nnose
                    lconnec_ok = lconnec_ok .and. &
                                 ((cnse(ise, in) .gt. 0 .and. cnse(ise, in) .le. nno+nnc) .or. &
                                  (cnse(ise, in) .gt. 1000 .and. cnse(ise, in) .le. 1000+ncompp) &
                                  .or. &
                                  (cnse(ise, in) .gt. 2000 .and. cnse(ise, in) .le. 2000+ncompm))
                    zi(jcnset-1+nnose*(ise2-1)+in) = cnse(ise, in)
                end do
                if (.not. lconnec_ok) then
                    call utmess('F', 'XFEMPRECOND_8')
                end if
            end do
!
        end do
        nit = cpt
!
    end do
!
!     ARCHIVAGE DE LONCHAM SOUS ELEMENTS
    if (nit*nnose .gt. ncompc) then
        call utmess('F', 'XFEM_55')
    end if
    zi(jlonch-1+1) = nit
!
!     REMPLISSAGE DE PJONCNO POUR LES ELEMENTS QUI CONTIENNENT LA JONCTION
    jonc = .false.
    do i = 1, 8
        jonact(i) = 0
    end do
    if (condition_joncno) then
        do i = 1, nnos
            if (zi(joncno-1+i) .eq. 1) jonc = .true.
        end do
        if (jonc) then
            call conare(typma, ar, nbar)
            do i = 1, nbar
                lsn(1) = zr(jlsn-1+(ar(i, 1)-1)*nfiss+1)
                lsn(2) = zr(jlsn-1+(ar(i, 2)-1)*nfiss+1)
                lsn(3) = zr(jlsn-1+(ar(i, 3)-1)*nfiss+1)
                if (abs(lsn(1)) .le. r8prem()) then
                    jonact(ar(i, 1)) = 1
                else if (abs(lsn(2)) .le. r8prem()) then
                    jonact(ar(i, 2)) = 1
                else if (abs(lsn(3)) .le. r8prem()) then
                    jonact(ar(i, 1)) = 1
                    jonact(ar(i, 2)) = 1
                else if (lsn(1)*lsn(2) .lt. 0.d0) then
                    jonact(ar(i, 1)) = 1
                    jonact(ar(i, 2)) = 1
                end if
            end do
            do i = 1, nnos
                if (jonact(i) .eq. 0) zi(joncno-1+i) = 0
            end do
            do i = 1, nnos
                if (zi(joncno-1+i) .ge. 1) zi(joncno-1+i) = jonact(i)
            end do
        end if
    end if
!
!     ARCHIVAGE DE LONCHAM POINTS D'INTERSECTION
    zi(jlonch-1+2) = npi
!
!     SERT UNIQUEMENT POUR LA VISU
!
!     POUR COMPTER LES NOEUDS EN DOUBLE (EN POST-TRAITEMENT)
!     NCOMB COMBINAISONS POSSIBLES
    ncomb = 2**nfiss
!     ON VA JUSQU'A NCOMPP, NOMBRE MAX DE POINTS D'INTERSECTION
    ASSERT(ncompp .le. ninmax)
    do i = 1, ninmax*(2**nfimax)
        ndoub2(i) = 0
    end do
!     POUR COMPTER LES PT D'INTER EN DOUBLE (EN POST-TRAITEMENT)
!     ON VA JUSQU'A NNO+1 POUR PRENDRE EN COMPTE LE 9EME NOEUD DU QUAD8
    ASSERT((nno+nnc) .le. ninmax)
    do i = 1, ninmax*(2**nfimax)
        ndoubl(i) = 0
    end do
!     POUR COMPTER LES PT MILIEU EN DOUBLE (EN POST-TRAITEMENT)
!     ON VA JUSQU'A NNO+NNC POUR PRENDRE EN COMPTE LES NOEUDS CENTRAUX
    do i = 1, nmmax*(2**nfimax)
        ndoub3(i) = 0
    end do
!
    do ise = 1, nit
        do in = 1, nnose
            i = zi(jcnset-1+nnose*(ise-1)+in)
            if (i .lt. 1000) then
! ----- ON SE PLACE EN BASE 2, LA DIMENSSION DU PB EST NFISS
! ----- POUR CHAQUE FISSURE H=-1 OU 0=>0
!                           H=+1     =>1
                iad = ncomb*(i-1)
                do ifiss = 1, nfiss
                    if (zi(jheavt-1+ncomph*(ifiss-1)+ise) .eq. 1) then
                        iad = iad+2**(nfiss-ifiss)
                    end if
                end do
                ndoubl(iad+1) = 1
            else if (i .gt. 1000 .and. i .lt. 2000) then
                iad = ncomb*(i-1001)
                do ifiss = 1, nfiss
                    if (zi(jheavt-1+ncomph*(ifiss-1)+ise) .eq. 1) then
                        iad = iad+2**(nfiss-ifiss)
                    end if
                end do
                ndoub2(iad+1) = 1
            else if (i .gt. 2000 .and. i .lt. 3000) then
                iad = ncomb*(i-2001)
                do ifiss = 1, nfiss
                    if (zi(jheavt-1+ncomph*(ifiss-1)+ise) .eq. 1) then
                        iad = iad+2**(nfiss-ifiss)
                    end if
                end do
                ndoub3(iad+1) = 1
            end if
        end do
    end do
!
!     NOMBRE DE NOUVEAUX POINTS : NNN
    nnn = 0
!  -  ON AJOUTE LES PT D'INTER QUI NE SONT PAS DES NOEUDS
!  -  AUTANT DE FOIS QU'IL Y A DE COMBINAISON D'HEAVISIDE DIFFERENTES
    do i = 1, ncomb*ncompp
        nnn = nnn+ndoub2(i)
    end do
!
!  -  ON AJOUTE LES NOEUDS, AUTANT DE FOIS QU'IL Y A DE COMBINAISON
!  -  D'HEAVISIDE DIFFERENTES.
!     ON VA JUSQU'A NNO+NNC POUR PRENDRE EN COMPTE LES NOEUDS CENTRAUX
    do i = 1, ncomb*(nno+nnc)
        nnn = nnn+ndoubl(i)
    end do
!
    if (.not. iselli(elp)) then
!  - CAS QUADRATIQUE, ON AJOUTE LES NOUVEAUX NOEUDS MILIEUX DES SE
!  -  AUTANT DE FOIS QU'IL Y A DE COMBINAISON D'HEAVISIDE DIFFERENTES
!        nnn=nnn+nmfis+npm
        do i = 1, ncomb*nmmax
            nnn = nnn+ndoub3(i)
        end do
    end if
!
    zi(jlonch-1+3) = nnn
    zi(jlonch-1+4) = 0
!   stockage des noeuds centraux servant a definir la connectivite sur sous-tetra courant
    xg(:) = 0.d0
!   pour chaque noeud central
    do j = 1, nnc
        deja = .false.
!     recuperation des coordonnees de refrence du noeud central courant
        do k = 1, ndim
            xg(k) = nmil(k, j)
        end do
!     pour chaque noeud de chaque sous-tetra
        do i = 1, nit*nnose
!       si le noeud courant du sous-tetra courant est le noeud central courant
            if (zi(jcnset-1+i) .eq. (nno+j)) then
!         s'il n'a pas deja ete ajoute, on ajoute le noeud central courant
                if (.not. deja) then
                    deja = .true.
                    npm = npm+1
                    do k = 1, ndim
                        zr(jpmilt+(npm-1)*ndim+k-1) = xg(k)
                    end do
                end if
!         NOMBRE TOTAL DE POINTS MILIEUX LIMITE A PMMAX
                if (npm .gt. ncompm) then
                    call utmess('F', 'XFEM_55')
                end if
!         stockage du noeud central courant
                zi(jcnset-1+i) = 3000+npm
            end if
        end do
    end do
!
!     ARCHIVAGE DE LONCHAM POINTS MILIEUX
!
    if (.not. iselli(elp)) then
        zi(jlonch-1+4) = npm
    end if
!
!   Deallocation of temporary arrays on heap
    deallocate (ndoubl)
    deallocate (ndoub2)
    deallocate (ndoub3)
!
end subroutine
