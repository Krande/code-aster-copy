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
subroutine modeau(melflu, noma, geom, fsvr, base, &
                  freqi, nbm, nuor, vicoq, torco, &
                  tcoef, amor, masg, fact, amfr, &
                  vecpr, maj)
    implicit none
!  CONFIGURATION DE TYPE "COQUES CYLINDRIQUES COAXIALES"
!  CALCUL DES MODES EN EAU AU REPOS
!  APPELANT : FLUST4
!-----------------------------------------------------------------------
!  IN : MELFLU : NOM DU CONCEPT DE TYPE MELASFLU PRODUIT
!  IN : NOMA   : NOM DU CONCEPT DE TYPE MAILLAGE
!  IN : GEOM   : VECTEUR DE GRANDEURS GEOMETRIQUES CARACTERISTIQUES
!  IN : FSVR   : OBJET .FSVR DU CONCEPT TYPE_FLUI_STRU
!  IN : BASE   : NOM DU CONCEPT DE TYPE MODE_MECA DEFINISSANT LA BASE
!                MODALE DU SYSTEME AVANT PRISE EN COMPTE DU COUPLAGE
!  IN : FREQI  : FREQUENCES MODALES AVANT PRISE EN COMPTE DU COUPLAGE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : NUOR   : LISTE DES NUMEROS D'ORDRE DES MODES SELECTIONNES POUR
!                LE COUPLAGE (PRIS DANS LE CONCEPT MODE_MECA)
!  IN : VICOQ  : VECTEUR D'INDICES CARACTERISANT LE MOUVEMENT DES COQUES
!                POUR CHAQUE MODE PRIS EN COMPTE POUR LE COUPLAGE
!                VICOQ(IMOD)=1 COQUE INTERNE SEULE EN MVT
!                VICOQ(IMOD)=2 COQUE EXTERNE SEULE EN MVT
!                VICOQ(IMOD)=3 COQUES INTERNE + EXTERNE EN MVT
!  IN : TORCO  : TABLEAU DES ORDRES DE COQUE ET DEPHASAGES
!  IN : TCOEF  : TABLEAU DES COEFFICIENTS DES DEFORMEES AXIALES DES
!                MODES PRIS EN COMPTE
!  IN : AMOR   : LISTE DES AMORTISSEMENTS REDUITS MODAUX INITIAUX
!  OUT: MASG   : MASSES GENERALISEES DES MODES PERTURBES
!                = MASSES MODALES EN EAU AU REPOS
!  OUT: AMFR   : AMORTISSEMENTS MODAUX ET FREQUENCES EN EAU AU REPOS
!  OUT: VECPR  : DEFORMEES MODALES EN EAU AU REPOS, DECOMPOSEES SUR
!                LA BASE MODALE DU SYSTEME AVANT PRISE EN COMPTE DU
!                COUPLAGE
!  OUT: MAJ    : MASSES MODALES AJOUTEES PAR LE FLUIDE (DANS LA BASE EN
!                EAU AU REPOS)
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/bmocca.h"
#include "asterfort/dismoi.h"
#include "asterfort/extmod.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pmavec.h"
#include "asterfort/prmama.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/vphqrp.h"
#include "asterfort/wkvect.h"
#include "blas/ddot.h"
    character(len=19) :: melflu
    character(len=8) :: noma
    real(kind=8) :: geom(9), fsvr(7)
    character(len=8) :: base
    real(kind=8) :: freqi(*)
    integer(kind=8) :: nbm, nuor(nbm), vicoq(nbm)
    real(kind=8) :: torco(4, nbm), tcoef(10, nbm), amor(nbm), fact(nbm)
    real(kind=8) :: masg(nbm), amfr(nbm, 2), vecpr(nbm, nbm), maj(nbm)
!
    integer(kind=8) :: iddl(6)
    real(kind=8) :: mcf0, mi, mk, ki
    character(len=8) :: k8b
    character(len=14) :: numddl
    character(len=24) :: matria, nomcha
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icalc, idec1, idec2, idpla, idple, ier
    integer(kind=8) :: ifact, ifr, imat1, imat2, imata, imatm, imatz
    integer(kind=8) :: imax, imod, ivale, ivapr
    integer(kind=8) :: ivec, ivecw, iwrk2, jmod, k, kmod, lfacx
    integer(kind=8) :: lmasg, nbm2, nbnoe, neq, nitqr, numod
    real(kind=8) :: cf0, ck, fi, fim, fk, fre, omegai
    real(kind=8) :: pi, rmax, rtamp, s0, tole, u0
    blas_int :: b_incx, b_incy, b_n
!
!-----------------------------------------------------------------------
    data iddl/1, 2, 3, 4, 5, 6/
!
!-----------------------------------------------------------------------
    call jemarq()
!
!
!-----1.INITIALISATIONS ET CREATION DE VECTEURS DE TRAVAIL
!
    pi = r8pi()
    tole = 1.d-8
!
    u0 = 0.d0
    cf0 = 0.d0
    mcf0 = 1.d0
    s0 = 0.d0
!
    call wkvect('&&MODEAU.TEMP.MAT1', 'V V C', nbm*nbm, imat1)
    call wkvect('&&MODEAU.TEMP.MAT2', 'V V R', nbm*nbm, imat2)
    call wkvect('&&MODEAU.TEMP.MATM', 'V V R', nbm*nbm, imatm)
!
    call wkvect('&&MODEAU.TEMP.MATA', 'V V R', nbm*nbm, imata)
    call wkvect('&&MODEAU.TEMP.VAPR', 'V V R', nbm, ivapr)
    nbm2 = 2*nbm
    call wkvect('&&MODEAU.TEMP.VECW', 'V V R', nbm2, ivecw)
    call wkvect('&&MODEAU.TEMP.MATZ', 'V V R', nbm2*nbm2, imatz)
    call wkvect('&&MODEAU.TEMP.WRK2', 'V V R', nbm2, iwrk2)
!
    call wkvect('&&MODEAU.TEMP.VEC ', 'V V R', nbm, ivec)
    call wkvect('&&MODEAU.TEMP.FACT', 'V V R', 3*nbm, ifact)
!
!
!-----2.CALCUL DE LA MATRICE DE MASSE AJOUTEE  => MAT1 COMPLEXE
!       NB : LA MATRICE CALCULEE CORRESPOND A -MAJ
!
    call bmocca(u0, geom, cf0, mcf0, fsvr, &
                nbm, vicoq, torco, tcoef, s0, &
                s0, zc(imat1))
!
!
!-----3.ON SYMETRISE LA MATRICE DE MASSE AJOUTEE  => MAT2 REELLE
!       DEDOUBLEMENT DE MAT2 DANS MATM
!
    do jmod = 1, nbm
        do imod = 1, jmod
            idec1 = nbm*(jmod-1)+imod-1
            idec2 = nbm*(imod-1)+jmod-1
            zr(imat2+idec1) = -0.5d0*(dble(zc(imat1+idec1))+dble(zc(imat1+idec2)))
            zr(imatm+idec1) = zr(imat2+idec1)
        end do
    end do
    if (nbm .gt. 1) then
        do jmod = 1, nbm-1
            do imod = jmod+1, nbm
                idec1 = nbm*(jmod-1)+imod-1
                idec2 = nbm*(imod-1)+jmod-1
                zr(imat2+idec1) = zr(imat2+idec2)
                zr(imatm+idec1) = zr(imat2+idec1)
            end do
        end do
    end if
!
500 format('MAJ(', i3, ',', i3, ') = ', g23.16)
    ifr = iunifi('RESULTAT')
    write (ifr, *) '<MODEAU>'
    write (ifr, *)
    write (ifr, *) 'CALCUL DE LA MATRICE DE MASSES AJOUTEES DANS LA ',&
     &             'BASE MODALE EN AIR'
    write (ifr, *)
    do jmod = 1, nbm
        do imod = 1, nbm
            write (ifr, 500) imod, jmod, zr(imat2+nbm*(jmod-1)+imod-1)
        end do
        write (ifr, *)
    end do
!
!
!-----4.CALCUL DE LA MATRICE DE MASSE COMPLETE  => MATM
!       SIMULTANEMENT ON CALCULE LA MATRICE DU PROBLEME MODAL
!       GENERALISE EN EAU AU REPOS  => MATA
!
    do imod = 1, nbm
        numod = nuor(imod)
        call rsadpa(base, 'L', 1, 'FACT_PARTICI_DX', numod, &
                    0, sjv=lfacx, styp=k8b)
        zr(ifact+imod-1) = zr(lfacx)
        zr(ifact+nbm+imod-1) = zr(lfacx+1)
        zr(ifact+2*nbm+imod-1) = zr(lfacx+2)
        fi = freqi(numod)
        call rsadpa(base, 'L', 1, 'MASS_GENE', numod, &
                    0, sjv=lmasg, styp=k8b)
        mi = zr(lmasg)
        ki = 4.d0*pi*pi*fi*fi*mi
        zr(imatm+nbm*(imod-1)+imod-1) = zr(imatm+nbm*(imod-1)+imod-1)+mi
        do jmod = 1, nbm
            zr(imata+nbm*(jmod-1)+imod-1) = -zr(imatm+nbm*(jmod-1)+imod-1)/ki
        end do
    end do
!
!
!-----5.RESOLUTION DU PROBLEME MODAL GENERALISE EN EAU AU REPOS
!
    icalc = 1
    call vphqrp(zr(imata), nbm, nbm, icalc, zr(ivecw), &
                zr(imatz), nbm, zr(iwrk2), 30, ier, &
                nitqr)
    if (ier .ne. 0) then
        call utmess('F', 'ALGELINE_99')
    end if
!
    do imod = 1, nbm
        fre = dble(abs(zr(ivecw+2*(imod-1))))
        fim = dble(abs(zr(ivecw+2*(imod-1)+1)))
        if (fim .gt. (tole*fre)) then
            call utmess('F', 'ALGELINE2_2')
        end if
        zr(ivapr+imod-1) = zr(ivecw+2*(imod-1))
        do jmod = 1, nbm
            vecpr(imod, jmod) = zr(imatz+2*nbm*(jmod-1)+2*(imod-1))
        end do
    end do
!
!
!-----6.ON REORDONNE LES VALEURS PROPRES PAR VALEURS ABSOLUES
!       DECROISSANTES (VALEURS PROPRES CALCULEES = -1/OMEGA2)
!       SIMULTANEMENT ON EFFECTUE LES PERMUTATIONS DES COLONNES DE LA
!       MATRICE DES VECTEURS PROPRES
!
    if (nbm .gt. 1) then
        do imod = 1, nbm-1
            rmax = dble(abs(zr(ivapr+imod-1)))
            imax = imod
            do jmod = imod+1, nbm
                if (dble(abs(zr(ivapr+jmod-1))) .gt. rmax) then
                    rmax = dble(abs(zr(ivapr+jmod-1)))
                    imax = jmod
                end if
            end do
            zr(ivapr+imax-1) = zr(ivapr+imod-1)
            zr(ivapr+imod-1) = rmax
            if (rmax .eq. 0.d0) then
                call utmess('F', 'ALGELINE2_3')
            end if
            do kmod = 1, nbm
                rtamp = vecpr(kmod, imax)
                vecpr(kmod, imax) = vecpr(kmod, imod)
                vecpr(kmod, imod) = rtamp
            end do
        end do
    end if
    zr(ivapr+nbm-1) = dble(abs(zr(ivapr+nbm-1)))
!
!
!-----7.DECOMPOSITION DES DEFORMEES MODALES EN EAU AU REPOS SUR LA
!       BASE PHYSIQUE
!
    call dismoi('REF_RIGI_PREM', base, 'RESU_DYNA', repk=matria)
    call dismoi('NOM_NUME_DDL', matria, 'MATR_ASSE', repk=numddl)
    call dismoi('NB_EQUA', matria, 'MATR_ASSE', repi=neq)
!
    call jelira(noma//'.COORDO    .VALE', 'LONMAX', nbnoe)
    nbnoe = nbnoe/3
    call wkvect('&&MODEAU.TEMP.DPLA', 'V V R', 6*nbnoe*nbm, idpla)
    call extmod(base, numddl, nuor, nbm, zr(idpla), &
                neq, nbnoe, iddl, 6)
!
    call wkvect('&&MODEAU.TEMP.DPLE', 'V V R', 6*nbnoe*nbm, idple)
    call prmama(1, zr(idpla), 6*nbnoe, 6*nbnoe, nbm, &
                vecpr, nbm, nbm, nbm, zr(idple), &
                6*nbnoe, 6*nbnoe, nbm, ier)
!
    nomcha(1:13) = melflu(1:8)//'.C01.'
    nomcha(17:24) = '001.VALE'
    do imod = 1, nbm
        numod = nuor(imod)
        write (nomcha(14:16), '(I3.3)') numod
        call jeveuo(nomcha, 'E', ivale)
        do k = 1, 6*nbnoe
            zr(ivale+k-1) = zr(idple+6*nbnoe*(imod-1)+k-1)
        end do
        call jelibe(nomcha)
    end do
!
!
!-----8.CALCULS SIMULTANES :
!        - DES FREQUENCES PROPRES EN EAU AU REPOS
!        - DES MASSES MODALES EN EAU AU REPOS
!        - DES MASSES MODALES AJOUTEES PAR LE FLUIDE (EN EAU AU REPOS)
!        - DES AMORTISSEMENTS MODAUX EN EAU AU REPOS
!
    do imod = 1, nbm
!-------FREQUENCES PROPRES
        omegai = 1.d0/dble(sqrt(zr(ivapr+imod-1)))
        amfr(imod, 2) = omegai/(2.d0*pi)
!-------MASSES MODALES
        call pmavec('ZERO', nbm, zr(imatm), vecpr(1, imod), zr(ivec))
        b_n = to_blas_int(nbm)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        masg(imod) = ddot(b_n, vecpr(1, imod), b_incx, zr(ivec), b_incy)
!-------FACTEURS DE PARTICIPATION
        b_n = to_blas_int(nbm)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        fact(3*(imod-1)+1) = ddot(b_n, zr(ivec), b_incx, zr(ifact), b_incy)
        b_n = to_blas_int(nbm)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        fact(3*(imod-1)+2) = ddot(b_n, zr(ivec), b_incx, zr(ifact+nbm), b_incy)
        b_n = to_blas_int(nbm)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        fact(3*(imod-1)+3) = ddot(b_n, zr(ivec), b_incx, zr(ifact+2*nbm), b_incy)
!-------MASSES MODALES AJOUTEES PAR LE FLUIDE
        call pmavec('ZERO', nbm, zr(imat2), vecpr(1, imod), zr(ivec))
        b_n = to_blas_int(nbm)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        maj(imod) = ddot(b_n, vecpr(1, imod), b_incx, zr(ivec), b_incy)
!-------AMORTISSEMENTS MODAUX
        amfr(imod, 1) = 0.d0
        do kmod = 1, nbm
            numod = nuor(kmod)
            fk = freqi(numod)
            call rsadpa(base, 'L', 1, 'MASS_GENE', numod, &
                        0, sjv=lmasg, styp=k8b)
            mk = zr(lmasg)
            ck = 4.d0*pi*fk*amor(kmod)*mk
            amfr(imod, 1) = amfr(imod, 1)+vecpr(kmod, imod)*ck*vecpr(kmod, imod)
        end do
    end do
!
! --- MENAGE
!
    call jedetr('&&MODEAU.TEMP.MAT1')
    call jedetr('&&MODEAU.TEMP.MAT2')
    call jedetr('&&MODEAU.TEMP.MATM')
    call jedetr('&&MODEAU.TEMP.MATA')
    call jedetr('&&MODEAU.TEMP.VAPR')
    call jedetr('&&MODEAU.TEMP.VECW')
    call jedetr('&&MODEAU.TEMP.MATZ')
    call jedetr('&&MODEAU.TEMP.WRK2')
    call jedetr('&&MODEAU.TEMP.VEC ')
    call jedetr('&&MODEAU.TEMP.FACT')
    call jedetr('&&MODEAU.TEMP.DPLA')
    call jedetr('&&MODEAU.TEMP.DPLE')
!
    call jedema()
!
end subroutine
