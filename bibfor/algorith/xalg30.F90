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
subroutine xalg30(ndim, elrefp, it, nnose, cnset, &
                  typma, ndime, geom, lsnelp, pmilie, &
                  ninter, ainter, ar, npts, nptm, &
                  pmmax, nmilie, mfis, lonref, pinref, &
                  pintt, pmitt, jonc, exit)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/detefa.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/xajpmi.h"
#include "asterfort/xmifis.h"
#include "asterfort/xmilar.h"
#include "asterfort/xmilfa.h"
#include "asterfort/xstudo.h"
#include "asterfort/xxmmvd.h"
    character(len=8) :: typma, elrefp
    integer(kind=8) :: ndim, ndime, it, nnose, cnset(*), exit(2)
    integer(kind=8) :: ninter, pmmax, npts, nptm, nmilie, mfis, ar(12, 3)
    real(kind=8) :: ainter(*), pmilie(*), lonref, lsnelp(*)
    real(kind=8) :: pinref(*), pintt(*), pmitt(*), geom(81)
    aster_logical :: jonc
!                    BUT :  TROUVER LES PTS MILIEUX DANS L ELEMENT COUPE
!
!     ENTREE
!       NDIM     : DIMENSION DE L ELEMENT
!       TYPMA    : TYPE DE MAILLE
!       TABCO    : COORDONNES DES NOEUDS DE LE ELEMENT PARENT
!       PINTER   : COORDONNES DES POINTS D INTERSECTION
!       PMILIE   : COORDONNES DES POINTS MILIEUX
!       NINTER   : NOMBRE DE POINTS D INTERSECTION
!       AINTER   : INFOS ARETE ASSOCIÉE AU POINTS D'INTERSECTION
!       AR       : CONNECTIVITE DU TETRA
!       PMMAX    : NOMBRE DE POINTS MILIEUX MAXIMAL DETECTABLE
!       NPTS     : NB DE PTS D'INTERSECTION COINCIDANT AVEC UN NOEUD SOMMET
!       LSNELP   : LSN AUX NOEUDS DE L'ELEMENT PARENT POUR LA FISSURE COURANTE
!       PINTT    : COORDONNEES REELLES DES POINTS D'INTERSECTION
!       PMITT    : COORDONNEES REELLES DES POINTS MILIEUX
!       JONC     : L'ELEMENT PARENT EST-IL TRAVERSE PAR PLUSIEURS FISSURES
!
!     SORTIE
!       NMILIE   : NOMBRE DE POINTS MILIEUX
!       PMILIE   : COORDONNES DES POINS MILIEUX
!     ----------------------------------------------------------------
!
    real(kind=8) :: milfi(3), milara(3), milarb(3)
    real(kind=8) :: milfa(3)
    real(kind=8) :: pmiref(12*ndime), ksia(ndime), ksib(ndime)
    integer(kind=8) :: n(3)
    integer(kind=8) :: i, ipm, k
    integer(kind=8) :: noeua
    integer(kind=8) :: j, r, ip, a2, a3, ip1(4), ip2(4), nbpi
    integer(kind=8) :: pm1a(4), pm1b(4), pm2(4)
    integer(kind=8) :: nm, inm, ia, ib, im, mfisloc
    integer(kind=8) :: zxain
    aster_logical :: ispm3, ispm2, ajout
!
! --------------------------------------------------------------------
!
    call jemarq()
!
    zxain = xxmmvd('ZXAIN')
!     COMPTEUR DES POINTS INTERSECTION NON CONFONDUS AVEC ND SOMMET
    ip = 0
!     COMPTEUR DE POINTS MILIEUX
    ipm = 0
    inm = 0
    mfisloc = 0
!
    pmilie(1:51) = 0.d0
!
    do i = 1, 4
        ip1(i) = 0
        ip2(i) = 0
        pm1a(i) = 0
        pm1b(i) = 0
        pm2(i) = 0
    end do
    call xstudo(ndime, ninter, npts, nptm, ainter, &
                nbpi, ip1, ip2, pm1a, pm1b, &
                pm2)
!    RECHERCHE DU NOEUD A
    noeua = 0
    a2 = nint(ainter(zxain*(2-1)+1))
    a3 = nint(ainter(zxain*(3-1)+1))
    do i = 1, 2
        do j = 1, 2
            if (ar(a2, i) .eq. ar(a3, j)) noeua = ar(a2, i)
        end do
    end do
    ASSERT(noeua .gt. 0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    RECHERCHE DU PREMIER TYPE DE POINT MILIEU    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do r = 1, ninter
        a2 = nint(ainter(zxain*(r-1)+1))
        ASSERT(a2 .ne. 0)
        ip = ip+1
        nm = ar(a2, 3)
        ia = 0
        ib = 0
!    ORDONANCEMENT DES NOEUDS MILIEUX SUR L ARETE : RECHERCHE DU NOEUD A SUR L ARETE A2
        do i = 1, 2
            if (ar(a2, i) .eq. noeua) then
                ia = cnset(nnose*(it-1)+ar(a2, 3-i))
                ib = cnset(nnose*(it-1)+ar(a2, i))
                im = cnset(nnose*(it-1)+ar(a2, 3))
            end if
        end do
        ASSERT((ia*ib) .gt. 0)
        milara(:) = 0.d0
        milarb(:) = 0.d0
        call xmilar(ndim, ndime, elrefp, geom, pinref, &
                    ia, ib, im, r, ksia, &
                    ksib, milara, milarb, pintt, pmitt)
!         STOCKAGE PMILIE
        call xajpmi(ndim, pmilie, pmmax, ipm, inm, &
                    milara, lonref, ajout)
        if (ajout) then
            do j = 1, ndime
                pmiref(ndime*(ipm-1)+j) = ksia(j)
            end do
        end if
        call xajpmi(ndim, pmilie, pmmax, ipm, inm, &
                    milarb, lonref, ajout)
        if (ajout) then
            do j = 1, ndime
                pmiref(ndime*(ipm-1)+j) = ksib(j)
            end do
        end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    RECHERCHE DU DEUXIEME TYPE DE POINT MILIEU    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1, nbpi
!      POINT MILIEU ENTRE IP1(R) ET IP2(R)
!
!      on ne calcule pas le premier type de point milieu si la
!      fissure coincide avec une arete
        ispm2 = (nint(ainter(zxain*(ip1(k)-1)+1)) .ne. 0) .or. &
                (nint(ainter(zxain*(ip2(k)-1)+1)) .ne. 0)
!
        if (ispm2) then
!        DETECTER LA COTE PORTANT LES DEUX POINTS D'INTERSECTIONS
            call detefa(nnose, ip1(k), ip2(k), it, typma, &
                        ainter, cnset, n)
!
!        CALCUL DU POINT MILIEU DE 101-102
!
            call xmifis(ndim, ndime, elrefp, geom, lsnelp, &
                        n, ip1(k), ip2(k), pinref, ksia, &
                        milfi, pintt, exit, jonc)
!
!        on incremente le nombre de points milieux sur la fissure
            mfisloc = mfisloc+1
!        STOCKAGE PMILIE
            call xajpmi(ndim, pmilie, pmmax, ipm, inm, &
                        milfi, lonref, ajout)
            if (ajout) then
                do j = 1, ndime
                    pmiref(ndime*(ipm-1)+j) = ksia(j)
                end do
            end if
        end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    RECHERCHE DU TROISIEME TYPE DE POINT MILIEU   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1, nbpi
!      POINT MILIEU ENTRE IP1(R) ET IP2(R)
!
!      on calcule le troisieme type de point milieu seulement si
!      aucun des deux points d'intersection n'est confondu avec
!      un noeud sommet
        ispm3 = ((nint(ainter(zxain*(ip1(k)-1)+1)) .ne. 0) .and. &
                 (nint(ainter(zxain*(ip2(k)-1)+1)) .ne. 0))
!
!
        if (ispm3) then
!        DETECTER LA COTE PORTANT LES DEUX POINTS D'INTERSECTIONS
            call detefa(nnose, ip1(k), ip2(k), it, typma, &
                        ainter, cnset, n)
!
            call xmilfa(elrefp, ndim, ndime, geom, cnset, &
                        nnose, it, ainter, ip1(k), ip2(k), &
                        pm2(k), typma, pinref, pmiref, ksia, &
                        milfa, pintt, pmitt)
!
            call xajpmi(ndim, pmilie, pmmax, ipm, inm, &
                        milfa, lonref, ajout)
            if (ajout) then
                do j = 1, ndime
                    pmiref(ndime*(ipm-1)+j) = ksia(j)
                end do
            end if
        end if
    end do
!
    ASSERT(ipm .eq. 12 .and. mfisloc .eq. 3)
!
    nmilie = ipm
    mfis = mfis+mfisloc
!
    call jedema()
end subroutine
