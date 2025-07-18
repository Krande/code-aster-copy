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
subroutine xalg42(ndim, elrefp, it, nnose, cnset, &
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
#include "asterfort/xstudo.h"
#include "asterfort/xxmmvd.h"
    character(len=8) :: typma, elrefp
    integer(kind=8) :: ndim, ndime, it, nnose, cnset(*), exit(2)
    integer(kind=8) :: ninter, pmmax, npts, nptm, nmilie, mfis, ar(12, 3)
    real(kind=8) :: lonref, ainter(*), pmilie(*), lsnelp(*)
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
    real(kind=8) :: pmiref(12), ksia(ndime), ksib(ndime)
    integer(kind=8) :: n(3)
    integer(kind=8) :: i, ipm, k, j
    integer(kind=8) :: r, a2, ip1(4), ip2(4), nbpi
    integer(kind=8) :: pm1a(4), pm1b(4), pm2(4)
    integer(kind=8) :: inm, ia, ib, im, mfisloc
    integer(kind=8) :: zxain
    aster_logical :: ispm2, ajout
!
! --------------------------------------------------------------------
!
    call jemarq()
!
    zxain = xxmmvd('ZXAIN')
!     COMPTEUR DES POINTS INTERSECTION NON CONFONDUS AVEC ND SOMMET
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    RECHERCHE DU PREMIER TYPE DE POINT MILIEU    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do r = 1, ninter
        a2 = nint(ainter(zxain*(r-1)+1))
! POUR EMPECHER L ALGO DE CALCULER TROP DE POINTS MILIEUX
! COMME LE NOEUD MILIEU EST EN 4EME POSITION QUAND NINTER=4,NPTS=2
        if (r .eq. 4) goto 300
        if (a2 .eq. 0) goto 300
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  REMARQUE IMPORTANTE ::                                                                         !
!  ON NE PEUT PAS DEFINIR UNE REFERENCE POUR L ORDONNACEMENT DANS CETTE CONFIGURATION DE DECOUPE  !
!  CAR LES DEUX DOMAINES DE DECOUPE SONT SYMETRIQUES (DEUX SOUS-TETRAS SYMETRIQUES)               !
!  L ORDONANCEMENT DES NOEUDS MILIEUX SUR L ARETE EST ALORS IMPLICITE                             !
!  IL DEPEND EN FAIT DE L ORDONANCEMMENT DES NOEUDS DANS LE TABLEAU AR                            !
!  ON A :                                                                                         !
!   IP1=CONNEC(AR(I,1))                                                                           !
!   IP2=CONNEC(AR(I,2))                                                                           !
!   PM1 EST MILIEU DE AR(A2,1) ET IP3                                                             !
!   PM2 EST MILIEU DE AR(A2,2) ET IP3                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ia = cnset(nnose*(it-1)+ar(a2, 1))
        ib = cnset(nnose*(it-1)+ar(a2, 2))
        im = cnset(nnose*(it-1)+ar(a2, 3))
        milara(:) = 0.d0
        milarb(:) = 0.d0
!    ORDONANCEMENT DES NOEUDS MILIEUX SUR L ARETE : RECHERCHE DU NOEUD A SUR L ARETE A2
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
300     continue
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
!
    ASSERT(ipm .eq. 4 .and. mfisloc .eq. 2)
!
    nmilie = ipm
    mfis = mfis+mfisloc
!
    call jedema()
end subroutine
