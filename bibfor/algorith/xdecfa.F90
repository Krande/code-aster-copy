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
! person_in_charge: daniele.colombo at ifpen.fr
! aslint: disable=W1306,W1504
!
subroutine xdecfa(elp, nno, igeom, jlsn, jlst, &
                  npi, npis, pinter, pinref, ainter, &
                  cooree, cooref, rainter, noeud, npts, &
                  nintar, lst, lonref, ndim, zxain, &
                  jgrlsn, mipos)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
#include "asterfort/elrfdf.h"
#include "asterfort/iselli.h"
#include "asterfort/padist.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/xcenfi.h"
#include "asterfort/xnewto.h"
#include "asterfort/xnormv.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: npi, noeud(9), npis
    integer(kind=8) :: igeom, jlsn, jlst, zxain
    integer(kind=8) :: nintar, npts, ndim, nno, jgrlsn
    real(kind=8) :: pinter(*), ainter(*), cooree(6, ndim), cooref(6, ndim)
    real(kind=8) :: rainter(3, 4), lst(6), lonref, pinref(43*ndim)
    character(len=8) :: elp
    aster_logical :: mipos
!
! ======================================================================
!            BUT :  TROUVER LES PTS DE LA FISSURE DANS LE TRIA
!
!     ENTREE
!       ELP      : MACRO ELEMENT PARENT
!       NNO      : NOMBRE DE NOEUD DE L'ELEMENT PARENT
!       NPI      : NOMBRE DE NOEUDS DISTINCTS DE LA FISSURE
!       PINTER   : COORDONNES DES POINTS DE LA FISSURE
!       AINTER   : INFOS ARETE ASSOCIÉE AU POINTS D'INTERSECTION
!       COOREE   : COOREDONNEES REELLES DES SOMMETS DU TRIA
!       COOREF   : COOREDONNEES DE REFERENCE DES SOMMETS DU TRIA
!       RAINTER  : INFO ARETES ASSOCIES AUX POINTS D'INTERSECTION
!       LST      : LST AUX NOEUDS SOMMETS DU TRIA
!       LONREF   : LONGUEUR CARACTERISTIQUE DE L'ELEMENT
!       NDIM     : DIMENSION DU MACRO ELEMENT
!       I        : SOUS ELEMENT COURANT
!       FACE     : FACE EN COURS DANS LE SOUS SOUS ELEMENT
!       F        : CONNECTIVITE DU SOUS SOUS ELEMENT
!       NNOSE    : NOMBRE DE NOEUDS PAR SOUS ELEMENT
!
!     SORTIE
!       NPI      : NOMBRE DE NOEUDS DISTINCTS DE LA FISSURE
!       PINTER   : COORDONNES DES POINTS DE LA FISSURE
!       AINTER   : INFOS ARETE ASSOCIÉE AU POINTS D'INTERSECTION
!       NOEUD    : NUMERO DES NOEUDS DE LA FISSURE DANS L'ELEMENT
!       NPTS     : NOMBRE DE NOEUDS SOMMETS DU TRIA TELS QUE LST<=0
!       NINTAR   : NOMBRE d'ARETES DU TRIA STRICTEMENT COUPEES PAR LST
!     ----------------------------------------------------------------
!
    real(kind=8) :: p(ndim), newpt(ndim), newptref(ndim), cenref(ndim)
    real(kind=8) :: norme, geom(ndim*nno), ff(27), cenfi(ndim), tabls(20)
    real(kind=8) :: x(ndim), xref(ndim), miref(ndim), mifis(ndim), ptxx(3*ndim)
    real(kind=8) :: vectn(ndim), ksi(ndim), dff(3, 27)
    real(kind=8) :: epsmax, cridist, a, b, c, ab(ndim), bc(ndim), gradlsn(ndim)
    real(kind=8) :: normfa(ndim), det, tempo, temp1(ndim), temp2(ndim), temp3(4)
    integer(kind=8) :: k, ii, jj, j, ni, kk, ibid, num(8)
    integer(kind=8) :: n(3), kkk, nn(4), exit(2)
    integer(kind=8) :: itemax
    aster_logical :: deja, jonc
    blas_int :: b_incx, b_incy, b_n
    parameter(cridist=1.d-7)
!
! --------------------------------------------------------------------
!
!
!   INITIALISATION DU NOMBRE D'ARETES DU TRIA INTERSECTEES
!   STRICTEMENT PAR LST
    nintar = 0
!
!   INITIALISATION DU NOMBRE DE SOMMETS DU TRIA TELS QUE LST<=0
    npts = 0
!
!   BOUCLE SUR LES SOMMETS DU TRIA
    newpt(:) = 0.d0
    tabls(:) = 0.d0
    geom(:) = 0.d0
    ksi(:) = 0.d0
    do j = 1, 9
        noeud(j) = 0
    end do
!
!   NECESSITE D'INVERSER LA CONNECTIVITE DE LA FACE A DECOUPER SI ELLE N'EST PAS
!   ORIENTEE SUIVANT GRADLSN
    do j = 1, ndim
        ab(j) = cooree(2, j)-cooree(1, j)
        bc(j) = cooree(3, j)-cooree(1, j)
        gradlsn(j) = zr(jgrlsn-1+j)
        normfa(j) = 0.d0
    end do
    call provec(ab, bc, normfa)
    call xnormv(ndim, normfa, norme)
    call xnormv(ndim, gradlsn, norme)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    det = ddot(b_n, gradlsn, b_incx, normfa, b_incy)
    if (det .lt. 0.d0) then
        tempo = lst(2)
        lst(2) = lst(3)
        lst(3) = tempo
        do j = 1, ndim
            temp1(j) = cooree(2, j)
            temp2(j) = cooref(2, j)
            cooree(2, j) = cooree(3, j)
            cooref(2, j) = cooref(3, j)
            cooree(3, j) = temp1(j)
            cooref(3, j) = temp2(j)
        end do
        do j = 1, 4
            temp3(j) = rainter(2, j)
            rainter(2, j) = rainter(3, j)
            rainter(3, j) = temp3(j)
        end do
        if (.not. iselli(elp)) then
            tempo = lst(4)
            lst(4) = lst(6)
            lst(6) = tempo
            do j = 1, ndim
                temp1(j) = cooree(4, j)
                temp2(j) = cooref(4, j)
                cooree(4, j) = cooree(6, j)
                cooref(4, j) = cooref(6, j)
                cooree(6, j) = temp1(j)
                cooref(6, j) = temp2(j)
            end do
        end if
    end if
!
!   BOUCLE SUR LES SOMMETS DU TRIA
    do k = 1, 3
        if (lst(k) .le. 0.d0) then
            npts = npts+1
            do ii = 1, ndim
                newpt(ii) = cooree(k, ii)
            end do
!      VERIF SI DEJA
            deja = .false.
            do ii = 1, npi
                do j = 1, ndim
                    p(j) = pinter(ndim*(ii-1)+j)
                end do
                if (padist(ndim, p, newpt) .lt. (lonref*cridist)) then
                    deja = .true.
                    ni = ii
                end if
            end do
!      ON ARCHIVE LES NOEUDS SOMMETS DU TRIA TELS QUE LST<=0
            if (.not. deja) then
                npi = npi+1
                npis = npis+1
                do j = 1, ndim
                    pinter(ndim*(npi-1)+j) = newpt(j)
                    pinref(ndim*(npi-1)+j) = cooref(k, j)
                end do
                do j = 1, zxain-1
                    ainter(zxain*(npi-1)+j) = rainter(k, j)
                end do
                noeud(npts) = npi
            else
                noeud(npts) = ni
            end if
        end if
    end do
    do ii = 1, nno
        tabls(ii) = zr(jlst-1+ii)
        do jj = 1, ndim
            geom((ii-1)*ndim+jj) = zr(igeom-1+ndim*(ii-1)+jj)
        end do
    end do
!
!      ON BOUCLE SUR LES ARETES DU TRIA POUR RECUPERER LES POINTS DU FOND DE
!      FISSURE SUR CES ARETES
    do k = 1, 3
        if (k .eq. 1) then
            kk = 2
            kkk = 3
        else if (k .eq. 2) then
            kk = 3
            kkk = 1
        else if (k .eq. 3) then
            kk = 1
            kkk = 2
        end if
!      SI L'ARETE EST COUPE PAR LE FOND DE FISSURE
        if (lst(k)*lst(kk) .lt. 0.d0) then
            nintar = nintar+1
!      RECHERCHE DU FOND DE FISSURE SUR L'ARETE
            do jj = 1, ndim
                ptxx(jj) = cooref(k, jj)
                ptxx(jj+ndim) = cooref(kk, jj)
            end do
            ASSERT(abs(lst(k)-lst(kk)) .gt. 1.d0/r8gaem())
            if (.not. iselli(elp)) then
!      RECUPERATION DU NOEUD MILIEU DE L'ARETE
                do jj = 1, ndim
                    ptxx(2*ndim+jj) = cooref(k+ndim, jj)
                end do
!      INITIALISATION DU NEWTON
                a = (lst(k)+lst(kk)-2.d0*lst(k+ndim))/2.d0
                b = (lst(kk)-lst(k))/2.d0
                c = lst(k+ndim)
                ASSERT(b**2 .ge. (4*a*c))
                if (abs(a) .lt. 1.d-8) then
                    ksi(1) = lst(k)/(lst(k)-lst(kk))
                else
                    ksi(1) = (-b-sqrt(b**2-4.d0*a*c))/(2.d0*a)
                    if (abs(ksi(1)) .gt. 1.d0) ksi(1) = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)
                    ASSERT(abs(ksi(1)) .le. 1.d0)
                    ksi(1) = (ksi(1)+1.d0)/2.d0
                end if
            else
                ksi(1) = lst(k)/(lst(k)-lst(kk))
                do jj = 1, ndim
                    ptxx(2*ndim+jj) = (ptxx(jj)+ptxx(ndim+jj))/2.d0
                end do
            end if
            do ii = 1, 3
                n(ii) = ii
            end do
            epsmax = 1.d-8
            itemax = 100
!      ALGORITHME DE NEWTON POUR TROUVER LE FOND DE FISSURE
            call xnewto(elp, 'XINTER', n, ndim, ptxx, &
                        ndim, geom, tabls, ibid, ibid, &
                        itemax, epsmax, ksi)
            xref(:) = 0.d0
            do ii = 1, ndim
                xref(ii) = 2.d0*( &
                           1.d0-ksi(1))*(5.d-1-ksi(1))*ptxx(ii)+4.d0*ksi(1)*(1.d0-ksi(1))*ptxx(i&
                           &i+2*ndim)+2.d0*ksi(1)*(ksi(1)-5.d-1)*ptxx(ii+ndim &
                           )
            end do
            ff(:) = 0.d0
            call elrfvf(elp, xref, ff, nno)
!      CALCUL DES COORDONNEES REELES DU FOND DE FISSURE
            x(:) = 0.d0
            do ii = 1, ndim
                do j = 1, nno
                    x(ii) = x(ii)+zr(igeom-1+ndim*(j-1)+ii)*ff(j)
                end do
            end do
!      VERIF SI DEJA
            deja = .false.
            do ii = 1, npi
                do j = 1, ndim
                    p(j) = pinter(ndim*(ii-1)+j)
                end do
                if (padist(ndim, p, x) .lt. (lonref*cridist)) then
                    deja = .true.
                    ni = ii
                end if
            end do
            if (.not. deja) then
                npi = npi+1
                npis = npis+1
                do j = 1, ndim
                    pinter(ndim*(npi-1)+j) = x(j)
                    pinref(ndim*(npi-1)+j) = xref(j)
                end do
                do j = 1, zxain-1
                    ainter(zxain*(npi-1)+j) = 0.d0
                end do
                noeud(2+nintar) = npi
            else
                noeud(2+nintar) = ni
            end if
!      DANS LE CAS QUADRATIQUE ON CHERCHE LE POINT MILEU CORRESP0NDANT (SUR
!      l'ARETE INTERSECTEE)
            if (.not. iselli(elp)) then
                if (lst(k) .lt. 0.d0) then
                    ksi(1) = ksi(1)/2.d0
                else
                    ksi(1) = (1.d0+ksi(1))/2.d0
                end if
                do ii = 1, ndim
                    miref(ii) = 2.d0*( &
                                1.d0-ksi(1))*(5.d-1-ksi(1))*ptxx(ii)+4.d0*ksi(1)*(1.d0-ksi(1))*p&
                                &txx(ii+2*ndim)+2.d0*ksi(1)*(ksi(1)-5.d-1)*ptxx(ii+ndim &
                                )
                end do
                call elrfvf(elp, miref, ff, nno)
                mifis(:) = 0.d0
                do ii = 1, ndim
                    do j = 1, nno
                        mifis(ii) = mifis(ii)+zr(igeom-1+ndim*(j-1)+ii)*ff(j)
                    end do
                end do
!      VERIF SI DEJA
                deja = .false.
                do ii = 1, npi
                    do j = 1, ndim
                        p(j) = pinter(ndim*(ii-1)+j)
                    end do
                    if (padist(ndim, p, mifis) .lt. (lonref*cridist)) then
                        deja = .true.
                        ni = ii
                    end if
                end do
                if (.not. deja) then
                    npi = npi+1
                    do j = 1, ndim
                        pinter(ndim*(npi-1)+j) = mifis(j)
                        pinref(ndim*(npi-1)+j) = miref(j)
                    end do
                    do j = 1, zxain-1
                        ainter(zxain*(npi-1)+j) = 0.d0
                    end do
                    noeud(4+nintar) = npi
                else
                    noeud(4+nintar) = ni
                end if
            end if
!      DANS LE CAS QUADRATIQUE ON CHERCHE LES MILIEUX D'ARETES LST<=0
        else if (lst(k) .le. 0.d0 .and. lst(kk) .le. 0.d0) then
            if (.not. iselli(elp)) then
                ASSERT(npts .eq. 2)
                do j = 1, ndim
                    newpt(j) = cooree(k+ndim, j)
                end do
                call reeref(elp, nno, zr(igeom), newpt, ndim, &
                            newptref, ff)
!      VERIF SI DEJA
                deja = .false.
                do ii = 1, npi
                    do j = 1, ndim
                        p(j) = pinter(ndim*(ii-1)+j)
                    end do
                    if (padist(ndim, p, newpt) .lt. (lonref*cridist)) then
                        deja = .true.
                        ni = ii
                    end if
                end do
                if (.not. deja) then
                    npi = npi+1
                    do j = 1, ndim
                        pinter(ndim*(npi-1)+j) = newpt(j)
                        pinref(ndim*(npi-1)+j) = newptref(j)
                    end do
                    do j = 1, zxain-1
                        ainter(zxain*(npi-1)+j) = 0.d0
                    end do
                    noeud(7) = npi
                else
                    noeud(7) = ni
                end if
            end if
        end if
    end do
!   ON DETERMINE MAINTENANT DANS LE CAS QUADRATIQUE LE NOEUD MILIEU
!   ENTRE LES DEUX POINTS DU FOND DE FISSURE
!
    if (.not. iselli(elp)) then
        if (mipos .and. npts .eq. 2) then
            do j = 1, ndim
                xref(j) = (pinref((noeud(3)-1)*ndim+j)+pinref((noeud(4)-1)*ndim+j))/2.d0
            end do
            x(:) = 0.d0
            call elrfvf(elp, xref, ff, nno)
            do j = 1, ndim
                do ii = 1, nno
                    x(j) = x(j)+zr(igeom-1+ndim*(ii-1)+j)*ff(ii)
                end do
            end do
        else
            if (nintar .eq. 2) then
                do j = 1, ndim
                    xref(j) = (pinref((noeud(4)-1)*ndim+j)+pinref((noeud(3)-1)*ndim+j))/2.d0
                end do
            else if (lst(1) .eq. 0.d0) then
                do j = 1, ndim
                    xref(j) = (pinref((noeud(3)-1)*ndim+j)+pinref((noeud(1)-1)*ndim+j))/2.d0
                end do
            else
                do j = 1, ndim
                    xref(j) = (pinref((noeud(3)-1)*ndim+j)+pinref((noeud(2)-1)*ndim+j))/2.d0
                end do
            end if
!   ON RECHERCHE SUR LA MEDIATRICE DU SEGMENT IP1IP2 PORTEE PAR GRADLST
            vectn(:) = 0.d0
            call elrfdf(elp, xref, dff, nno, ndim)
            do ii = 1, ndim
                do j = 1, nno
                    vectn(ii) = vectn(ii)+dff(ii, j)*tabls(j)
                end do
            end do
            call xnormv(ndim, vectn, norme)
            do j = 1, ndim
                ptxx(j) = vectn(j)
                ptxx(ndim+j) = xref(j)
            end do
            call xnewto(elp, 'XMIFIS', n, ndim, ptxx, &
                        ndim, geom, tabls, ibid, ibid, &
                        itemax, epsmax, ksi)
            do j = 1, ndim
                xref(j) = xref(j)+ksi(1)*ptxx(j)
            end do
            call elrfvf(elp, xref, ff, nno)
            x(:) = 0.d0
            do ii = 1, ndim
                do j = 1, nno
                    x(ii) = x(ii)+zr(igeom-1+ndim*(j-1)+ii)*ff(j)
                end do
            end do
        end if
!   ON ARCHIVE CE POINT
        npi = npi+1
        do j = 1, ndim
            pinter(ndim*(npi-1)+j) = x(j)
            pinref(ndim*(npi-1)+j) = xref(j)
        end do
        noeud(8) = npi
        do j = 1, zxain-1
            ainter(zxain*(npi-1)+j) = 0.d0
        end do
!   DANS LE CAS NPTS=2 ET NINTAR=2, IL NOUS RESTE ENCORE UN POINT MILIEU A
!   DETERMINER (DANS LE QUAD TEL QUE LST<=0)
        if (npts .eq. 2 .and. nintar .eq. 2.) then
!   ORDONANCEMENT DES POINTS D'INTERSECTION SUR LA FACE QUADRANGLE
            if (lst(3) .gt. 0.d0) then
                num(1) = noeud(1)
                num(2) = noeud(2)
                num(3) = noeud(4)
                num(4) = noeud(3)
                num(5) = noeud(7)
                num(6) = noeud(6)
                num(7) = noeud(8)
                num(8) = noeud(5)
            else if (lst(2) .gt. 0.d0) then
                num(1) = noeud(2)
                num(2) = noeud(1)
                num(3) = noeud(4)
                num(4) = noeud(3)
                num(5) = noeud(7)
                num(6) = noeud(6)
                num(7) = noeud(8)
                num(8) = noeud(5)
            else
                num(1) = noeud(1)
                num(2) = noeud(2)
                num(3) = noeud(3)
                num(4) = noeud(4)
                num(5) = noeud(7)
                num(6) = noeud(5)
                num(7) = noeud(8)
                num(8) = noeud(6)
            end if
            nn(1:4) = 0
            jonc = .true.
            exit(1:2) = 0
            call xcenfi(elp, ndim, ndim, nno, geom, &
                        zr(jlsn), pinref, pinref, cenref, cenfi, &
                        nn, exit, jonc, num)
!   ON ARCHIVE CE POINT
            npi = npi+1
            do j = 1, ndim
                pinter(ndim*(npi-1)+j) = cenfi(j)
                pinref(ndim*(npi-1)+j) = cenref(j)
            end do
            noeud(9) = npi
            do j = 1, zxain-1
                ainter(zxain*(npi-1)+j) = 0.d0
            end do
        end if
    end if
!
!
end subroutine
