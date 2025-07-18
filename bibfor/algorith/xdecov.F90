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
subroutine xdecov(ndim, elp, nnop, nnose, it, &
                  pintt, cnset, heavt, ncomp, lsn, &
                  fisco, igeom, nfiss, ifiss, pinter, &
                  ninter, npts, ainter, nse, cnse, &
                  heav, nfisc, nsemax)
! aslint: disable=W1306,W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/xpente.h"
#include "asterfort/xxmmvd.h"
#include "blas/ddot.h"
    real(kind=8) :: lsn(*), pintt(*), pinter(*), ainter(*)
    integer(kind=8) :: ndim, nnop, nnose, it, cnset(*), heavt(*), ncomp, igeom
    integer(kind=8) :: ninter, npts, nfiss, ifiss, nse, cnse(6, 10), fisco(*)
    integer(kind=8) :: nsemax, nfisc
    real(kind=8) :: heav(*)
    character(len=8) :: elp
! person_in_charge: samuel.geniaut at edf.fr
!                      DÃCOUPER LE TETRA EN NSE SOUS-TETRAS
!
!     ENTREE
!       NNOSE    : NOMBRE DE NOEUDS DU SOUS TETRA
!       IT       : INDICE DU TETRA EN COURS
!       CNSET    : CONNECTIVITÃ DES NOEUDS DU TETRA
!       HEAVT    : FONCTION HEAVYSIDE DES TETRAS
!       LSN      : VALEURS DE LA LEVEL SET NORMALE
!       IGEOM    : ADRESSE DES COORDONNÃES DES NOEUDS DE L'ELT PARENT
!       PINTER   : COORDONNÃES DES POINTS D'INTERSECTION
!       NINTER   : NB DE POINTS D'INTERSECTION
!       NPTS     : NB DE PTS D'INTERSECTION COINCIDANT AVEC UN NOEUD
!       AINTER   : INFOS ARETE CORRESPONDATE AU PT INTERSECTION
!       NSEMAX   : NOMBRE DE SOUS-ELT MAX (TETRAS)
!
!     SORTIE
!       NSE      : NOMBRE DE SOUS-ELTS (TETRAS)
!       CNSE     : CONNECTIVITE DES SOUS-ÃLÃMENTS (TETRAS)
!       HEAV     : FONCTION HEAVYSIDE CONSTANTE SUR CHAQUE SOUS-ELT
!     ------------------------------------------------------------------
!
    real(kind=8) :: xyz(4, 3), ab(3), ac(3), ad(3), vn(3), ps, geom(3)
    real(kind=8) :: somlsn(nfisc+1), ff(nnop), rbid2(ndim)
    integer(kind=8) :: in, inh, i, j, ar(12, 3), nbar, ise, npent(18)
    integer(kind=8) :: a1, a2, a3, a4, a, b, c, d, ndime
    character(len=8) :: typma, elrese(3)
    integer(kind=8) :: zxain, mxstac
    aster_logical :: axi
    blas_int :: b_incx, b_incy, b_n
    parameter(mxstac=1000)
!
    data elrese/'SEG2', 'TRIA3', 'TETRA4'/
! ----------------------------------------------------------------------
!
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(nnop .le. mxstac)
    ASSERT(nfisc .le. mxstac)
    ASSERT(ndim .le. mxstac)
!
!
    call elrefe_info(fami='RIGI', ndim=ndime)
    zxain = xxmmvd('ZXAIN')
!
    axi = lteatt('AXIS', 'OUI')
!     ATTENTION, NE PAS CONFONDRE NDIM ET NDIME  !!
!     NDIM EST LA DIMENSION DU MAILLAGE
!     NDIME EST DIMENSION DE L'ELEMENT FINI
!     PAR EXEMPLE, POUR LES ELEMENT DE BORDS D'UN MAILLAGE 3D :
!     NDIME = 2 ALORS QUE NDIM = 3
!
    do in = 1, 6
        do j = 1, 6
            cnse(in, j) = 0
        end do
    end do
!
    typma = elrese(ndime)
    call conare(typma, ar, nbar)
!
!-----------------------------------------------------------------------
!     REMPLISSAGE DE LA CONNECTIVITÃ DES SOUS-ELEMENTS TÃTRAS
!                  ALGO BOOK III (26/04/04)
!-----------------------------------------------------------------------
!
    if (ndime .eq. 2) then
!
!
        if (ninter .lt. 2) then
!         INTER DOUTEUSE
            ASSERT(npts .eq. ninter)
!         1 SEUL ELEMENT
            nse = 1
            do in = 1, nnose
                cnse(1, in) = cnset(nnose*(it-1)+in)
            end do
        else if (ninter .eq. 2) then
            a1 = nint(ainter(zxain*(1-1)+1))
            a2 = nint(ainter(zxain*(2-1)+1))
            if (npts .eq. 2) then
!           1 SEUL ELEMENT
                nse = 1
                do in = 1, nnose
                    cnse(1, in) = cnset(nnose*(it-1)+in)
                end do
            else if (npts .eq. 1) then
!           2 ELEMENTS
                nse = 2
                ASSERT(a1 .eq. 0 .and. a2 .ne. 0)
                cnse(1, 1) = nint(ainter(zxain*(npts-1)+2))
                cnse(1, 2) = 102
                cnse(1, 3) = cnset(nnose*(it-1)+ar(a2, 1))
                cnse(2, 1) = nint(ainter(zxain*(npts-1)+2))
                cnse(2, 2) = 102
                cnse(2, 3) = cnset(nnose*(it-1)+ar(a2, 2))
            else
!           3 ELEMENTS
                nse = 3
                ASSERT(a1 .ne. 0)
!           101 ET 102 LES 2 POINTS D'INTERSECTION
!           ON SE PLACE DANS LA CONF DE REF (VOIR ALGO)
                do i = 1, 2
                    do j = 1, 2
                        if (ar(a1, i) .eq. ar(a2, j)) then
                            a = ar(a1, i)
                            b = ar(a1, 3-i)
                            c = ar(a2, 3-j)
                        end if
                    end do
                end do
                cnse(1, 1) = 101
                cnse(1, 2) = 102
                cnse(1, 3) = cnset(nnose*(it-1)+a)
                cnse(2, 1) = 101
                cnse(2, 2) = 102
                cnse(2, 3) = cnset(nnose*(it-1)+c)
                cnse(3, 1) = 101
                cnse(3, 2) = cnset(nnose*(it-1)+b)
                cnse(3, 3) = cnset(nnose*(it-1)+c)
            end if
        else if (ninter .eq. 3) then
!         L'INTERFACE COINCIDE AVEC LE TRIA
            ASSERT(npts .eq. ninter)
!         1 SEUL ELEMENT
            nse = 1
            do in = 1, nnose
                cnse(1, in) = cnset(nnose*(it-1)+in)
            end do
        else
!         TROP DE POINTS D'INTERSECTION
            ASSERT(ninter .le. 3)
        end if
!
    else if (ndime .eq. 1) then
!
        if (ninter .lt. 1) then
!         INTER DOUTEUSE
            ASSERT(npts .eq. ninter)
!         1 SEUL ELEMENT
            nse = 1
            do in = 1, 2
                cnse(1, in) = cnset(nnose*(it-1)+in)
            end do
        else if (ninter .eq. 1) then
            a1 = nint(ainter(zxain*(1-1)+1))
            if (npts .eq. 1) then
!           1 SEUL ELEMENT
                nse = 1
                do in = 1, 2
                    cnse(1, in) = cnset(nnose*(it-1)+in)
                end do
            else if (npts .eq. 0) then
!           2 ELEMENTS
                nse = 2
                ASSERT(a1 .ne. 0)
                a = ar(a1, 1)
                b = ar(a1, 2)
!
!           101 ET 102 LES 2 POINTS D'INTERSECTION
!           ON SE PLACE DANS LA CONF DE REF (VOIR ALGO)
                cnse(1, 1) = 101
                cnse(1, 2) = cnset(nnose*(it-1)+a)
                cnse(2, 1) = 101
                cnse(2, 2) = cnset(nnose*(it-1)+b)
            end if
        else if (ninter .eq. 2) then
!         L'INTERFACE COINCIDE AVEC LE SEG
            ASSERT(npts .eq. ninter)
!         1 SEUL ELEMENT
            nse = 1
            do in = 1, 2
                cnse(1, in) = cnset(nnose*(it-1)+in)
            end do
        else
!         TROP DE POINTS D'INTERSECTION
            ASSERT(ninter .le. 2)
        end if
!
!
!
    else if (ndime .eq. 3) then
!
        do i = 1, 18
            npent(i) = 0
        end do
        if (ninter .lt. 3) then
!
!       1Â°) AVEC MOINS DE TROIS POINTS D'INTERSECTION
!       ---------------------------------------------
!
!         INTER DOUTEUSE
            ASSERT(npts .eq. ninter)
!         ON A UN SEUL ELEMENT
            nse = 1
            do in = 1, nnose
                cnse(1, in) = cnset(nnose*(it-1)+in)
            end do
!
        else if (ninter .eq. 3) then
!
!         2Â°) AVEC TROIS POINTS D'INTERSECTION
!         ------------------------------------
            a1 = nint(ainter(zxain*(1-1)+1))
            a2 = nint(ainter(zxain*(2-1)+1))
            a3 = nint(ainter(zxain*(3-1)+1))
!
            if (npts .eq. 3) then
!           ON A UN SEUL ELEMENT
                nse = 1
                do in = 1, nnose
                    cnse(1, in) = cnset(nnose*(it-1)+in)
                end do
!
            else if (npts .eq. 2) then
!           ON A DEUX SOUS-ELEMENTS
                nse = 2
                ASSERT(a1 .eq. 0 .and. a2 .eq. 0)
                ASSERT(nint(ainter(2)) .gt. 0 .and. nint(ainter(zxain+2)) .gt. 0)
!
!           CONNECTIVITE DES NSE PAR RAPPORT AU NUM DE NOEUDS DU PARENT
!           AVEC 101, 102 ET 103 LES 3 PTS D'INTERSECTION
!           ON REMPLACE 101 ET 102 PAR LES NUMEROS DES NOEUDS COUPÉS
                cnse(1, 1) = nint(ainter(2))
                cnse(1, 2) = nint(ainter(zxain+2))
                cnse(1, 3) = 103
                cnse(1, 4) = cnset(nnose*(it-1)+ar(a3, 1))
                cnse(2, 1) = nint(ainter(2))
                cnse(2, 2) = nint(ainter(zxain+2))
                cnse(2, 3) = 103
                cnse(2, 4) = cnset(nnose*(it-1)+ar(a3, 2))
!
            else if (npts .eq. 1) then
!           ON A TROIS SOUS-ELEMENTS
                nse = 3
                ASSERT(a1 .eq. 0 .and. a2 .ne. 0)
                ASSERT(nint(ainter(2)) .gt. 0)
!           ON SE PLACE DANS LA CONF DE REF (VOIR ALGO)
                do i = 1, 2
                    do j = 1, 2
                        if (ar(a2, i) .eq. ar(a3, j)) then
                            a = ar(a2, i)
                            b = ar(a2, 3-i)
                            c = ar(a3, 3-j)
                        end if
                    end do
                end do
!           ON REMPLACE 101 PAR LE NUMERO DU NOEUD COUPÉ
                cnse(1, 1) = nint(ainter(2))
                cnse(1, 2) = 102
                cnse(1, 3) = 103
                cnse(1, 4) = cnset(nnose*(it-1)+a)
                cnse(2, 1) = nint(ainter(2))
                cnse(2, 2) = 102
                cnse(2, 3) = 103
                cnse(2, 4) = cnset(nnose*(it-1)+c)
                cnse(3, 1) = nint(ainter(2))
                cnse(3, 2) = 102
                cnse(3, 3) = cnset(nnose*(it-1)+b)
                cnse(3, 4) = cnset(nnose*(it-1)+c)
!
            else if (npts .eq. 0) then
!           ON A QUATRE SOUS-ELEMENTS
                nse = 4
                a = 0
                b = 0
                c = 0
                d = 0
                do i = 1, 2
                    do j = 1, 2
                        if (ar(a1, i) .eq. ar(a2, j)) then
                            a = ar(a1, i)
                            b = ar(a1, 3-i)
                            c = ar(a2, 3-j)
                        end if
                    end do
                end do
                do i = 1, 2
                    if (ar(a3, i) .ne. a) then
                        d = ar(a3, i)
                    end if
                end do
                ASSERT((a*b*c*d) .gt. 0)
!
                cnse(1, 1) = 101
                cnse(1, 2) = 102
                cnse(1, 3) = 103
                cnse(1, 4) = cnset(nnose*(it-1)+a)
!
                npent(1) = 103
                npent(2) = 101
                npent(3) = 102
                npent(4) = cnset(nnose*(it-1)+d)
                npent(5) = cnset(nnose*(it-1)+b)
                npent(6) = cnset(nnose*(it-1)+c)
                call xpente(2, cnse, npent)
            end if
!
        else if (ninter .eq. 4) then
!
!         2Â°) AVEC QUATRE POINTS D'INTERSECTION
!          -------------------------------------
            a1 = nint(ainter(zxain*(1-1)+1))
            a2 = nint(ainter(zxain*(2-1)+1))
            a3 = nint(ainter(zxain*(3-1)+1))
            a4 = nint(ainter(zxain*(4-1)+1))
!
!         ON A SIX SOUS-ELEMENTS (DANS TOUS LES CAS ?)
            nse = 6
            ASSERT((a1*a2*a3*a4) .gt. 0)
            a = 0
            b = 0
            c = 0
            d = 0
            do i = 1, 2
                do j = 1, 2
                    if (ar(a1, i) .eq. ar(a2, j)) then
                        a = ar(a1, i)
                        b = ar(a1, 3-i)
                        c = ar(a2, 3-j)
                    end if
                    if (ar(a3, i) .eq. ar(a4, j)) then
                        d = ar(a3, i)
                    end if
                end do
            end do
            ASSERT((a*b*c*d) .gt. 0)
            npent(1) = 104
            npent(2) = 102
            npent(3) = cnset(nnose*(it-1)+c)
            npent(4) = 103
            npent(5) = 101
            npent(6) = cnset(nnose*(it-1)+b)
            call xpente(1, cnse, npent)
            npent(1) = cnset(nnose*(it-1)+a)
            npent(2) = 101
            npent(3) = 102
            npent(4) = cnset(nnose*(it-1)+d)
            npent(5) = 103
            npent(6) = 104
            call xpente(4, cnse, npent)
        end if
    end if
!
!-----------------------------------------------------------------------
!     VÃRIFICATION DU SENS DES SOUS-ÃLÃMENTS TETRA
!                  ALGO BOOK III (28/04/04)
!-----------------------------------------------------------------------
!
    if (ndime .eq. 3) then
!
        do ise = 1, nse
            do in = 1, 4
                inh = cnse(ise, in)
                if (inh .lt. 100) then
                    do j = 1, 3
                        xyz(in, j) = zr(igeom-1+ndim*(inh-1)+j)
                    end do
                else if (inh .gt. 100 .and. inh .lt. 1000) then
                    do j = 1, 3
                        xyz(in, j) = pinter(ndim*(inh-100-1)+j)
                    end do
                else
                    do j = 1, 3
                        xyz(in, j) = pintt(ndim*(inh-1001)+j)
                    end do
                end if
            end do
!
            do j = 1, 3
                ab(j) = xyz(2, j)-xyz(1, j)
                ac(j) = xyz(3, j)-xyz(1, j)
                ad(j) = xyz(4, j)-xyz(1, j)
            end do
!
            call provec(ab, ac, vn)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            ps = ddot(b_n, vn, b_incx, ad, b_incy)
!
            if (ps .lt. 0) then
!          MAUVAIS SENS DU TETRA, ON INVERSE LES NOEUDS 3 ET 4
                inh = cnse(ise, 3)
                cnse(ise, 3) = cnse(ise, 4)
                cnse(ise, 4) = inh
            end if
!
        end do
!
!
    end if
!
!-----------------------------------------------------------------------
!             MATRICE DES COORDONNÃES ET FONCTION HEAVYSIDE
!             ALGO BOOK III (28/04/04)
!-----------------------------------------------------------------------
    ASSERT(nse .le. nsemax)
    do ise = 1, nse
        do i = 1, ifiss-1
! ----- ON RECOPIE LES VALEURS PRÉCÉDENTES
            heav(ifiss*(ise-1)+i) = heavt(ncomp*(i-1)+it)
        end do
! ----- ON TRAITE LA FISSURE COURANTE
        somlsn(:) = 0.d0
        do in = 1, nnose
            inh = cnse(ise, in)
            if (inh .lt. 100) then
                do i = 1, nfisc
                    somlsn(i) = somlsn(i)+lsn((inh-1)*nfiss+fisco(2*i-1))
                end do
                somlsn(nfisc+1) = somlsn(nfisc+1)+lsn((inh-1)*nfiss+ifiss)
            else
!           RECUP DE LA GÉOMETRIE
                if (inh .gt. 1000) then
                    do j = 1, ndim
                        geom(j) = pintt(ndim*(inh-1001)+j)
                    end do
                else if (inh .lt. 1000) then
                    do j = 1, ndim
                        geom(j) = pinter(ndim*(inh-101)+j)
                    end do
                end if
!           CALCUL DES FF
!
!
                call reeref(elp, nnop, zr(igeom), geom, ndim, &
                            rbid2, ff)
!
                do j = 1, nnop
                    do i = 1, nfisc
                        somlsn(i) = somlsn(i)+ff(j)*lsn((j-1)*nfiss+fisco(2*i-1))
                    end do
                    somlsn(nfisc+1) = somlsn(nfisc+1)+ff(j)*lsn((j-1)*nfiss+ifiss)
                end do
            end if
        end do
!
!       MISE À ZERO POUR LA FONCTION JONCTION AU NIVEAU DU BRANCHEMENT
!
        do i = 1, nfisc
            if (fisco(2*i)*somlsn(i) .gt. 0.d0) goto 300
        end do
!
        if (somlsn(nfisc+1) .lt. 0.d0) then
            heav(ifiss*ise) = -1.d0
        else if (somlsn(nfisc+1) .gt. 0.d0) then
            heav(ifiss*ise) = +1.d0
        else
!       REMARQUE IMPORTANTE :
!       SI ON EST SUR UN ELEMENT DE BORD COINCIDANT AVEC L'INTERCE
!       (NDIME = NDIM - 1 ET NPTS = NINTER = NDIM) ALORS ON NE PEUT PAS
!       DÃTERMINER DE QUEL COTÃ DE L'INTERFACE ON SE TROUVE, CAR
!       ON EST TOUJOURS SUR L'INTERFACE. LA VALEUR DE HEAV(ISE)
!       EST DONC FAUSSE DANS CE CAS : ON MET 99.
!       UNE CORRECTION EST FAITE DANS XORIPE LORS DE L'ORIENTATION DES
!       NORMALES, OU ON EN PROFITE POUR CORRIGER AUSSI HEAV(ISE)
            ASSERT(ndime .eq. ndim-1 .and. npts .eq. ndim .and. nse .eq. 1)
            heav(ifiss*ise) = 99.d0
        end if
!
300     continue
    end do
!
end subroutine
