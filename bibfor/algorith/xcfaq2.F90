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
subroutine xcfaq2(jlsn, jlst, jgrlsn, igeom, noma, &
                  nmaabs, pinter, ainter, nface, nptf, &
                  cface, nbtot, nfiss, ifiss)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/abscvf.h"
#include "asterfort/abscvl.h"
#include "asterfort/assert.h"
#include "asterfort/conare.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/elrfvf.h"
#include "asterfort/loncar.h"
#include "asterfort/padist.h"
#include "asterfort/tecael.h"
#include "asterfort/xajpin.h"
#include "asterfort/xcfacf.h"
#include "asterfort/xintar.h"
#include "asterfort/xinvac.h"
#include "asterfort/xmilfi.h"
#include "asterfort/xxmmvd.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: jgrlsn, igeom, nface, cface(30, 6), jlsn, jlst
    integer(kind=8) :: nfiss, ifiss, nptf, nbtot, nmaabs
    real(kind=8) :: pinter(*), ainter(*)
    character(len=8) :: noma
!                     TROUVER LES PTS D'INTERSECTION ENTRE LES ARETES,
!                     ET LE PLAN DE FISSURE, DÉCOUPAGE EN FACETTES,
!                     POINT MILIEU DE FISSURE (UNIQUEMENT 2D)
!     ENTREE
!       LSN      : VALEURS DE LA LEVEL SET NORMALE
!       LST      : VALEURS DE LA LEVEL SET TANGENTE
!       JGRLSN   : ADRESSE DU GRADIENT DE LA LEVEL SET NORMALE
!       IGEOM    : ADRESSE DES COORDONNEES DES NOEUDS DE L'ELT PARENT
!       NOMA     : NOM DU MAILLAGE
!       NMAABS   : INDICE DE LA MAILLE
!
!     SORTIE
!       PTINTER  : COORDONNEES DES POINTS D'INTERSECTION
!       NINTER  : NOMBRE DE POINTS D'INTERSECTION
!       AINTER  : INFOS ARETE ASSOCIEE AU POINTS D'INTERSECTION
!       NFACE   : NOMBRE DE FACETTES
!       NPTF    : NOMBRE DE POINTS PAR FACETTE
!       CFACE   : CONNECTIVITE DES NOEUDS DES FACETTES
!
!     ----------------------------------------------------------------
!
    real(kind=8) :: a(3), b(3), c(3), lsna, lsnb, longar, tampor(4)
    real(kind=8) :: alpha, nd(3), coor2d(9)
    real(kind=8) :: ab(2), lsta, lstb, lstc, abprim(2), prec, lonref, cridist
    real(kind=8) :: eps
    real(kind=8) :: ff(27), ksic(3), sc, tabar(9), minlsn, maxlsn
    real(kind=8) :: m(3), lsnm, lstm, ksi, milfi(3), smilfi, lsnabs
    integer(kind=8) :: j, ar(12, 3), nbar, na, nb, ins, n(3)
    integer(kind=8) :: ia, i, ipt, nno, k
    integer(kind=8) :: iadzi, iazk24, ndim, ptmax
    integer(kind=8) :: zxain
    integer(kind=8) :: inm, inc, nm, ninter
    parameter(cridist=1.d-7)
    aster_logical :: cut, ajout, arete
    character(len=8) :: typma, elp, elc
    blas_int :: b_incx, b_incy, b_n
!
    parameter(ptmax=4, elc='SE3')
! --------------------------------------------------------------------
!
!
    eps = -1.0d-10
!
!     PREC PERMET D'EVITER LES ERREURS DE PRECISION CONDUISANT
!     A IA=IN=0 POUR LES MAILLES DU FRONT
    prec = 1000*r8prem()
    zxain = xxmmvd('ZXAIN')
    call elref1(elp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno)
    ASSERT(ndim .eq. 2)
!
!     1) RECHERCHE DES POINTS D'INTERSECTION
!     --------------------------------------
!
!     VECTEUR REEL À ZXAIN COMPOSANTES, POUR CHAQUE PT D'INTER :
!     - NUMÉRO ARETE CORRESPONDANTE         (0 SI C'EST UN NOEUD SOMMET)
!     - NUMÉRO NOEUD SI NOEUD SOMMET        (0 SINON)
!     - LONGUEUR DE L'ARETE
!     - POSITION DU PT SUR L'ARETE          (0 SI C'EST UN NOEUD SOMMET)
!     - ARETE VITALE                        (0 SI NON)
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
    call loncar(ndim, typma, zr(igeom), lonref)
!
!     L'ELEMENT EST-IL TRAVERSE STRICTEMENT PAR LSN=0?
    nbtot = 0
    cut = .false.
    i = 1
1   continue
!     (1) RECHERCHE D'UN NOEUD PIVOT (LSN NON NULLE)
    if (zr(jlsn-1+(i-1)*nfiss+ifiss) .ne. 0.d0 .and. i .lt. nno) then
        do k = i+1, nno
!     (2) PRODUIT DE CE PIVOT PAR LES AUTRES LSN
            if (zr(jlsn-1+(i-1)*nfiss+ifiss)*zr(jlsn-1+(k-1)*nfiss+ifiss) .lt. 0.d0) cut = &
                .true.
        end do
    else if (i .lt. nno) then
        i = i+1
        goto 1
    end if
!     RECHERCHE DE MINLSN
    minlsn = 0.d0
    maxlsn = 0.d0
    do i = 1, nno
        minlsn = min(zr(jlsn-1+(i-1)*nfiss+ifiss), minlsn)
        maxlsn = max(zr(jlsn-1+(i-1)*nfiss+ifiss), maxlsn)
    end do
!
!     ON NE PREND QUE CERTAINS ELEMENTS POUR NE PAS AVOIR DE "DOUBLONS"
    arete = .false.
    lsnabs = 0.d0
    if (.not. cut) then
        call conare(typma, ar, nbar)
        do i = 1, nbar
            lsnabs = abs( &
                     zr(jlsn-1+(ar(i, 1)-1)*nfiss+ifiss))+abs(zr(jlsn-1+(ar(i, 2)-1)*nfiss+ifiss) &
                                                              )
            if (lsnabs .le. cridist*lonref) arete = .true.
        end do
        if (.not. arete) goto 999
        if (arete .and. minlsn .ge. 0.d0) goto 999
    end if
    ipt = 0
!     COMPTEUR DE POINT INTERSECTION = NOEUD SOMMET
    ins = 0
!     COMPTEUR DE POINT INTERSECTION = POINT ARETE
    inc = 0
!     COMPTEUR DE POINT INTERSECTION = NOEUD MILIEU
    inm = 0
    call conare(typma, ar, nbar)
!
!     BOUCLE SUR LES ARETES POUR DETERMINER LES POINTS D'INTERSECTION
    do ia = 1, nbar
!
!       NUM NO DE L'ELEMENT
        na = ar(ia, 1)
        nb = ar(ia, 2)
        nm = ar(ia, 3)
        lsna = zr(jlsn-1+na)
        lsnb = zr(jlsn-1+nb)
        lsnm = zr(jlsn-1+nm)
        lsta = zr(jlst-1+na)
        lstb = zr(jlst-1+nb)
        lstm = zr(jlst-1+nm)
        do i = 1, ndim
            a(i) = zr(igeom-1+ndim*(na-1)+i)
            b(i) = zr(igeom-1+ndim*(nb-1)+i)
            m(i) = zr(igeom-1+ndim*(nm-1)+i)
        end do
        if (ndim .lt. 3) then
            a(3) = 0.d0
            b(3) = 0.d0
            c(3) = 0.d0
            m(3) = 0.d0
        end if
!
        ksi = 1.d0
        longar = padist(ndim, a, b)
!
        if ((lsna*lsnb) .le. 0.d0) then
            if ((lsna .eq. 0.d0) .and. (lsta .le. prec)) then
!           ON AJOUTE A LA LISTE LE POINT A
                if (lsta .ge. 0.d0) then
                    call xajpin(ndim, pinter, ptmax, ipt, ins, &
                                a, longar, ainter, 0, 0, &
                                0.d0, ajout)
                else
                    call xajpin(ndim, pinter, ptmax, ipt, ins, &
                                a, longar, ainter, 0, na, &
                                0.d0, ajout)
                end if
            end if
            if (lsnb .eq. 0.d0 .and. lstb .le. prec) then
!           ON AJOUTE A LA LISTE LE POINT B
                if (lstb .ge. 0.d0) then
                    call xajpin(ndim, pinter, ptmax, ipt, ins, &
                                b, longar, ainter, 0, 0, &
                                0.d0, ajout)
                else
                    call xajpin(ndim, pinter, ptmax, ipt, ins, &
                                b, longar, ainter, 0, nb, &
                                0.d0, ajout)
                end if
            end if
!
            if (lsnm .eq. 0.d0 .and. lstm .le. prec) then
!           ON AJOUTE A LA LISTE LE POINT M
                alpha = padist(ndim, a, m)
                if (lstm .ge. 0.d0) then
                    call xajpin(ndim, pinter, ptmax, ipt, inm, &
                                m, longar, ainter, 0, 0, &
                                0.d0, ajout)
                else
                    if (cut) then
                        call xajpin(ndim, pinter, ptmax, ipt, inc, &
                                    m, longar, ainter, ia, 0, &
                                    alpha, ajout)
                    else if (.not. cut) then
                        call xajpin(ndim, pinter, ptmax, ipt, inm, &
                                    m, longar, ainter, 0, nm, &
                                    alpha, ajout)
                    end if
                end if
            end if
!
            if (lsna .ne. 0.d0 .and. lsnb .ne. 0.d0 .and. lsnm .ne. 0) then
!           INTERPOLATION DES COORDONNÉES DE C
                call xintar(lsna, lsnb, lsnm, a, b, &
                            m, ndim, c)
!           POSITION DU PT D'INTERSECTION SUR L'ARETE
                alpha = padist(ndim, a, c)
                do i = 1, ndim
                    tabar(i) = b(i)
                    tabar(ndim+i) = a(i)
                    tabar(2*ndim+i) = m(i)
                end do
!          CALCUL DES FF DU SE3 (REEREF N'ACCEPTE PAS NDIM=2 & NNO=3)
                call abscvl(ndim, tabar, c, sc)
                call xinvac(elp, ndim, tabar, sc, ksic)
                ASSERT(ksic(1) .ge. -1 .and. ksic(1) .le. 1)
                call elrfvf(elc, ksic(1), ff)
                lstc = ff(1)*lstb+ff(2)*lsta+ff(3)*lstm
                if (lstc .le. prec) then
                    if (lstc .ge. 0.d0) then
                        call xajpin(ndim, pinter, ptmax, ipt, inc, &
                                    c, longar, ainter, 0, 0, &
                                    0.d0, ajout)
                    else
                        call xajpin(ndim, pinter, ptmax, ipt, inc, &
                                    c, longar, ainter, ia, 0, &
                                    alpha, ajout)
                    end if
                end if
            end if
!
        end if
!
    end do
!
!     RECHERCHE SPECIFIQUE POUR LES ELEMENTS EN FOND DE FISSURE
    call xcfacf(pinter, ptmax, ipt, ainter, zr(jlsn), &
                zr(jlst), igeom, nno, ndim, typma, &
                noma, nmaabs)
!
    ninter = ins+inc
    nbtot = ninter+inm
!
    if (cut .and. ninter .eq. 2 .and. inm .ne. 1) then
        do j = 1, 3
            n(j) = 0
        end do
! LE NEWTON NE PEUT PAS CONVERGER QUAND NDIM=3
! IL NOUS MANQUE L INFORMATION SUR LA FACE DE L ELT PARENT
! TRANSPORTEE PAS N(3)
        ASSERT(ndim .lt. 3)
!       RECHERCHE POINT MILIEU FISSURE
        call xmilfi(elp, n, ndim, nno, pinter, &
                    ndim, igeom, jlsn, 1, 2, &
                    milfi)
        do j = 1, ndim
            coor2d(j) = pinter(j)
            coor2d(ndim+j) = pinter(ndim+j)
            coor2d(2*ndim+j) = milfi(j)
        end do
        ksi = 0.d0
        call abscvf(ndim, coor2d, ksi, smilfi)
!       ON AJOUTE A LA LISTE LE POINT MILFI
        call xajpin(ndim, pinter, ptmax, ipt, nbtot, &
                    milfi, smilfi*2, ainter, 5, 0, &
                    smilfi, ajout)
    end if
!
!     2) DECOUPAGE EN FACETTES TRIANGULAIRES DE LA SURFACE DEFINIE
!     ------------------------------------------------------------
!
!                  (BOOK IV 09/09/04)
!
!     CAS 2D
    if (ndim .eq. 2) then
!
        do i = 1, 30
            do j = 1, 6
                cface(i, j) = 0
            end do
        end do
        if (ninter .eq. 2) then
!         NORMALE A LA FISSURE (MOYENNE DE LA NORMALE AUX NOEUDS)
            nd(:) = 0.d0
            do i = 1, nno
                do j = 1, 2
                    nd(j) = nd(j)+zr(jgrlsn-1+2*(i-1)+j)/nno
                end do
            end do
!
            do j = 1, 2
                a(j) = pinter(j)
                b(j) = pinter(2+j)
                ab(j) = b(j)-a(j)
            end do
!
            abprim(1) = -ab(2)
            abprim(2) = ab(1)
!
            b_n = to_blas_int(2)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            if (ddot(b_n, abprim, b_incx, nd, b_incy) .lt. 0.d0) then
                do k = 1, 2
                    tampor(k) = pinter(k)
                    pinter(k) = pinter(2+k)
                    pinter(2+k) = tampor(k)
                end do
                do k = 1, 4
                    tampor(k) = ainter(k)
                    ainter(k) = ainter(zxain+k)
                    ainter(zxain+k) = tampor(k)
                end do
            end if
            nface = 1
            nptf = 3
            cface(1, 1) = 1
            cface(1, 2) = 2
            cface(1, 3) = 3
        else
            nptf = 0
            nface = 0
        end if
!
    else
!       PROBLEME DE DIMENSION : NI 2D, NI 3D
        ASSERT(ndim .eq. 2)
    end if
!
999 continue
end subroutine
