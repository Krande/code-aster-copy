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
subroutine jacobi(nbvec, nperm, tol, toldyn, ar, &
                  br, vecpro, valpro, valaux, nitjac, &
                  type, iordre)
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/vpordo.h"
    integer(kind=8) :: nbvec, nperm, nitjac, type, iordre
    real(kind=8) :: ar(nbvec), br(nbvec), vecpro(nbvec, nbvec), valpro(nbvec)
    real(kind=8) :: valaux(nbvec), tol, toldyn
!
! BUT : RESOLUTION DU PROBLEME REDUIT AUX VALEURS PROPRES PAR
!       LA DECOMPOSITION DE JACOBI GENERALISEE
!
!     IN  : NBVEC  : NOMBRE DE VECTEURS
!     IN  : NPERM  : NOMBRE MAX D'ITERATIONS DE LA METHODE DE JACOBI
!     IN  : TOL    : PRECISION DE CONVERGENCE
!     IN  : TOLDYN : PRECISION DE PETITESSE DYNAMIQUE
!     IN / OUT  : AR : MATRICE DE RAIDEUR PROJETEE
!     IN / OUT  : BR : MATRICE DE MASSE PROJETEE
!                      CES DEUX MATRICES SONT SYMETRIQUES ET ON A
!                      DONC STOCKE SOUS FORME D'UN VECTEUR SEULEMENT
!                      LA MOITIE : A11 A12  A1N A22  A2N
!         ATTENTION : CES DEUX MATRICES SONT MODIFIEES ET DONC
!                     INUTILISABLES EN SORTIE.
!     OUT : VECPRO : VECTEURS PROPRES DE L ITERATION
!     OUT : VALPRO : VALEURS PROPRES DE L ITERATION
!     OUT : NITJAC : NOMBRE D'ITERATIONS DE JACOBI
!     IN  : TYPE  : TYPE DU TRI SUR LES VALEURS PROPRES (ET VECT. PRO.)
!        * SI TYPE = 0  TRI EN VALEUR RELATIVE
!        * SI TYPE = 1  TRI EN VALEUR ABSOLUE
!        * SI TYPE = 2  PAS DE TRI
! IN  IORDRE : IS : ORDRE DU TRI SUR LES VALEURS.
!        * SI IORDRE = 0  TRI PAR ORDRE CROISSANT
!        * SI IORDRE = 1  TRI PAR ORDRE DECROISSANT
!        * SI IORDRE = 2  PAS DE TRI
! ----------------------------------------------------------------------
    aster_logical :: iconv
! ----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!     ------- INITIALISATION DES VALEURS ET VECTEURS PROPRES -----------
!     ------------------------------------------------------------------
!
!     ---       INITIALISATION DES VALEURS PROPRES      ---
!     --- TERME DIAGONAL RAIDEUR / TERME DIAGONAL MASSE ---
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ii, ij, ik, im1, j, ji
    integer(kind=8) :: jj, jk, jm1, jp1, k, ki, kk
    integer(kind=8) :: km1, kp1, lji, ljk, lki, niter
    real(kind=8) :: ab, aj, ajj, ak, akk, bb, bj
    real(kind=8) :: bk, ca, cg, compa, compb, d1, d2
    real(kind=8) :: den, dif, epcoma, epcomb, eps, epsa, epsb
    real(kind=8) :: eptola, eptolb, raci, rtol, verif, xj, xk
!
!-----------------------------------------------------------------------
    ii = 1
    do i = 1, nbvec
        if (br(ii) .eq. 0.0d0) then
            call utmess('F', 'ALGELINE4_19')
        end if
        valaux(i) = ar(ii)/br(ii)
        valpro(i) = valaux(i)
        ii = ii+nbvec+1-i
    end do
!
!     --- INITIALISATION DES VECTEURS PROPRES (MATRICE IDENTITE) ---
!
!
    vecpro(:, :) = 0.d0
!
    do i = 1, nbvec
        vecpro(i, i) = 1.0d0
    end do
!
!     ------------------------------------------------------------------
!     ------------------- ALGORITHME DE JACOBI -------------------------
!     ------------------------------------------------------------------
!
    niter = 0
    nitjac = 0
!
30  continue
!
    niter = niter+1
    eps = (toldyn**niter)**2
!
!     --- BOUCLE SUR LES LIGNES ---
    do j = 1, nbvec-1
        jp1 = j+1
        jm1 = j-1
        ljk = jm1*nbvec-jm1*j/2
        jj = ljk+j
!        ---- BOUCLE SUR LES COLONNES ---
        do k = jp1, nbvec
            kp1 = k+1
            km1 = k-1
            jk = ljk+k
            kk = km1*nbvec-km1*k/2+k
!           --- CALCUL DES COEFFICIENTS DE LA ROTATION DE GIVENS ---
            eptola = abs((ar(jk)*ar(jk)))
            epcoma = abs(ar(jj)*ar(kk))*eps
            eptolb = abs((br(jk)*br(jk)))
            epcomb = abs(br(jj)*br(kk))*eps
            if ((eptola .eq. 0.d0) .and. (eptolb .eq. 0.d0)) goto 41
            if ((eptola .le. epcoma) .and. (eptolb .le. epcomb)) goto 41
            akk = ar(kk)*br(jk)-br(kk)*ar(jk)
            ajj = ar(jj)*br(jk)-br(jj)*ar(jk)
            ab = ar(jj)*br(kk)-ar(kk)*br(jj)
            verif = (ab*ab+4.0d0*akk*ajj)/4.0d0
            if (verif .ge. 0.0d0) then
                raci = sqrt(verif)
                d1 = ab/2.0d0+raci
                d2 = ab/2.0d0-raci
            else
                goto 41
            end if
            den = d1
            if (abs(d2) .gt. abs(d1)) den = d2
            if (den .eq. 0.0d0) then
                ca = 0.d0
                cg = -ar(jk)/ar(kk)
            else
                ca = akk/den
                cg = -ajj/den
            end if
!           --- TRANSFORMATION DES MATRICES DE RAIDEUR ET DE MASSE ---
            if (nbvec-2 .ne. 0) then
                if (jm1-1 .ge. 0) then
                    do i = 1, jm1
                        im1 = i-1
                        ij = im1*nbvec-im1*i/2+j
                        ik = im1*nbvec-im1*i/2+k
                        aj = ar(ij)
                        bj = br(ij)
                        ak = ar(ik)
                        bk = br(ik)
                        ar(ij) = aj+cg*ak
                        br(ij) = bj+cg*bk
                        ar(ik) = ak+ca*aj
                        br(ik) = bk+ca*bj
                    end do
                end if
                if (kp1-nbvec .le. 0) then
                    lji = jm1*nbvec-jm1*j/2
                    lki = km1*nbvec-km1*k/2
                    do i = kp1, nbvec
                        ji = lji+i
                        ki = lki+i
                        aj = ar(ji)
                        bj = br(ji)
                        ak = ar(ki)
                        bk = br(ki)
                        ar(ji) = aj+cg*ak
                        br(ji) = bj+cg*bk
                        ar(ki) = ak+ca*aj
                        br(ki) = bk+ca*bj
                    end do
                end if
                if (jp1-km1 .le. 0) then
                    lji = jm1*nbvec-jm1*j/2
                    do i = jp1, km1
                        ji = lji+i
                        im1 = i-1
                        ik = im1*nbvec-im1*i/2+k
                        aj = ar(ji)
                        bj = br(ji)
                        ak = ar(ik)
                        bk = br(ik)
                        ar(ji) = aj+cg*ak
                        br(ji) = bj+cg*bk
                        ar(ik) = ak+ca*aj
                        br(ik) = bk+ca*bj
                    end do
                end if
            end if
            ak = ar(kk)
            bk = br(kk)
            ar(kk) = ak+2.0d0*ca*ar(jk)+ca*ca*ar(jj)
            br(kk) = bk+2.0d0*ca*br(jk)+ca*ca*br(jj)
            ar(jj) = ar(jj)+2.0d0*cg*ar(jk)+cg*cg*ak
            br(jj) = br(jj)+2.0d0*cg*br(jk)+cg*cg*bk
            ar(jk) = 0.0d0
            br(jk) = 0.0d0
!           --- TRANSFORMATION DES VECTEURS PROPRES ---
            do i = 1, nbvec
                xj = vecpro(i, j)
                xk = vecpro(i, k)
                vecpro(i, j) = xj+cg*xk
                vecpro(i, k) = xk+ca*xj
            end do
!
41          continue
        end do
    end do
!
!     --- CALCUL DES NOUVELLES VALEURS PROPRES ---
!
    ii = 1
    do i = 1, nbvec
        if (br(ii) .eq. 0.0d0) then
            call utmess('F', 'ALGELINE4_19')
        end if
        valpro(i) = ar(ii)/br(ii)
        ii = ii+nbvec+1-i
    end do
!
!     --- TEST DE CONVERGENCE SUR LES VALEURS PROPRES ---
!
    iconv = .true.
    do i = 1, nbvec
        rtol = tol*valaux(i)
        dif = abs(valpro(i)-valaux(i))
        if (dif .gt. abs(rtol)) then
            iconv = .false.
            goto 9998
        end if
    end do
!
!     ---    CALCUL DES FACTEURS DE COUPLAGE   ---
!     --- TEST DE CONVERGENCE SUR CES FACTEURS ---
!
    eps = tol**2
    do j = 1, nbvec-1
        jm1 = j-1
        jp1 = j+1
        ljk = jm1*nbvec-jm1*j/2
        jj = ljk+j
        do k = jp1, nbvec
            km1 = k-1
            jk = ljk+k
            kk = km1*nbvec-km1*k/2+k
            epsa = abs(ar(jk)*ar(jk))
            compa = eps*abs(ar(jj)*ar(kk))
            epsb = abs(br(jk)*br(jk))
            compb = eps*abs(br(jj)*br(kk))
            if (epsa .ge. compa .or. epsb .ge. compb) then
                iconv = .false.
                goto 9998
            end if
        end do
    end do
!
9998 continue
!
!     ---  SI ON N'A PAS CONVERGE ---
!
    if (.not. iconv) then
!
!        --- TRANSLATION DES VALEURS PROPRES ---
!
        do i = 1, nbvec
            valaux(i) = valpro(i)
        end do
!
!        --- TEST SUR LE NOMBRE D'ITERATIONS ---
!
        if (niter .lt. nperm) goto 30
!
    end if
!
!     --- SI CONVERGENCE OU NOMBRE MAX D'ITERATIONS ATTEINT ---
!     ---         MISE A JOUR DES VECTEURS PROPRES          ---
!
    ii = 1
    do i = 1, nbvec
        if (br(ii) .ge. 0.0d0) then
            bb = sqrt(br(ii))
        else
            bb = -sqrt(abs(br(ii)))
        end if
        do k = 1, nbvec
            vecpro(k, i) = vecpro(k, i)/bb
        end do
        ii = ii+nbvec+1-i
    end do
!
    nitjac = niter
!
!     --- CLASSEMENT DES MODES
    if ((type .eq. 2) .or. (iordre .eq. 2)) then
    else
        call vpordo(type, iordre, nbvec, valpro, vecpro, &
                    nbvec)
    end if
end subroutine
