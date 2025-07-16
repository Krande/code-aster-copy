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
subroutine rcjaco(ar, valpro)
    implicit none
#include "asterf_types.h"
#include "asterfort/utmess.h"
    real(kind=8) :: ar(*), valpro(3)
! BUT : ROUTINE SIMPLIFIEE DE JACOBI POUR PERFORMANCE
!
! ----------------------------------------------------------------------
    integer(kind=8) :: nperm, i, ii, niter, j, jp1, jm1, ljk, jj, kp1, km1, jk, kk, im1
    integer(kind=8) :: ij, ik, lji, lki, ji, ki, k
    real(kind=8) :: tol, toldyn, valaux(3), eps, akk, ajj, ab, verif, br(6)
    real(kind=8) :: eptola, epcoma, eptolb, epcomb, raci, d1, d2, den, ca, cg
    real(kind=8) :: aj, bj, ak, bk, rtol, dif, epsa, compa, epsb, compb
    aster_logical :: iconv
    data nperm, tol, toldyn/12, 1.d-10, 1.d-2/
!
! ----------------------------------------------------------------------
!
!     ---       INITIALISATION DE LA MATRICE DE MASSE      ---
!
    br(1) = 1.d0
    br(2) = 0.d0
    br(3) = 0.d0
    br(4) = 1.d0
    br(5) = 0.d0
    br(6) = 1.d0
!
! ----------------------------------------------------------------------
!
!     ---       INITIALISATION DES VALEURS PROPRES      ---
!     --- TERME DIAGONAL RAIDEUR / TERME DIAGONAL MASSE ---
!
!
    ii = 1
    do i = 1, 3
        if (br(ii) .eq. 0.0d0) then
            call utmess('F', 'ALGELINE4_19')
        end if
        valaux(i) = ar(ii)/br(ii)
        valpro(i) = valaux(i)
        ii = ii+3+1-i
    end do
!
!     ------------------------------------------------------------------
!     ------------------- ALGORITHME DE JACOBI -------------------------
!     ------------------------------------------------------------------
!
    niter = 0
!
30  continue
!
    niter = niter+1
    eps = (toldyn**niter)**2
!
!     --- BOUCLE SUR LES LIGNES ---
    do j = 1, 3-1
        jp1 = j+1
        jm1 = j-1
        ljk = jm1*3-jm1*j/2
        jj = ljk+j
!        ---- BOUCLE SUR LES COLONNES ---
        do k = jp1, 3
            kp1 = k+1
            km1 = k-1
            jk = ljk+k
            kk = km1*3-km1*k/2+k
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
                d1 = ab*0.5d0+raci
                d2 = ab*0.5d0-raci
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
            if (jm1-1 .ge. 0) then
                do i = 1, jm1
                    im1 = i-1
                    ij = im1*3-im1*i/2+j
                    ik = im1*3-im1*i/2+k
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
            if (kp1-3 .le. 0) then
                lji = jm1*3-jm1*j/2
                lki = km1*3-km1*k/2
                do i = kp1, 3
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
                lji = jm1*3-jm1*j/2
                do i = jp1, km1
                    ji = lji+i
                    im1 = i-1
                    ik = im1*3-im1*i/2+k
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
            ak = ar(kk)
            bk = br(kk)
            ar(kk) = ak+2.0d0*ca*ar(jk)+ca*ca*ar(jj)
            br(kk) = bk+2.0d0*ca*br(jk)+ca*ca*br(jj)
            ar(jj) = ar(jj)+2.0d0*cg*ar(jk)+cg*cg*ak
            br(jj) = br(jj)+2.0d0*cg*br(jk)+cg*cg*bk
            ar(jk) = 0.0d0
            br(jk) = 0.0d0
!
41          continue
        end do
    end do
!
!     --- CALCUL DES NOUVELLES VALEURS PROPRES ---
!
    ii = 1
    do i = 1, 3
        if (br(ii) .eq. 0.0d0) then
            call utmess('F', 'ALGELINE4_19')
        end if
        valpro(i) = ar(ii)/br(ii)
        ii = ii+3+1-i
    end do
!
!     --- TEST DE CONVERGENCE SUR LES VALEURS PROPRES ---
!
    iconv = .true.
    do i = 1, 3
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
    do j = 1, 3-1
        jm1 = j-1
        jp1 = j+1
        ljk = jm1*3-jm1*j/2
        jj = ljk+j
        do k = jp1, 3
            km1 = k-1
            jk = ljk+k
            kk = km1*3-km1*k/2+k
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
        do i = 1, 3
            valaux(i) = valpro(i)
        end do
!
!        --- TEST SUR LE NOMBRE D'ITERATIONS ---
!
        if (niter .lt. nperm) goto 30
!
    end if
!
end subroutine
