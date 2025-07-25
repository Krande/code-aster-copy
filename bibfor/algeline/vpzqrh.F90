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
subroutine vpzqrh(h, neq, ih, k, l, &
                  wr, wi, z, iz, mxiter, &
                  ier, nitqr)
    implicit none
#include "asterf_types.h"
#include "asterc/r8prem.h"
    integer(kind=8) :: neq, ih, k, l, iz, ier, nitqr
    real(kind=8) :: h(ih, neq), wr(neq), wi(neq), z(iz, neq)
!     RECHERCHE DES VALEURS PROPRES PAR LA METHODE QR SUR UNE MATRICE
!     MISE SOUS LA FORME DE HESSENBERG
!     ------------------------------------------------------------------
!     REFERENCE: F.L. BAUER - J.H. WILKINSON - C. REINSCH
!        HANDBOOK FOR AUTOMATIC COMPUTATION - LINEAR ALGEBRA - VOL.2
!        PAGE XXX
!     ------------------------------------------------------------------
    integer(kind=8) :: i, ien, ienm2, npl, ll, lb, naml, mm, m, mp2, ka, na
    integer(kind=8) :: iter, j, jj
    real(kind=8) :: epsmac, t, x, y, w, s, zz
    real(kind=8) :: r, p, q, rnorm, ra, sa, scale, vr, vi
    complex(kind=8) :: z3
    aster_logical :: notlas
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ii, mxiter, nn
!-----------------------------------------------------------------------
    ier = 0
    epsmac = r8prem()
    nitqr = 0
!
!     --- STOCKER LES RACINES TROUVEES PAR VPZBAL (DIT LA BALANCE) ---
    rnorm = 0.0d0
    ka = 1
    do i = 1, neq
        do j = ka, neq
            rnorm = rnorm+abs(h(i, j))
        end do
        ka = i
        if (i .ge. k .and. i .le. l) goto 10
        wr(i) = h(i, i)
        wi(i) = 0.d0
10      continue
    end do
    ien = l
    t = 0.d0
!
!     --- RECHERCHE DES VALEURS PROPRES SUIVANTES ---
15  continue
    if (ien .lt. k) goto 145
    iter = 0
    na = ien-1
    ienm2 = na-1
!
!     --- RECHERCHE DU PLUS PETIT ELEMENT (SIMPLE) SUR LA SUR-DIAGONALE
20  continue
    npl = ien+k
    do ll = k, ien
        lb = npl-ll
        if (lb .eq. k) goto 30
        s = abs(h(lb-1, lb-1))+abs(h(lb, lb))
        if (s .eq. 0.0d0) s = rnorm
        if (abs(h(lb, lb-1)) .le. epsmac*s) goto 30
    end do
!
30  continue
    x = h(ien, ien)
    if (lb .eq. ien) goto 110
    y = h(na, na)
    w = h(ien, na)*h(na, ien)
    if (lb .eq. na) goto 115
    if (iter .eq. mxiter) goto 250
!
!     --- FORMER UN DECALAGE TOUT LES DIX COUPS ---
    iter = iter+1
    if (iter .gt. nitqr) then
        nitqr = iter
    end if
    if (mod(iter, 10) .eq. 0) then
        t = t+x
        do i = k, ien
            h(i, i) = h(i, i)-x
        end do
        s = abs(h(ien, na))+abs(h(na, ienm2))
        x = 0.75d0*s
        y = x
        w = -0.4375d0*s*s
    end if
!
!     --- RECHERCHE DES 2 PLUS PETITS ELEMENTS SUR LA SUR-DIAGONALE
    naml = ienm2+lb
    do mm = lb, ienm2
        m = naml-mm
        zz = h(m, m)
        r = x-zz
        s = y-zz
        p = (r*s-w)/h(m+1, m)+h(m, m+1)
        q = h(m+1, m+1)-zz-r-s
        r = h(m+2, m+1)
        s = abs(p)+abs(q)+abs(r)
        p = p/s
        q = q/s
        r = r/s
        if (m .eq. lb) goto 50
        if (abs(h(m, m-1))*(abs(q)+abs(r)) .le. &
            epsmac*abs(p)*(abs(h(m-1, m-1))+abs(zz)+abs(h(m+1, m+1)))) goto 50
    end do
50  continue
    mp2 = m+2
    do i = mp2, ien
        h(i, i-2) = 0.d0
        if (i .eq. mp2) goto 55
        h(i, i-3) = 0.d0
55      continue
    end do
!
!     EN AVANT POUR LE "DOUBLE QR" SUR LA SOUS MATRICE
!              LIGNES DE "L" A "EN" ET COLONNES DE "M" A "EN"
    do ka = m, na
        notlas = ka .ne. na
        if (ka .eq. m) goto 60
        p = h(ka, ka-1)
        q = h(ka+1, ka-1)
        r = 0.d0
        if (notlas) r = h(ka+2, ka-1)
        x = abs(p)+abs(q)+abs(r)
        if (x .eq. 0.d0) goto 105
        p = p/x
        q = q/x
        r = r/x
60      continue
        s = sign(sqrt(p*p+q*q+r*r), p)
        if (ka .eq. m) then
            if (lb .ne. m) h(ka, ka-1) = -h(ka, ka-1)
        else
            h(ka, ka-1) = -s*x
        end if
        p = p+s
        x = p/s
        y = q/s
        zz = r/s
        q = q/p
        r = r/p
!        --- ALTERATION DES LIGNES ---
        do j = ka, neq
            p = h(ka, j)+q*h(ka+1, j)
            if (notlas) then
                p = p+r*h(ka+2, j)
                h(ka+2, j) = h(ka+2, j)-p*zz
            end if
            h(ka+1, j) = h(ka+1, j)-p*y
            h(ka, j) = h(ka, j)-p*x
        end do
        j = min(ien, ka+3)
!        --- ALTERATION DES COLONNES ---
        do i = 1, j
            p = x*h(i, ka)+y*h(i, ka+1)
            if (.not. notlas) goto 85
            p = p+zz*h(i, ka+2)
            h(i, ka+2) = h(i, ka+2)-p*r
85          continue
            h(i, ka+1) = h(i, ka+1)-p*q
            h(i, ka) = h(i, ka)-p
        end do
        if (iz .ge. neq) then
!           --- ON GARDE LES TRANSFORMATIONS POUR LES VECTEURS ----
            do i = k, l
                p = x*z(i, ka)+y*z(i, ka+1)
                if (notlas) then
                    p = p+zz*z(i, ka+2)
                    z(i, ka+2) = z(i, ka+2)-p*r
                end if
                z(i, ka+1) = z(i, ka+1)-p*q
                z(i, ka) = z(i, ka)-p
            end do
        end if
105     continue
    end do
    goto 20
!
!     ---------------------------- AU CAS PAR CAS  ---------------------
!
!     --- UNE RACINE TROUVEE ---
110 continue
    h(ien, ien) = x+t
    wr(ien) = h(ien, ien)
    wi(ien) = 0.d0
    ien = na
    goto 15
!
!     --- DEUX RACINES TROUVEES ---
115 continue
    p = (y-x)*0.5d0
    q = p*p+w
    zz = sqrt(abs(q))
    h(ien, ien) = x+t
    x = h(ien, ien)
    h(na, na) = y+t
    if (q .lt. 0.d0) goto 135
!
!     --- RACINES DOUBLES REELLES ---
    zz = p+sign(zz, p)
    wr(na) = x+zz
    wr(ien) = wr(na)
    if (zz .ne. 0.d0) wr(ien) = x-w/zz
    wi(na) = 0.d0
    wi(ien) = 0.d0
    x = h(ien, na)
!
!     --- SI X ET ZZ TROP PETIT ALORS ON FAIT UNE NORMALISATION   ---
!     --- MISE A L'ECHELLE OU RECADRAGE C'EST SELON VOTRE CULTURE ---
    scale = abs(x)+abs(zz)
    r = scale*sqrt((x/scale)**2+(zz/scale)**2)
    p = x/r
    q = zz/r
!
!     --- ALTERATION DES LIGNES ---
    do j = na, neq
        zz = h(na, j)
        h(na, j) = q*zz+p*h(ien, j)
        h(ien, j) = q*h(ien, j)-p*zz
    end do
!
!     --- ALTERATION DES COLONNES ---
    do i = 1, ien
        zz = h(i, na)
        h(i, na) = q*zz+p*h(i, ien)
        h(i, ien) = q*h(i, ien)-p*zz
    end do
!
    if (iz .ge. neq) then
!        --- STOCKER LES TRANSFORMATIONS POUR LES VECTEURS ---
        do i = k, l
            zz = z(i, na)
            z(i, na) = q*zz+p*z(i, ien)
            z(i, ien) = q*z(i, ien)-p*zz
        end do
    end if
    goto 140
!
!     --- VALEURS COMPLEXES CONJUGUEES ---
135 continue
    wr(na) = x+p
    wr(ien) = x+p
    wi(na) = zz
    wi(ien) = -zz
140 continue
    ien = ienm2
    goto 15
!
!
!     --- TOUTES LES RACINES SONT TROUVEES, ON DEBUTE LA REMONTEE ---
145 continue
    if (iz .lt. neq) goto 999
!
!     ---- ON S'OCCUPE MAINTENANT DES VECTEURS ---
!
    if (rnorm .eq. 0.d0) goto 999
    do nn = 1, neq
        ien = neq+1-nn
        p = wr(ien)
        q = wi(ien)
        na = ien-1
        if (q .gt. 0.d0) goto 220
        if (q .lt. 0.d0) goto 180
!
!        --- VECTEUR REEL ---
        m = ien
        h(ien, ien) = 1.d0
        if (na .eq. 0) goto 220
        do ii = 1, na
            i = ien-ii
            w = h(i, i)-p
            r = h(i, ien)
            do j = m, na
                r = r+h(i, j)*h(j, ien)
            end do
            if (wi(i) .ge. 0.d0) goto 160
            zz = w
            s = r
            goto 175
160         continue
            m = i
            if (wi(i) .ne. 0.d0) goto 165
            t = w
            if (w .eq. 0.d0) t = epsmac*rnorm
            h(i, ien) = -r/t
            goto 175
!
!           RESOLUTION DANS LE CAS REEL ---
165         continue
            x = h(i, i+1)
            y = h(i+1, i)
            q = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)
            t = (x*s-zz*r)/q
            h(i, ien) = t
            if (abs(x) .le. abs(zz)) then
                h(i+1, ien) = (-s-y*t)/zz
            else
                h(i+1, ien) = (-r-w*t)/x
            end if
175         continue
        end do
        goto 220
!
!        --- CAS OU LE DERNIER VECTEUR EST IMAGINAIRE ---
180     continue
        m = na
!        --- VECTEUR COMPLEXE ---
        if (abs(h(ien, na)) .le. abs(h(na, ien))) goto 185
        h(na, na) = q/h(ien, na)
        h(na, ien) = -(h(ien, ien)-p)/h(ien, na)
        goto 190
185     continue
        z3 = dcmplx(0.d0, -h(na, ien))/dcmplx(h(na, na)-p, q)
        h(na, na) = dble(z3)
        h(na, ien) = dimag(z3)
190     continue
        h(ien, na) = 0.d0
        h(ien, ien) = 1.d0
        ienm2 = na-1
        if (ienm2 .eq. 0) goto 220
        do ii = 1, ienm2
            i = na-ii
            w = h(i, i)-p
            ra = 0.d0
            sa = h(i, ien)
            do j = m, na
                ra = ra+h(i, j)*h(j, na)
                sa = sa+h(i, j)*h(j, ien)
            end do
            if (wi(i) .lt. 0.d0) then
                zz = w
                r = ra
                s = sa
            else if (wi(i) .eq. 0.d0) then
                m = i
                z3 = dcmplx(-ra, -sa)/dcmplx(w, q)
                h(i, na) = dble(z3)
                h(i, ien) = dimag(z3)
            else
!
!              --- RESOUDRE LES EQUATIONS (EN COMPLEXE)
                m = i
                x = h(i, i+1)
                y = h(i+1, i)
                vr = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)-q*q
                vi = (wr(i)-p)*q
                vi = vi+vi
                if (vr .eq. 0.d0 .and. vi .eq. 0.d0) vr = epsmac*rnorm*( &
                                                          abs(w)+abs(q)+abs(x)+abs(y)+abs(zz))
                z3 = dcmplx(x*r-zz*ra+q*sa, x*s-zz*sa-q*ra)/dcmplx(vr, vi)
                h(i, na) = dble(z3)
                h(i, ien) = dimag(z3)
                if (abs(x) .le. abs(zz)+abs(q)) then
                    z3 = dcmplx(-r-y*h(i, na), -s-y*h(i, ien))/dcmplx(zz, q)
                    h(i+1, na) = dble(z3)
                    h(i+1, ien) = dimag(z3)
                else
                    h(i+1, na) = (-ra-w*h(i, na)+q*h(i, ien))/x
                    h(i+1, ien) = (-sa-w*h(i, ien)-q*h(i, na))/x
                end if
            end if
        end do
220     continue
    end do
!     --- FIN DE LA REMONTEE ---
!
!     --- VECTEURS DES RACINES ISOLEES ---
    do i = 1, neq
        if (i .ge. k .and. i .le. l) goto 230
        do j = i, neq
            z(i, j) = h(i, j)
        end do
230     continue
    end do
    if (l .eq. 0) goto 999
!
!     APPLICATION DES TRANSFORMATIONS ---
    do jj = k, neq
        j = neq+k-jj
        m = min(j, l)
        do i = k, l
            zz = 0.d0
            do ka = k, m
                zz = zz+z(i, ka)*h(ka, j)
            end do
            z(i, j) = zz
        end do
    end do
    goto 999
!
!     PAS DE CONVERGEBCE APRES "MXITER" INTERATION
!         ==>  IER = L'INDICE DE LA VALEUR PROPRE COURANTE
!
250 continue
    ier = ien
    do i = 1, ien
        wr(i) = 0.d0
        wi(i) = 0.d0
    end do
    if (iz .ge. neq) then
        do i = 1, neq
            do j = 1, neq
                z(i, j) = 0.d0
            end do
        end do
    end if
!     --- SORTIE ---
999 continue
end subroutine
