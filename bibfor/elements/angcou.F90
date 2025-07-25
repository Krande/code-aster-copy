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
subroutine angcou(coor, zk1, izk, icoude, zk2, &
                  rayon, theta, angl1, angl2, angl3, &
                  pgl1, pgl2, pgl3, omega, dn1n2, &
                  epsi, crit, zk3)
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/angvxy.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    real(kind=8) :: coor(9), rayon, theta, epsi, t1(3), t2(3), coo1(3), coo2(3)
    real(kind=8) :: coo3(3), norme1, norme2, normez, x3(3), x1(3), x2(3), ct2
    real(kind=8) :: st2
    real(kind=8) :: y1(3), y2(3), y3(3), pgl1(3, 3), pgl2(3, 3), pgl3(3, 3)
    real(kind=8) :: a(3)
    real(kind=8) :: angl1(3), angl2(3), angl3(3), cosome, sinome, nx1, zcoud(3)
    real(kind=8) :: zk1(3), zk2(3), axe(3), zzk1(3), zzzk1(3), zkini(3), zk3(3)
    real(kind=8) :: costet, ct, st, omega, dn1n2, nx2, omega2, norme3, epsi2
    real(kind=8) :: themax
    real(kind=8) :: psca, valr(2)
    character(len=8) :: crit
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DE LA GEOMETRIE TUYAU DROIT OU COURBE
!    - ARGUMENTS
!        DONNEES:      COOR       -->  CORDONNEES DES NOEUDS
!                      ZK1        -->  VECTEUR ZK A UN NOEUD EXTREMITE
!                      IZK        -->  INDICE DE PARCOURS
!                      EPSI       -->  PRECISION POUR LE CALCUL
!                      CRIT       -->  RELATIF OU ABSOLU
!         SORTIE:      ICOUDE     -->  =0 DROIT =1 COUDE
!                      ZK2        -->  VECTEUR ZK A L'AUTRE EXTREMITE
!                      RAYON      -->  RAYON DU COUDE SI C'EN EST UN
!                      THETA      -->  ANGLE DU COUDE
!                      ANGL1,2,3  -->  ANGLES NAUTIQUES EN CHAQUE NOEUD
!                      PGL1,2,3   -->  MATRICE DE CHAGMENT DE REPERE
!                      OMEGA      -->  ANGLE ENTRE N ET LA GENERATRICE
!                      DN1N2      -->  DN1N2  DISTANCE ENTRE EXTREMITES
!                      ZK2        -->  VECTEUR ZK AU NOEUD MILIEU
! ......................................................................
!
    integer(kind=8) :: icoude, i, izk
    real(kind=8) :: test
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
!
!     NOEUD 3 = NOEUD MILIEU
    do i = 1, 3
        angl1(i) = 0.d0
        angl2(i) = 0.d0
        angl3(i) = 0.d0
        coo1(i) = coor(i)
        coo2(i) = coor(3+i)
        coo3(i) = coor(6+i)
        zkini(i) = zk1(i)
    end do
!
!     POUR VERIFICATIONS (PAS TRES EXIGEANTES) SUR LA GEOMETRIE
!     SINON, IL FAUDRAIT INTRODUIRE UN AUTRE MOT CLE PRECISON2
!     DIFFERENT DE PRECISION
!
    epsi2 = 1.d-4
    themax = r8pi()/8.d0*(1.d0+epsi2)
!
    dn1n2 = sqrt((coo1(1)-coo2(1))**2+(coo1(2)-coo2(2))**2+(coo1(3)-coo2(3))**2)
!
    t1 = coo3-coo1
    t2 = coo2-coo3
    call normev(t1, norme1)
    call normev(t2, norme2)
    call provec(t2, t1, zcoud)
    call normev(zcoud, normez)
!
!     VERIF QUE LE NOEUD MILIEU EST BIEN LE TROISIEME
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    psca = ddot(b_n, t2, b_incx, t1, b_incy)
    if (psca .le. 0.d0) then
        call utmess('F', 'ELEMENTS_5')
    end if
!
!     EPSI EST CELUI DONNE PAR LE MOT CLE PRECISION
!
    if (crit .eq. 'RELATIF') then
!     CONTRAIREMENT AUX APPARENCES C'EST UN CRITERE RELATIF
!     PUISQUE T1 ET T2 SONT NORMES
        test = epsi
    else if (crit .eq. 'ABSOLU') then
        test = epsi/norme1/norme2
    end if
!
!     TUYAU DROIT OU COURBE ?
!
    if (normez .le. test) then
        icoude = 0
        rayon = 0.d0
        theta = 0.d0
        omega = 0.d0
        x1 = coo2-coo1
!
!        ON VEUT UN ZK1 PERPENDICULAIRE A X1
!        ON PROJETTE ZK1 SUR LE PLAN NORMAL A X1
!
        call provec(x1, zkini, a)
        call normev(a, norme1)
!
!        TEST DE NON COLINEARITE
!
        call normev(zkini, norme3)
        test = epsi2*dn1n2*norme3
!
        if (norme1 .le. test) then
            call utmess('F', 'ELEMENTS_6')
        end if
!
        call provec(a, x1, zk1)
        call normev(zk1, norme2)
        do i = 1, 3
            zk2(i) = zk1(i)
            zk3(i) = zk1(i)
        end do
        call provec(x1, zk1, y1)
        call angvxy(x1, y1, angl1)
        do i = 1, 3
            angl2(i) = angl1(i)
            angl3(i) = angl1(i)
        end do
!
        call matrot(angl1, pgl1)
        call matrot(angl2, pgl2)
        call matrot(angl3, pgl3)
!
    else
        icoude = 1
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        costet = ddot(b_n, t1, b_incx, t2, b_incy)
        theta = 2.d0*atan2(normez, costet)
        if (theta .gt. themax) then
            valr(1) = theta
            valr(2) = themax
            call utmess('A', 'ELEMENTS_7', nr=2, valr=valr)
        end if
        rayon = dn1n2/2.d0/normez
!        CALCUL DES REPERES LOCAUX EN CHAQUE NOEUD
        x3 = coo2-coo1
        call provec(x3, zcoud, y3)
        ct2 = cos(theta/2.d0)
        st2 = sin(theta/2.d0)
        do i = 1, 3
            x1(i) = x3(i)*ct2-y3(i)*st2
            x2(i) = x3(i)*ct2+y3(i)*st2
        end do
        call provec(x1, zcoud, y1)
        call provec(x2, zcoud, y2)
        call angvxy(x1, y1, angl1)
        call angvxy(x2, y2, angl2)
        call angvxy(x3, y3, angl3)
        call matrot(angl1, pgl1)
        call matrot(angl2, pgl2)
        call matrot(angl3, pgl3)
!
!        CALCUL DE L'ANGLE OMEGA ENTRE N ET LA GENERATRICE
!        SI ON CONSTRUIT LA GENERATRICE DANS LE SENS DE
!        DESCRIPTION DU MAILLAGE, ON FAIT UNE ROTATION
!        DE THETA AUTOUR DE -Z. SINON, AUTOUR DE Z
!
        if (izk .lt. 0) then
            do i = 1, 3
                axe(i) = zcoud(i)
            end do
        else
            do i = 1, 3
                axe(i) = -zcoud(i)
            end do
        end if
!
!        ON VEUT UN ZK1 PERPENDICULAIRE A X1
!        ON PROJETTE ZK1 SUR LE PLAN NORMAL A X1
!
        if (izk .gt. 0) then
            call provec(x1, zkini, a)
            call normev(a, norme1)
            call provec(a, x1, zk1)
            call normev(zk1, norme2)
        else
            call provec(x2, zkini, a)
            call normev(a, norme1)
            call provec(a, x2, zk1)
            call normev(zk1, norme2)
        end if
!
!        TEST DE NON COLINEARITE
!
        test = epsi2*dn1n2*norme2
        if (norme1 .le. test) then
            call utmess('F', 'ELEMENTS_6')
        end if
!
        call provec(axe, zk1, zzk1)
        call provec(axe, zzk1, zzzk1)
        ct = cos(theta)
        st = sin(theta)
!
        do i = 1, 3
            zk2(i) = zk1(i)+st*zzk1(i)+(1.d0-ct)*zzzk1(i)
        end do
        do i = 1, 3
            zk3(i) = zk1(i)+st2*zzk1(i)+(1.d0-ct2)*zzzk1(i)
        end do
        call normev(zk1, norme1)
        call normev(zk2, norme2)
        call normev(zk3, norme3)
!
!        OMEGA ANGLE ENTRE Z ET ZK1 ET AUSSI ENTRE Z ET ZK2
!
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cosome = ddot(b_n, zcoud, b_incx, zk1, b_incy)
        call provec(zk1, zcoud, a)
        if (izk .lt. 0) then
            call normev(x2, nx2)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinome = ddot(b_n, x2, b_incx, a, b_incy)
        else
            call normev(x1, nx1)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinome = ddot(b_n, x1, b_incx, a, b_incy)
        end if
        omega2 = atan2(sinome, cosome)
!
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        cosome = ddot(b_n, zcoud, b_incx, zk2, b_incy)
        call provec(zk2, zcoud, a)
        if (izk .gt. 0) then
            call normev(x2, nx2)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinome = ddot(b_n, x2, b_incx, a, b_incy)
        else
            call normev(x1, nx1)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinome = ddot(b_n, x1, b_incx, a, b_incy)
        end if
        omega = atan2(sinome, cosome)
!
        test = 0.d0
        if (crit .eq. 'RELATIF') then
            if (abs(omega) .lt. r8prem()) then
                test = r8prem()
            else
                test = epsi2*abs(omega)
            end if
        else if (crit .eq. 'ABSOLU') then
            test = epsi2
        end if
!
        if (abs(omega2-omega) .gt. test) then
            valr(1) = omega
            valr(2) = omega2
            call utmess('F', 'ELEMENTS_1', nr=2, valr=valr)
        end if
    end if
end subroutine
