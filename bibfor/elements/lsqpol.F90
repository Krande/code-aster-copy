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
subroutine lsqpol(ordre, e1, npt, xx, yy, &
                  ordok, poly, sigma)
! aslint: disable=W1306
    implicit none
!
! ------------------------------------------------------------------
!
!              REGRESSION POLYNOMIALE DE TYPE MOINDRE CARRES
!                  (LEAST SQUARES POLYNOMIAL FITTING)
!
!              REFERENCE: BASIC SCIENTIFIC SUBROUTINES, VOL. II
!                         F.R. RUCKDESCHEL, BYTE/MCGRAWW-HILL, 1981
!                         ET D'APRES J-P MOREAU, PARIS
!
!                                                   P.KOECHLIN 02-2004
!
! ------------------------------------------------------------------
!
!   ENTREE
!       ORDRE : ORDRE DE LA REGRESSION (ORDRE>0)
!               = DEGRE MAX DU POLYNOME QUE L'ON CHERCHE
!       E1    : SI E1/=0 : PARAMETRE SERVANT DE CRITERE POUR PASSER A
!               L'ORDRE SUPERIEUR
!               ON COMPARE E1 A LA DIMINUTION DE LA DEVIATION STANDARD
!       NPT   : NOMBRE DE POINTS (NPT>1)
!       XX    : ABSCISSE DES POINTS (XX DISTINCTS)
!       YY    : ORDONNEES DES POINTS
!
!   SORTIE
!       ORDOK : ORDRE FINALEMENT OBTENU >=0
!                  (  AVANT VERIFICATION DU DEGRE:
!                        ORDOK <= ORDRE  SI E1/=0
!                        ORDOK = ORDRE   SI E1=0     )
!                  ORDOK=-1 CORRESPOND A UNE ERREUR
!       POLY  : POLYNOME OBTENU
!               POLY=POLY(1) + POLY(2)*X + POLY(2)*X^2 +...+ POLY(N)*X^N
!       SIGMA : DEVIATION STANDARD
!
!
!      LA ROUTINE N'EST VALABLE QUE POUR ORDRE >=1
!      MAIS EN SORTIE, ON A ORDOK >=0
!
    real(kind=8) :: grand
    parameter(grand=1.0d20)
!
    integer(kind=8) :: ordre, npt
    real(kind=8) :: xx(npt), yy(npt)
    real(kind=8) :: poly(ordre+1), sigma
!
    real(kind=8) :: aa(ordre+1), bb(ordre+1), ff(ordre+1), poly2(ordre+1)
    real(kind=8) :: vv(npt), dd(npt), ee(npt)
    integer(kind=8) :: i, j, ordok, ordok2, ordloo
    real(kind=8) :: a1, a2, b1, b2, poly1, d1, e1, f1, f2, v1, v2
!
    do i = 1, ordre+1
        poly(i) = 0.d0
    end do
    ordok = -1
!
    if (ordre .lt. 1) then
        goto 999
    end if
!
!      ON A TOUJOURS DES POINTS DISTINCTS (LES XX SONT DISTINCTS)
!      ON A TOUJOURS NPT>1
!
    if (npt .eq. 2) then
!         UNIQUEMENT POUR + DE PRECISION: MAIS LA ROUTINE MARCHE AUSSI
!         POUR NPT=2
        ordok = 1
        poly(2) = (yy(2)-yy(1))/(xx(2)-xx(1))
        poly(1) = (yy(1)+yy(2)-poly(2)*(xx(1)+xx(2)))/2.d0
        sigma = 0.d0
        goto 999
    end if
!
    v1 = grand
!
! --- INITIALISATION ------------------------------
!
    do i = 1, ordre+1
        aa(i) = 0.d0
        bb(i) = 0.d0
        ff(i) = 0.d0
    end do
    do i = 1, npt
        vv(i) = 0.d0
        dd(i) = 0.d0
    end do
!
    d1 = sqrt(npt*1.0d0)
!
    do i = 1, npt
        ee(i) = 1.d0/d1
    end do
    f1 = d1
!
    a1 = 0.d0
    do i = 1, npt
        a1 = a1+xx(i)*ee(i)*ee(i)
    end do
!
    poly1 = 0.d0
    do i = 1, npt
        poly1 = poly1+yy(i)*ee(i)
    end do
!
    bb(1) = 1.d0/f1
    ff(1) = bb(1)*poly1
!
    do i = 1, npt
        vv(i) = vv(i)+poly1*ee(i)
    end do
!
! --- DEBUT BOUCLE ----------------------------------------
!
    do ordloo = 1, ordre
!
! SAVE LATEST RESULTS
        do i = 1, ordre+1
            poly2(i) = poly(i)
        end do
        ordok2 = ordok
        v2 = v1
        f2 = f1
        a2 = a1
!
        f1 = 0.d0
        do i = 1, npt
            b1 = ee(i)
            ee(i) = (xx(i)-a2)*b1-f2*dd(i)
            dd(i) = b1
            f1 = f1+ee(i)*ee(i)
        end do
!
        f1 = sqrt(f1)
        do i = 1, npt
            ee(i) = ee(i)/f1
        end do
        a1 = 0.d0
        do i = 1, npt
            a1 = a1+xx(i)*ee(i)*ee(i)
        end do
!
        poly1 = 0.d0
        do i = 1, npt
            poly1 = poly1+yy(i)*ee(i)
        end do
!
        do i = 0, ordloo
            j = ordloo-i+1
            b2 = bb(j)
            d1 = 0.d0
            if (j .gt. 1) d1 = bb(j-1)
            d1 = d1-a2*bb(j)-f2*aa(j)
            bb(j) = d1/f1
            aa(j) = b2
        end do
!
        do i = 1, npt
            vv(i) = vv(i)+poly1*ee(i)
        end do
!
        do i = 1, ordre+1
            ff(i) = ff(i)+poly1*bb(i)
        end do
!
        do i = 1, ordre+1
            poly(i) = ff(i)
        end do
!
        ordok = ordloo
!
        sigma = 0.d0
        do i = 1, npt
            sigma = sigma+(vv(i)-yy(i))*(vv(i)-yy(i))
        end do
!
!        NOTE THE DIVISION IS BY THE NUMBER OF DEGREES OF FREEDOM
        if (npt .gt. ordloo+1) then
!          SIGMA = SQRT(SIGMA / DFLOAT(NPT - ORDLOO - 1))
            sigma = sqrt(sigma/(npt-ordloo-1))
        else
            goto 999
        end if
!
        if (e1 .gt. 0.d0) then
!          TEST FOR MINIMAL IMPROVEMENT OR IF ERROR IS LARGER, QUIT
            if ((abs(v1-sigma) .lt. (e1*sigma)) .or. (e1*sigma .gt. e1*v1)) then
!           ABORTED SEQUENCE, RECOVER LAST VALUES
                ordok = ordok2
                sigma = v2
                do i = 1, ordre+1
                    poly(i) = poly2(i)
                end do
                goto 999
            end if
        end if
!
        v1 = sigma
!
    end do
999 continue
end subroutine
