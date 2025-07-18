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
subroutine diagp3(tens, vecp, valp)
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/zerop3.h"
    real(kind=8) :: tens(6), valp(3), vecp(3, 3)
!
! ----------------------------------------------------------------------
!  DIAGONALISATION MATRICE 3x3 SYMETRIQUE PAR UNE METHODE DIRECTE
!    IN    TENS   : TENSEUR SOUS LA FORME
!                     (XX XY XZ YY YZ ZZ)
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, nrac, ind
    aster_logical :: invvp, tnull
    real(kind=8) :: trace, x(6), y(6), det(4), rtemp
    real(kind=8) :: a, b, c, theta
    real(kind=8) :: f, g
!
!
!
    call r8inir(9, 0.d0, vecp, 1)
!
!      write (6,*) '-----------------------------------------------'
!      write (6,*) TENS(1),TENS(4),TENS(6),TENS(2),TENS(3),TENS(5)
! -- PASSAGE AU DEVIATEUR
    trace = (tens(1)+tens(4)+tens(6))/3.d0
!      write (6,*) TRACE,R8PREM()
    tens(1) = tens(1)-trace
    tens(4) = tens(4)-trace
    tens(6) = tens(6)-trace
!      write (6,*) TENS(1),TENS(4),TENS(6),TENS(2),TENS(3),TENS(5)
!
! -- CALCUL DES COEFFICIENTS DU POLYNOME P3
!
    det(1) = tens(4)*tens(6)-tens(5)**2
    det(2) = tens(1)*tens(6)-tens(3)**2
    det(3) = tens(1)*tens(4)-tens(2)**2
!
!
    det(4) = tens(1)*det(1)-tens(2)*(tens(2)*tens(6)-tens(5)*tens(3))&
     &       +tens(3)*(tens(2)*tens(5)-tens(4)*tens(3))
!
    call zerop3(0.d0, det(1)+det(2)+det(3), -det(4), valp, nrac)
!
! -- CAS DES RACINES DOUBLES QUI PASSENT MAL NUMERIQUEMENT
    if (nrac .eq. 1) then
        rtemp = (det(1)+det(2)+det(3))
        if (rtemp .le. 0.d0) then
            valp(2) = sqrt(-rtemp/3.d0)
            valp(3) = -sqrt(-rtemp/3.d0)
            if (abs(valp(2)**3+rtemp*valp(2)-det(4)) .lt. &
                abs(valp(3)**3+rtemp*valp(3)-det(4))) then
                valp(3) = valp(2)
            else
                valp(2) = valp(3)
            end if
        else
            tnull = .true.
            do i = 1, 6
                if (abs(tens(i)) .gt. (r8prem()*100*abs(trace))) then
                    tnull = .false.
                end if
            end do
            if (tnull) then
                valp(1) = 0.d0
                valp(2) = 0.d0
                valp(3) = 0.d0
            else
                call utmess('F', 'ALGORITH2_79')
            end if
        end if
        if (valp(2) .gt. valp(1)) then
            rtemp = valp(1)
            valp(1) = valp(2)
            valp(3) = rtemp
        end if
    end if
!
    rtemp = valp(1)
    valp(1) = valp(3)
    valp(3) = rtemp
!
!      write (6,*) 'VALEURS PROPRES :',VALP(1),VALP(2),VALP(3)
!
! -- ON RECHERCHE LA VALEUR PROPRE LA PLUS SEPAREE
!    ON LA MET DANS VALP(1)
    if ((valp(2)-valp(1)) .gt. (valp(3)-valp(1))/2.d0) then
        invvp = .false.
    else
        invvp = .true.
        rtemp = valp(1)
        valp(1) = valp(3)
        valp(3) = rtemp
    end if
!
! -- VECP DE LA VAL PROPRE LA PLUS SEPAREE
!      write (6,*) 'VAL PROPRE APRES TRI : ',VALP(1),VALP(2),VALP(3)
! -- ON MULTIPLIE LES 3 VECT DE BASE PAR (A-LAMBDA_3.ID)(A-LAMBDA_2.ID)
!      ON PRENDRA CELUI DONT LA NORME EST LA PLUS GRANDE
    do ind = 1, 3
        if (ind .eq. 1) then
            x(1) = tens(1)-valp(2)
            x(2) = tens(2)
            x(3) = tens(3)
        else if (ind .eq. 2) then
            x(1) = tens(2)
            x(2) = tens(4)-valp(2)
            x(3) = tens(5)
        else
            x(1) = tens(3)
            x(2) = tens(5)
            x(3) = tens(6)-valp(2)
        end if
        vecp(1, ind) = (tens(1)-valp(3))*x(1)+tens(2)*x(2)+tens(3)*x(3)
        vecp(2, ind) = tens(2)*x(1)+(tens(4)-valp(3))*x(2)+tens(5)*x(3)
        vecp(3, ind) = tens(3)*x(1)+tens(5)*x(2)+(tens(6)-valp(3))*x(3)
        y(ind) = (vecp(1, ind))**2+(vecp(2, ind))**2+(vecp(3, ind))**2
    end do
    rtemp = y(1)
    ind = 1
    if (y(2) .gt. y(1)) then
        if (y(3) .gt. y(2)) then
            ind = 3
        else
            ind = 2
        end if
    else if (y(3) .gt. y(1)) then
        ind = 3
    end if
    a = sqrt(y(ind))
! -- CAS DE 3 VALEURS PROPRES EGALES
    if (a .lt. r8miem()) then
        call r8inir(9, 0.d0, vecp, 1)
        vecp(1, 1) = 1.d0
        vecp(2, 2) = 1.d0
        vecp(3, 3) = 1.d0
        goto 999
    end if
    do i = 1, 3
        vecp(i, 1) = vecp(i, ind)/a
    end do
!
! -- AUTRES VECTEURS PROPRES : ON PASSE DANS LE SOUS-ESPACE
!    ORTHOGONAL AU PREMIER VECTEUR PROPRE
!
! -- ON CHERCHE DEUX VECTEURS ORTHOGONAUX AU PREMIER VECT PROPRE
!    ON COMMENCE PAR PRENDRE CELUI ORTHOGONAL A UN VECT DE BASE
!    DONT LA NORME EST LA PLUS GRANDE (POUR LE CAS OU LE PREMIER
!    VECTEUR PROPRE SERAIT UN VECTEUR DE BASE)
    y(1) = vecp(3, 1)**2+vecp(2, 1)**2
    y(2) = vecp(3, 1)**2+vecp(1, 1)**2
    y(3) = vecp(1, 1)**2+vecp(2, 1)**2
    rtemp = y(1)
    ind = 1
    if (y(2) .gt. y(1)) then
        if (y(3) .gt. y(2)) then
            ind = 3
        else
            ind = 2
        end if
    else if (y(3) .gt. y(1)) then
        ind = 3
    end if
    a = sqrt(y(ind))
    if (ind .eq. 1) then
        vecp(1, 2) = 0.d0
        vecp(2, 2) = -vecp(3, 1)/a
        vecp(3, 2) = vecp(2, 1)/a
    else if (ind .eq. 2) then
        vecp(1, 2) = vecp(3, 1)/a
        vecp(2, 2) = 0.d0
        vecp(3, 2) = -vecp(1, 1)/a
    else
        vecp(1, 2) = -vecp(2, 1)/a
        vecp(2, 2) = vecp(1, 1)/a
        vecp(3, 2) = 0.d0
    end if
    vecp(1, 3) = vecp(2, 1)*vecp(3, 2)-&
     &                 vecp(3, 1)*vecp(2, 2)
    vecp(2, 3) = vecp(3, 1)*vecp(1, 2)-&
     &                 vecp(1, 1)*vecp(3, 2)
    vecp(3, 3) = vecp(1, 1)*vecp(2, 2)-&
     &                      vecp(2, 1)*vecp(1, 2)
!
! -- ON PASSE DANS LE SOUS-ESPACE ORTHOGONAL
!
    a = vecp(1, 2)*(vecp(1, 2)*tens(1)+2.d0*vecp(2, 2)*tens(2))&
     & +vecp(2, 2)*(vecp(2, 2)*tens(4)+2.d0*vecp(3, 2)*tens(5))&
     & +vecp(3, 2)*(vecp(3, 2)*tens(6)+2.d0*vecp(1, 2)*tens(3))
    b = vecp(1, 3)*(vecp(1, 3)*tens(1)+2.d0*vecp(2, 3)*tens(2))&
     & +vecp(2, 3)*(vecp(2, 3)*tens(4)+2.d0*vecp(3, 3)*tens(5))&
     & +vecp(3, 3)*(vecp(3, 3)*tens(6)+2.d0*vecp(1, 3)*tens(3))
    c = vecp(1, 2)*(vecp(1, 3)*tens(1)+vecp(2, 3)*tens(2)+vecp(3, 3)*tens(3)&
     &)+vecp(2, 2)*(vecp(1, 3)*tens(2)+vecp(2, 3)*tens(4)+vecp(3, 3)*tens(5)&
     &)+vecp(3, 2)*(vecp(1, 3)*tens(3)+vecp(2, 3)*tens(5)+vecp(3, 3)*tens(6)&
     &)
! -- ON CHERCHE L'ANGLE DONT EST TOURNE LE REPERE PROPRE
!
    f = 2.d0*c
    g = a-b
    if (abs(f) .lt. r8miem()) then
        theta = 0.d0
    else
        theta = atan2(f, g)
        theta = theta/2.d0
    end if
!
! -- EST-CE THETA OU THETA+PI/2 ?
!
    rtemp = (a-b)*cos(2.d0*theta)+2.d0*c*sin(2.d0*theta)
    if (rtemp*(valp(2)-valp(3)) .lt. 0.d0) theta = theta+r8pi()/2.d0
!        write (6,*) 'THETA =',THETA
    a = cos(theta)
    b = sin(theta)
    y(1) = vecp(1, 2)*a+vecp(1, 3)*b
    y(2) = vecp(2, 2)*a+vecp(2, 3)*b
    y(3) = vecp(3, 2)*a+vecp(3, 3)*b
    vecp(1, 2) = y(1)
    vecp(2, 2) = y(2)
    vecp(3, 2) = y(3)
!
    vecp(1, 3) = vecp(2, 1)*vecp(3, 2)-&
     &                vecp(3, 1)*vecp(2, 2)
    vecp(2, 3) = vecp(3, 1)*vecp(1, 2)-&
     &                vecp(1, 1)*vecp(3, 2)
    vecp(3, 3) = vecp(1, 1)*vecp(2, 2)-&
     &                      vecp(2, 1)*vecp(1, 2)
!
    if (invvp) then
        do i = 1, 3
            rtemp = vecp(i, 1)
            vecp(i, 1) = vecp(i, 3)
            vecp(i, 3) = -rtemp
        end do
        rtemp = valp(1)
        valp(1) = valp(3)
        valp(3) = rtemp
    end if
!
!        Y(1)=TENS(1)
!        Y(2)=TENS(4)
!        Y(3)=TENS(6)
!        Y(4)=TENS(2)
!        Y(5)=TENS(3)
!        Y(6)=TENS(5)
!        CALL BGTOBP(Y,X,VECP)
!
!        IF ((ABS(X(1)-VALP(1)).GT.TOL).OR.
!     &      (ABS(X(2)-VALP(2)).GT.TOL).OR.
!     &      (ABS(X(3)-VALP(3)).GT.TOL)) THEN
!          write (6,*) 'X(1) = ',X(1),' ; VALP(1) = ',VALP(1)
!          write (6,*) 'X(2) = ',X(2),' ; VALP(2) = ',VALP(2)
!          write (6,*) 'X(3) = ',X(3),' ; VALP(3) = ',VALP(3)
!        ENDIF
!
!
999 continue
!
    do i = 1, 3
        valp(i) = valp(i)+trace
    end do
!
!      IF (MOD(INT(TPS(2)),100000).EQ.0) THEN
!        write (6,*) 'NB APP = ',TPS(2),' ; TOT = ',TPS(3),
!     &              ' ; MOY = ',TPS(4)
!      ENDIF
end subroutine
