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

subroutine splin1(x, y, d2y, n, ptx, &
                  dyptx, iret)
    implicit none
! DESCRIPTION : INTERPOLATION SPLINE CUBIQUE
! -----------
!               ETANT DONNEE LA TABULATION DE LA FONCTION Y(I) = F(X(I))
!               EN N POINTS DE DISCRETISATION X(I)
!               ETANT DONNEE LA TABULATION DE LA DERIVEE SECONDE DE LA
!               FONCTION INTERPOLEE D2Y(I), CALCULEE EN AMONT PAR LA
!               ROUTINE SPLINE
!               ETANT DONNE UN POINT PTX
!               CETTE ROUTINE CALCULE LA VALEUR DYPTX DE L'INTERPOLATION
!               SPLINE CUBIQUE DE LA DERIVEE PREMIERE DE LA FONCTION AU
!               POINT PTX
!
! IN     : X     : REAL*8 , VECTEUR DE DIMENSION N
!                  CONTIENT LES POINTS DE DISCRETISATION X(I)
! IN     : Y     : REAL*8 , VECTEUR DE DIMENSION N
!                  CONTIENT LES VALEURS DE LA FONCTION AUX POINTS X(I)
! IN     : D2Y   : REAL*8 , VECTEUR DE DIMENSION N
!                  CONTIENT LES VALEURS DE LA DERIVEE SECONDE DE LA
!                  FONCTION INTERPOLEE AUX POINTS X(I)
! IN     : N     : INTEGER , SCALAIRE
!                  NOMBRE DE POINTS DE DISCRETISATION
! IN     : PTX   : REAL*8 , SCALAIRE
!                  VALEUR DU POINT OU L'ON SOUHAITE CALCULER LA DERIVEE
!                  PREMIERE DE LA FONCTION INTERPOLEE
! OUT    : DYPTX : REAL*8 , SCALAIRE
!                  VALEUR DE LA DERIVEE PREMIERE DE LA FONCTION
!                  INTERPOLEE AU POINT PTX
! OUT    : IRET  : INTEGER , SCALAIRE , CODE RETOUR
!                  IRET = 0  OK
!                  IRET = 1  DEUX POINTS CONSECUTIFS DE LA
!                            DISCRETISATION X(I) SONT EGAUX
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
    real(kind=8) :: x(*), y(*), d2y(*), ptx, dyptx
    integer(kind=8) :: n, iret
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: k, kinf, ksup
    real(kind=8) :: a, b, h
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    iret = 0
!
    kinf = 1
    ksup = n
10  continue
    if (ksup-kinf .gt. 1) then
        k = (ksup+kinf)/2
        if (x(k) .gt. ptx) then
            ksup = k
        else
            kinf = k
        end if
        goto 10
    end if
!
    h = x(ksup)-x(kinf)
    if (h .eq. 0.0d0) then
        iret = 1
        goto 9999
    end if
    a = (x(ksup)-ptx)/h
    b = (ptx-x(kinf))/h
    dyptx = (y(ksup)-y(kinf))/h+((1.0d0-3.0d0*a*a)*d2y(kinf)+(3.0d0*b*b-1.0d0)*d2y(ksup) &
                                 )*h/6.0d0
!
9999 continue
!
! --- FIN DE SPLIN1.
end subroutine
