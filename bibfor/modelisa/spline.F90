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
subroutine spline(x, y, n, dy1, dyn, &
                  d2y, iret)
    implicit none
! DESCRIPTION : INTERPOLATION SPLINE CUBIQUE
! -----------
!
!               ETANT DONNEE LA TABULATION DE LA FONCTION Y(I) = F(X(I))
!               EN N POINTS DE DISCRETISATION X(I) TELS QUE
!                               X(1) < X(2) < ... < X(N)
!               ETANT DONNEES LES VALEURS DE LA DERIVEE PREMIERE DY1 ET
!               DYN AUX PREMIER ET DERNIER POINTS X(1) ET X(N)
!               CETTE ROUTINE CALCULE LES VALEURS DE LA DERIVEE SECONDE
!               D2Y(I) DE LA FONCTION INTERPOLEE AUX N POINTS X(I)
!
!               SI DY1 ET/OU DYN DEPASSENT EN VALEUR ABSOLUE LE PLUS
!               GRAND NOMBRE ACCESSIBLE PAR LA MACHINE, LA SOLUTION
!               CALCULEE EST TELLE QUE LA DERIVEE SECONDE EST NULLE
!               AUX BORNES DE L'INTERVALLE
!
! IN     : X    : REAL*8 , VECTEUR DE DIMENSION N
!                 CONTIENT LES POINTS DE DISCRETISATION X(I)
! IN     : Y    : REAL*8 , VECTEUR DE DIMENSION N
!                 CONTIENT LES VALEURS DE LA FONCTION AUX POINTS X(I)
! IN     : N    : INTEGER , SCALAIRE
!                 NOMBRE DE POINTS DE DISCRETISATION
! IN     : DY1  : REAL*8 , SCALAIRE
!                 VALEUR DE LA DERIVEE PREMIERE DE LA FONCTION
!                 AU POINT X1
! IN     : DYN  : REAL*8 , SCALAIRE
!                 VALEUR DE LA DERIVEE PREMIERE DE LA FONCTION
!                 AU POINT XN
! OUT    : D2Y  : REAL*8 , VECTEUR DE DIMENSION N
!                 CONTIENT LES VALEURS DE LA DERIVEE SECONDE
!                 DE LA FONCTION INTERPOLEE AUX POINTS X(I)
! OUT    : IRET : INTEGER , SCALAIRE , CODE RETOUR
!                 IRET = 0  OK
!                 IRET = 1  VALEUR DE N INVALIDE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    real(kind=8) :: x(*), y(*), dy1, dyn, d2y(*)
    integer(kind=8) :: n, iret
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: i
    real(kind=8) :: bignum, p, qn, sig, un
    real(kind=8), pointer :: work(:) => null()
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
    AS_ALLOCATE(vr=work, size=n)
!
    iret = 0
    if (n .lt. 3) then
        iret = 1
        goto 999
    end if
!
    bignum = 0.99d0/r8miem()
!
    if (dble(abs(dy1)) .gt. bignum) then
        d2y(1) = 0.0d0
        work(1) = 0.0d0
    else
        d2y(1) = -0.5d0
        work(1) = 3.0d0/(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1))-dy1)
    end if
!
    do i = 2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*d2y(i-1)+2.0d0
        d2y(i) = (sig-1.0d0)/p
        work(i) = ( &
                  6.0d0*( &
                  (y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1) &
                                                                            )-sig*work(1+i-2 &
                                                                                       ) &
                  )/p
    end do
!
    if (dble(abs(dyn)) .gt. bignum) then
        qn = 0.0d0
        un = 0.0d0
    else
        qn = 0.5d0
        un = 3.0d0/(x(n)-x(n-1))*(dyn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
!
    d2y(n) = (un-qn*work(1+n-2))/(qn*d2y(n-1)+1.0d0)
    do i = n-1, 1, -1
        d2y(i) = d2y(i)*d2y(i+1)+work(i)
    end do
!
999 continue
    AS_DEALLOCATE(vr=work)
    call jedema()
!
! --- FIN DE SPLINE.
end subroutine
