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
subroutine cjsci1(crit, mater, deps, sigd, i1f, &
                  tract, iret)
    implicit none
#include "asterf_types.h"
! CJS : CALCUL DE I1F  --> I1 A T+DT
!       RESOLUTION DE L'EQUATION SCALAIRE F(I1) = 0 DU COMPORTEMENT
!       ELASTIQUE NON LINEAIRE AVEC
!                                            I1F+QINIT
!       F(I1) = I1F - I1D - 3 KOE TR(DEPS) (----------)**N
!                                           3 PA
!
! ----------------------------------------------------------------------
! IN  CRIT  : CRITERES DE CONVERGENCE
! IN  MATER : COEFFICIENTS MATERIAU A T+DT
! IN  DEPS  : INCREMENT DE DEFORMATION
! IN  SIGD  : CONTRAINTE A T
! OUT I1    : TRACE DE SIG A T+DT
!     TRACT : VARIABLE LOGIQUE INDIQUANT LA TRACTION (I1F > QINIT)
! OUT IRET  : CODE RETOUR DE LORS DE LA RESOLUTION DE L'EQUATION
!             SCALAIRE
!                 IRET=0 => PAS DE PROBLEME
!                 IRET=1 => ECHEC
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ndt, ndi, imax, iret
    parameter(imax=60)
    real(kind=8) :: mater(14, 2), crit(*), deps(6), sigd(6), i1d, i1f
    real(kind=8) :: trdeps, coef, pa, n, multi
    real(kind=8) :: x0, x1, x2, oldx2, y0, y1, y2
    real(kind=8) :: zero, un, deux, trois, qinit
    aster_logical :: tract
    integer(kind=8) :: i
!
    common/tdim/ndt, ndi
!
!
    data zero/0.d0/
    data un/1.d0/
    data deux/2.d0/
    data trois/3.d0/
!
!-----------------------------------------------------------------------
!       METHODE DE LA SECANTE
!-----------------------------------------------------------------------
!
!
    qinit = mater(13, 2)
!--->   DETERMINATION DE TERME COEF = 3 KOE TR(DEPS)
!
    trdeps = zero
    do i = 1, ndi
        trdeps = trdeps+deps(i)
    end do
!
    coef = mater(1, 1)/(un-deux*mater(2, 1))*trdeps
    pa = mater(12, 2)
    n = mater(3, 2)
!
    i1d = zero
    do i = 1, ndi
        i1d = i1d+sigd(i)
    end do
    if ((i1d+qinit) .ge. 0.d0) then
        i1d = -qinit+1.d-12*pa
    end if
!
!
!--->  TRAITEMENT EXPLICITE DE L'EQUATION POUR LE NIVEAU CJS 1
! - CAS N.0: NIVEAU CJS 1
!   +++++++++++++++++++++
!
    tract = .false.
    if (n .eq. zero) then
        i1f = i1d+coef
        if (i1f .ge. (-qinit)) then
            tract = .true.
        end if
        goto 999
    end if
!
!--->  TRAITEMENT DE L'EQUATION EN FONCTION DE TRACE DE DEPS
! - CAS N.1: TRACE NULLE
!   ++++++++++++++++++++
    if (trdeps .eq. zero) then
        i1f = i1d
    end if
!
! - CAS N.2: TRACE NEGATIVE (CHARGEMENT)
!   ++++++++++++++++++++++++++++++++++++
!
    if (trdeps .lt. zero) then
!
!       DERTEMINATION DES BORNES DE L'INTERVALLE DE RECHERCHE
!       AVEC Y0>0 ET Y1<0
!
        x0 = i1d
        y0 = x0-i1d-coef*((x0+qinit)/trois/pa)**n
        multi = 2.d0
        x1 = x0+qinit
        do i = 1, imax
            x1 = multi*x1
            y1 = x1-i1d-coef*((x1+qinit)/trois/pa)**n
            if (y1 .lt. zero) goto 25
        end do
        iret = 1
        goto 999
25      continue
!
!
!
!       RECHERCHE DU ZERO DE LA FONCTION ENTRE (X0,Y0) ET (X1,Y1)
!
        oldx2 = zero
!
        do i = 1, int(abs(crit(1)))
            x2 = (x0*y1-x1*y0)/(y1-y0)
            y2 = x2-i1d-coef*((x2+qinit)/trois/pa)**n
!
            if (abs((x2-oldx2)/x2) .lt. crit(3) .or. y2 .eq. zero) goto 40
!
            oldx2 = x2
            if (y2 .gt. zero) then
                x0 = x2
                y0 = y2
            else
                x1 = x2
                y1 = y2
            end if
!
        end do
        iret = 1
        goto 999
40      continue
!
        i1f = x2
!
    end if
!
! - CAS N.3: TRACE POSITIVE (DECHARGEMENT)
!   ++++++++++++++++++++++++++++++++++++++
!
    if (trdeps .gt. zero) then
!
!       DERTEMINATION DES BORNES DE L'INTERVALLE DE RECHERCHE
!       AVEC Y0<0 ET Y1>0
!
        x0 = i1d
        y0 = x0-i1d-coef*((x0+qinit)/trois/pa)**n
        x1 = -qinit
        y1 = x1-i1d
!
!       RECHERCHE DU ZERO DE LA FONCTION ENTRE (X0,Y0) ET (X1,Y1)
!
        oldx2 = zero
!
        do i = 1, int(abs(crit(1)))
            x2 = (x0*y1-x1*y0)/(y1-y0)
            y2 = x2-i1d-coef*((x2+qinit)/trois/pa)**n
!
            if (abs((x2-oldx2)/x2) .lt. crit(3) .or. y2 .eq. zero) goto 70
!
            oldx2 = x2
            if (y2 .gt. zero) then
                x1 = x2
                y1 = y2
            else
                x0 = x2
                y0 = y2
            end if
!
        end do
        iret = 1
        goto 999
70      continue
!
        i1f = x2
!
    end if
!
999 continue
!
end subroutine
