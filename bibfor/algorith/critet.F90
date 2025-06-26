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
subroutine critet(epsp, epsd, eta, lambda, deuxmu, &
                  fpd, seuil, crit, critp)
!
!
    implicit none
#include "asterfort/bptobg.h"
#include "asterfort/diagp3.h"
#include "asterfort/r8inir.h"
#include "blas/ddot.h"
    real(kind=8) :: epsp(6), epsd(6)
    real(kind=8) :: lambda, deuxmu
    real(kind=8) :: eta, fpd, seuil
    real(kind=8) :: crit, critp
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS - ENDO_ISOT_BETON)
!
! CALCUL DU CRITERE DE F(ETA) ET DE SA DERIVEE
!
! ----------------------------------------------------------------------
!
!
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! IN  ETA    : PARAMETRE DE PILOTAGE
! IN  LAMBDA : COEFFICIENT DE LAME
! IN  DEUXMU : COEFFICIENT DE LAME
! IN  FPD    : PARTIE DEPENDANTE DE D DANS LA FORCE THERMO
!             VAUT (1+GAMMA)/(1+GAMMA*D)**2
! IN  SEUIL  : SEUIL DU CRITERE
! OUT CRIT   : VALEUR DU CRITERE POUR ETA DONNEE EN ENTREE
! OUT CRITP  : VALEUR DE LA DERIVEE DU CRITERE POUR ETA DONNEE EN ENTREE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: k, i
    real(kind=8) :: tr(6), vecp(3, 3)
    real(kind=8) :: epm(3), tre, rac2
    real(kind=8) :: treps, sigel(3), ppeps(6), dfde(6)
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    rac2 = sqrt(2.d0)
! -- ON DIAGONALISE LE TENSEUR DE DEFORMATION
    tr(1) = epsp(1)+eta*epsd(1)
    tr(2) = (epsp(4)+eta*epsd(4))/rac2
    tr(3) = (epsp(5)+eta*epsd(5))/rac2
    tr(4) = epsp(2)+eta*epsd(2)
    tr(5) = (epsp(6)+eta*epsd(6))/rac2
    tr(6) = epsp(3)+eta*epsd(3)
!
    call diagp3(tr, vecp, epm)
!
! -- CALCUL DU CRITERE
!
    treps = epm(1)+epm(2)+epm(3)
    if (treps .gt. 0.d0) then
        do k = 1, 3
            sigel(k) = lambda*treps
        end do
    else
        do k = 1, 3
            sigel(k) = 0.d0
        end do
    end if
    do k = 1, 3
        if (epm(k) .gt. 0.d0) then
            sigel(k) = sigel(k)+deuxmu*epm(k)
        end if
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    crit = fpd*0.5d0*ddot(b_n, epm, b_incx, sigel, b_incy)-seuil
!
    do i = 1, 3
        if (epm(i) .lt. 0.d0) then
            tr(i) = 0.d0
        else
            tr(i) = epm(i)
        end if
        tr(i+3) = 0.d0
    end do
!
! -- CALCUL DE LA DERIVEE DU CRITERE
!
    call bptobg(tr, ppeps, vecp)
    call r8inir(6, 0.d0, dfde, 1)
    tre = epm(1)+epm(2)+epm(3)
!
    if (tre .gt. 0.d0) then
        do i = 1, 3
            dfde(i) = fpd*lambda*tre
        end do
    end if
    do i = 1, 3
        dfde(i) = dfde(i)+deuxmu*fpd*ppeps(i)
    end do
    do i = 4, 6
        dfde(i) = deuxmu*fpd*ppeps(i)*rac2
    end do
!
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    critp = ddot(b_n, dfde, b_incx, epsd, b_incy)
!
!
end subroutine
