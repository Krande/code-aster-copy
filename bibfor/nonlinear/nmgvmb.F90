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
subroutine nmgvmb(ndim, nno1, nno2, npg, axi, &
                  geoi, vff1, vff2, idfde1, idfde2, &
                  iw, nddl, neps, b, w, &
                  ni2ldc)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
    aster_logical, intent(in) :: axi
    integer(kind=8), intent(in) :: ndim, nno1, nno2, npg, idfde1, idfde2, iw
    real(kind=8), intent(in) :: geoi(ndim, nno1)
    real(kind=8), intent(in) :: vff1(nno1, npg), vff2(nno2, npg)
    integer(kind=8), intent(out) :: nddl, neps
    real(kind=8), intent(out), allocatable :: b(:, :, :)
    real(kind=8), intent(out), allocatable :: w(:, :), ni2ldc(:, :)
! ----------------------------------------------------------------------
!  CALCUL DES ELEMENTS CINEMATIQUES POUR LA MODELISATION GRAD_VARI
! ----------------------------------------------------------------------
! IN  NDIM   DIMENSION DE L'ESPACE
! IN  NNO1   NOMBRE DE NOEUDS TOTAL (SUPPORT DES DEPLACEMENTS)
! IN  NNO2   NOMBRE DE NOEUDS SOMMETS (SUPPORT DE VI ET LAGRANGE)
! IN  NPG    NOMBRE DE POINTS DE GAUSS
! IN  AXI    .TRUE. SI MODELISATION AXIS
! IN  GEOM   COORDONNEES DES NOEUDS
! IN  VFF1   VALEUR DE LA FAMILLE DE FONCTIONS DE FORME NO 1
! IN  VFF2   VALEUR DE LA FAMILLE DE FONCTIONS DE FORME NO 2
! IN  IDFDE1 POINTEUR SUR LES DER. REFERENCE FAMILLE FCT FORME NO 1
! IN  IDFDE2 POINTEUR SUR LES DER. REFERENCE FAMILLE FCT FORME NO 2
! IN  IW     POINTEUR SUR LES POIDS DES PTS DE GAUSS DE REFERENCE
! OUT NDDL   NOMBRE DE DDL / ELEMENT
! OUT NEPS   NBR DE COMPOSANTE DE DEFORMATION (GENERALISEE)
! OUT B      MATRICE CINEMATIQUE EPS = B.U
! OUT W      POIDS DES POINTS DE GAUSS CONFIG INITIALE
! OUT NI2LDC CONVERSION CONTRAINTE STOCKEE -> CONTRAINTE LDC (AVEC RAC2)
! ----------------------------------------------------------------------
    real(kind=8), parameter :: rac2 = sqrt(2.d0), r2 = 0.5*rac2
    real(kind=8), parameter :: vrac2(6) = [1.d0, 1.d0, 1.d0, rac2, rac2, rac2]
! ----------------------------------------------------------------------
    integer(kind=8) :: g, n
    real(kind=8) :: r, unsurr, w0
    real(kind=8), allocatable :: dff1(:, :), dff2(:, :)
! ----------------------------------------------------------------------
!
#define to_aster_int(a) int(a, ASTER_INT_SIZE)
# define iu1(n,i) (n-1)*(ndim+2) + i
# define iu2(n,i) nno2*2 + (n-1)*ndim + i
# define ia(n) (n-1)*(ndim+2) + ndim + 1
# define il(n) (n-1)*(ndim+2) + ndim + 2
! ----------------------------------------------------------------------
!
    nddl = nno1*ndim+nno2*2
    neps = 3*ndim+2
    allocate (b(neps, npg, nddl), w(neps, npg), ni2ldc(neps, npg))
    allocate (dff1(nno1, ndim), dff2(nno2, ndim))
!
!
!
!
! - AFFECTATION DE LA MATRICE CINEMATIQUE B
!
    b = 0.d0
    do g = 1, npg
!
!       Derivee des fonctions de forme no 2 (r et w non utilise)
        call dfdmip(ndim, nno2, axi, geoi, g, &
                    iw, vff2(1, g), idfde2, r, w0, &
                    dff2)
!
!       Derivee des fonctions de forme no 1, rayon et poids
        call dfdmip(ndim, nno1, axi, geoi, g, &
                    iw, vff1(1, g), idfde1, r, w0, &
                    dff1)
        w(:, g) = w0
!
        if (ndim .eq. 2) then
            if (axi) then
                unsurr = 1/r
            else
                unsurr = 0
            end if
!
            do n = 1, nno2
                b(1, g, iu1(n, 1)) = dff1(n, 1)
                b(2, g, iu1(n, 2)) = dff1(n, 2)
                b(3, g, iu1(n, 1)) = vff1(n, g)*unsurr
                b(4, g, iu1(n, 1)) = r2*dff1(n, 2)
                b(4, g, iu1(n, 2)) = r2*dff1(n, 1)
                b(5, g, ia(n)) = vff2(n, g)
                b(6, g, il(n)) = vff2(n, g)
                b(7, g, ia(n)) = dff2(n, 1)
                b(8, g, ia(n)) = dff2(n, 2)
            end do
!
            do n = nno2+1, nno1
                b(1, g, iu2(n, 1)) = dff1(n, 1)
                b(2, g, iu2(n, 2)) = dff1(n, 2)
                b(3, g, iu2(n, 1)) = vff1(n, g)*unsurr
                b(4, g, iu2(n, 1)) = r2*dff1(n, 2)
                b(4, g, iu2(n, 2)) = r2*dff1(n, 1)
            end do
!
        else if (ndim .eq. 3) then
            do n = 1, nno2
                b(1, g, iu1(n, 1)) = dff1(n, 1)
                b(2, g, iu1(n, 2)) = dff1(n, 2)
                b(3, g, iu1(n, 3)) = dff1(n, 3)
                b(4, g, iu1(n, 1)) = r2*dff1(n, 2)
                b(4, g, iu1(n, 2)) = r2*dff1(n, 1)
                b(5, g, iu1(n, 1)) = r2*dff1(n, 3)
                b(5, g, iu1(n, 3)) = r2*dff1(n, 1)
                b(6, g, iu1(n, 2)) = r2*dff1(n, 3)
                b(6, g, iu1(n, 3)) = r2*dff1(n, 2)
                b(7, g, ia(n)) = vff2(n, g)
                b(8, g, il(n)) = vff2(n, g)
                b(9, g, ia(n)) = dff2(n, 1)
                b(10, g, ia(n)) = dff2(n, 2)
                b(11, g, ia(n)) = dff2(n, 3)
            end do
!
            do n = nno2+1, nno1
                b(1, g, iu2(n, 1)) = dff1(n, 1)
                b(2, g, iu2(n, 2)) = dff1(n, 2)
                b(3, g, iu2(n, 3)) = dff1(n, 3)
                b(4, g, iu2(n, 1)) = r2*dff1(n, 2)
                b(4, g, iu2(n, 2)) = r2*dff1(n, 1)
                b(5, g, iu2(n, 1)) = r2*dff1(n, 3)
                b(5, g, iu2(n, 3)) = r2*dff1(n, 1)
                b(6, g, iu2(n, 2)) = r2*dff1(n, 3)
                b(6, g, iu2(n, 3)) = r2*dff1(n, 2)
            end do
        end if
    end do
!
!
! - AFFECTATION DE LA FONCTION DE TRANSFERT SIGMA NICE --> SIGMA LDC
    do g = 1, npg
        ni2ldc(1:2*ndim, g) = vrac2(1:2*ndim)
    end do
    ni2ldc(2*ndim+1:neps, :) = 1.d0
!
!
    deallocate (dff1, dff2)
end subroutine
