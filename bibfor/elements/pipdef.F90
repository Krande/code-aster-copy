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
subroutine pipdef(ndim, nno, kpg, ipoids, ivf, &
                  idfde, geom, typmod, compor, deplm, &
                  ddepl, depl0, depl1, dfdi, fm, &
                  epsm, epsp, epsd)
!
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/nmgeom.h"
#include "blas/daxpy.h"
    integer(kind=8) :: ndim, nno, kpg
    integer(kind=8) :: ipoids, ivf, idfde
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*)
    real(kind=8) :: geom(ndim, *), deplm(*)
    real(kind=8) :: ddepl(*), depl0(*), depl1(*)
    real(kind=8) :: dfdi(*)
    real(kind=8) :: epsm(6), epsp(6), epsd(6)
    real(kind=8) :: fm(3, 3)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS/DEFORMATION)
!
! CALCUL DES DEFORMATIONS
!
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  IPOIDS : POIDS DES POINTS DE GAUSS
! IN  IVF    : VALEUR DES FONCTIONS DE FORME
! IN  IDFDE  : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM   : COORDONEES DES NOEUDS
! IN  TYPMOD : TYPE DE MODELISATION
! IN  COMPOR : COMPORTEMENT
! IN  DEPLM  : DEPLACEMENT EN T-
! IN  DDEPL  : INCREMENT DE DEPLACEMENT A L'ITERATION NEWTON COURANTE
! IN  DEPL0  : CORRECTION DE DEPLACEMENT POUR FORCES FIXES
! IN  DEPL1  : CORRECTION DE DEPLACEMENT POUR FORCES PILOTEES
! OUT DFDI   : DERIVEE DES FONCTIONS DE FORME
! OUT FM     : GRADIENT DE LA TRANSFORMATION AU TEMPS MOINS
! OUT EPSM   : DEFORMATIONS AU TEMPS MOINS
! OUT EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! OUT EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
!
!
!
    aster_logical :: axi, grand
    integer(kind=8) :: ndimsi
    real(kind=8) :: r, deps(6)
    real(kind=8) :: t9bid(3, 3)
    real(kind=8) :: poids
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
!
!
! --- INITIALISATIONS
!
    axi = typmod(1) .eq. 'AXIS'
    grand = compor(3) .ne. 'PETIT'
    ndimsi = 2*ndim
!
    if (typmod(2) .eq. 'DEPLA') then
!
! ----- CALCUL DE EPSM (LINEAIRE) OU EM (GREEN)  = EPS(UM)
        call nmgeom(ndim, nno, axi, grand, geom, &
                    kpg, ipoids, ivf, idfde, deplm, &
                    .true._1, poids, dfdi, fm, epsm, &
                    r)
!
! ----- REACTUALISATION DE LA GEOMETRIE SI GRANDES DEFS
        if (grand) then
            b_n = to_blas_int(ndim*nno)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, deplm, b_incx, geom, &
                       b_incy)
        end if
!
! ----- CALCUL DE DEPS = EPS(DU)
        call nmgeom(ndim, nno, axi, .false._1, geom, &
                    kpg, ipoids, ivf, idfde, ddepl, &
                    .true._1, poids, dfdi, t9bid, deps, &
                    r)
!
! ----- CALCUL DE EPSP (= DEPS + EPS(DU0) )
        call nmgeom(ndim, nno, axi, .false._1, geom, &
                    kpg, ipoids, ivf, idfde, depl0, &
                    .true._1, poids, dfdi, t9bid, epsp, &
                    r)
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, epsp, &
                   b_incy)
!
! ----- CALCUL DE EPSD (DEPS = EPSP + ETA EPSD)
        call nmgeom(ndim, nno, axi, .false._1, geom, &
                    kpg, ipoids, ivf, idfde, depl1, &
                    .true._1, poids, dfdi, t9bid, epsd, &
                    r)
!
    else
        ASSERT(.false.)
    end if
!
!
!
!
end subroutine
