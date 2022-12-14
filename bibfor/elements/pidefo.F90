! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine pidefo(ndim, npg, kpg, compor, fm,&
                  epsm, epsp, epsd, copilo)
!
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
    integer :: ndim, kpg, npg
    character(len=16) :: compor(*)
    real(kind=8) :: epsm(6), epsp(6), epsd(6)
    real(kind=8) :: fm(3, 3)
    real(kind=8) :: copilo(5, npg)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS/DEFORMATION)
!
! PILOTAGE PAR DEFORMATION
!
! ----------------------------------------------------------------------
!
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NPG    : NOMBRE DE POINTS DE GAUSS
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  COMPOR : COMPORTEMENT
! IN  FM     : GRADIENT DE LA TRANSFORMATION AU TEMPS MOINS
! IN  EPSM   : DEFORMATIONS AU TEMPS MOINS
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! OUT COPILO : COEFFICIENTS A0 ET A1 POUR CHAQUE POINT DE GAUSS
!
!
!
!
    aster_logical :: grand
    integer :: ndimsi
    integer :: indi(6), indj(6), prac(6)
    real(kind=8) :: ff
    real(kind=8) :: rac2
    real(kind=8) :: em(6), epsmno
    integer :: ij, kl, i, j, k, l
!
    data indi /1,2,3,2,3,3/
    data indj /1,2,3,1,1,2/
    data prac /0,0,0,1,1,1/
!
! ----------------------------------------------------------------------
!
!
!
! --- INITIALISATIONS
!
    rac2 = sqrt(2.d0)
    grand = compor(3) .ne. 'PETIT'
    ndimsi = 2*ndim
!
! --- TRANSPORT DU TENSEUR DES DEFORMATIONS E := F E FT
!
    if (grand) then
        call dcopy(ndimsi, epsm, 1, em, 1)
        call r8inir(ndimsi, 0.d0, epsm, 1)
!
        do ij = 1, ndimsi
            do kl = 1, ndimsi
                i = indi(ij)
                j = indj(ij)
                k = indi(kl)
                l = indj(kl)
                ff = (fm(i,k)*fm(j,l) + fm(i,l)*fm(j,k)) / 2
                ff = ff * rac2**prac(ij) * rac2**prac(kl)
                epsm(ij) = epsm(ij) + ff*em(kl)
            end do
        end do
    endif
!
! --- INCREMENT DE DEFORMATION PROJETE
!
    epsmno = dnrm2(ndimsi,epsm ,1)
    copilo(1,kpg) = ddot(ndimsi, epsm,1, epsp,1)/epsmno
    copilo(2,kpg) = ddot(ndimsi, epsm,1, epsd,1)/epsmno
!
!
end subroutine
