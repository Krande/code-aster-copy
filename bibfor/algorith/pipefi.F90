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
subroutine pipefi(npg, lgpg, mate, geom, vim,&
                  ddepl, deplm, ddepl0, ddepl1, dtau,&
                  typmod, compor, copilo)
!
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/nmfisa.h"
#include "asterfort/pipeba.h"
#include "asterfort/pipetu.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer :: mate, npg, lgpg
    real(kind=8) :: geom(2, 4), vim(lgpg, npg), ddepl(2, 4), deplm(2, 4)
    real(kind=8) :: ddepl0(2, 4), ddepl1(2, 4), dtau
    character(len=8) :: typmod(2)
    character(len=16) :: compor
    real(kind=8) :: copilo(5, npg)
!
!-----------------------------------------------------------------------
!  PILOTAGE PRED_ELAS POUR LES ELEMENTS DE JOINT 2D
!
! IN  : GEOM, MATE, VIM, DDEPL, DEPLM, DDEPL0, DDELP1, DTAU, NPG
! OUT : COPILO
!-----------------------------------------------------------------------
!
    aster_logical :: axi
    integer :: i, j, kpg
    real(kind=8) :: up(8), ud(8), sup(2), sud(2), b(2, 8), poids
!-----------------------------------------------------------------------
!
! INITIALISATION
!
    axi = typmod(1) .eq. 'AXIS'
!
! DEPLACEMENT U(ETA) = UP + ETA * UD
!
    call dcopy(8, deplm, 1, up, 1)
    call daxpy(8, 1.d0, ddepl, 1, up,&
               1)
    call daxpy(8, 1.d0, ddepl0, 1, up,&
               1)
    call dcopy(8, ddepl1, 1, ud, 1)
! BOUCLE SUR LES POINTS DE GAUSS :
!
    do kpg = 1, npg
!
!      SAUT AU POINT DE GAUSS : SU(ETA) = SUP + ETA * SUD
        call nmfisa(axi, geom, kpg, poids, b)
        do i = 1, 2
            sup(i) = 0.d0
            sud(i) = 0.d0
            do j = 1, 8
                sup(i) = sup(i) + b(i,j)*up(j)
                sud(i) = sud(i) + b(i,j)*ud(j)
            end do
        end do
!
!      INITIALISATION DES COEFFICIENTS DE PILOTAGE
        call r8inir(4, 0.d0, copilo(1, kpg), 1)
        copilo(5,kpg) = r8vide()
!
!      APPEL DU PILOTAGE PRED_ELAS SPECIFIQUE A LA LOI DE COMPORTEMENT
        if ((compor .eq. 'CZM_LIN_REG') .or. (compor .eq. 'CZM_EXP_REG')) then
            call pipeba(2, mate, sup, sud, vim(1, kpg),&
                        dtau, copilo(1, kpg))
        else if (compor .eq. 'CZM_TURON') then
            call pipetu(2, mate, sup, sud, vim(1, kpg),&
                        dtau, copilo(1, kpg))
        else
            call utmess('F', 'MECANONLINE_59')
        endif
!
    end do
!
end subroutine
