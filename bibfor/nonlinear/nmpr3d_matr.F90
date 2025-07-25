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

subroutine nmpr3d_matr(nno, npg, poidsg, vff, dff, &
                       geom, p, matc)
!
    implicit none
!
#include "asterfort/r8inir.h"
#include "asterfort/subaco.h"
#include "asterfort/subacv.h"
#include "asterfort/sumetr.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: nno
    integer(kind=8), intent(in) :: npg
    real(kind=8), intent(in) :: poidsg(npg)
    real(kind=8), intent(in) :: vff(nno, npg)
    real(kind=8), intent(in) :: dff(2, nno, npg)
    real(kind=8), intent(in) :: geom(3, nno)
    real(kind=8), intent(in) :: p(npg)
    real(kind=8), intent(out) :: matc(3, nno, 3, nno)
!
! --------------------------------------------------------------------------------------------------
!
! Loads computation
!
! Pressure for faces of 3D elements - Tangent matrix
!
! --------------------------------------------------------------------------------------------------
!
!
! In  nno       : number of nodes
! In  nng       : number of Gauss points
! In  poidsg    : weight of Gauss points
! In  vff       : shape functions at Gauss points
! In  dff       : derivative of shape functions at Gauss point point
! In  geom      : coordinates of nodes
! In  p         : pressure at Gauss points
! Out matc      : tangent matrix (following force)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpg, n, i, m, j
    real(kind=8) :: cova(3, 3), metr(2, 2), jac, cnva(3, 2)
    real(kind=8) :: t1, t2, t3, t, acv(2, 2)
!
! --------------------------------------------------------------------------------------------------
!
! - Initializations
!
    call r8inir(nno*nno*9, 0.d0, matc, 1)
!
! - Loop on Gauss points
!
    do kpg = 1, npg
!
! ----- Covariant basis
!
        call subaco(nno, dff(1, 1, kpg), geom, cova)
!
! ----- Metric tensor
!
        call sumetr(cova, metr, jac)
!
! ----- Contra-variant basis
!
        call subacv(cova, metr, jac, cnva, acv)
!
! ----- Tangent matrix
!
        do m = 1, nno
            do j = 1, 3
                do n = 1, nno
                    do i = 1, 3
                        t1 = (dff(1, m, kpg)*cnva(j, 1)+ &
                              dff(2, m, kpg)*cnva(j, 2))*vff(n, kpg)*cova(i, 3)
                        t2 = dff(1, m, kpg)*cova(j, 3)*vff(n, kpg)*cnva(i, 1)
                        t3 = dff(2, m, kpg)*cova(j, 3)*vff(n, kpg)*cnva(i, 2)
                        t = poidsg(kpg)*p(kpg)*jac*(t1-t2-t3)
                        matc(i, n, j, m) = matc(i, n, j, m)+t
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
