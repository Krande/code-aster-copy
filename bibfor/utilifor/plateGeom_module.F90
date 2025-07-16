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
! ==================================================================================================
!
! Module for geometry of plates
!
! ==================================================================================================
!
module plateGeom_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: checkPlaneity
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterfort/plateGeom_module.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! checkPlaneity
!
! Check planeity of quadrangular cell
!
! In  xyzg             : coordinates of vertices
! In  errorTole        : tolerance for check
! Out errorCode        : error return code
! Out distAbso         : distance to plane (absolute)
! Out distRela         : distance to plane (relative)
!
! --------------------------------------------------------------------------------------------------
    subroutine checkPlaneity(xyzg, errorTole, errorCode, distAbso, distRela)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: xyzg(3, 4)
        integer(kind=8), intent(out) :: errorCode
        real(kind=8), intent(in) :: errorTole
        real(kind=8), intent(out) :: distAbso, distRela
! ----- Local
        real(kind=8) :: x12, y12, z12, x13, y13, z13, x14, y14, z14
        real(kind=8) :: ux, uy, uz, pscal, normu, norm4, dist
!   ------------------------------------------------------------------------------------------------
!
        errorCode = BASE_NO_ERROR
        distAbso = 0.d0
        distRela = 0.d0

! ----- First vector
        x12 = xyzg(1, 2)-xyzg(1, 1)
        y12 = xyzg(2, 2)-xyzg(2, 1)
        z12 = xyzg(3, 2)-xyzg(3, 1)

! ----- Second vector
        x13 = xyzg(1, 3)-xyzg(1, 1)
        y13 = xyzg(2, 3)-xyzg(2, 1)
        z13 = xyzg(3, 3)-xyzg(3, 1)

! ----- Third vector
        x14 = xyzg(1, 4)-xyzg(1, 1)
        y14 = xyzg(2, 4)-xyzg(2, 1)
        z14 = xyzg(3, 4)-xyzg(3, 1)

! ----- Define plane on node 1-2-3
        ux = (y12*z13)-(y13*z12)
        uy = (z12*x13)-(z13*x12)
        uz = (x12*y13)-(x13*y12)

        pscal = (ux*x14)+(uy*y14)+(uz*z14)

! ----- Compute normal to plane
        normu = sqrt((ux*ux)+(uy*uy)+(uz*uz))
        if (normu .lt. r8miem()) then
! --------- Something wrong: degenerated normal
            errorCode = BASE_CELL_DEGE
        else
            norm4 = sqrt((x14*x14)+(y14*y14)+(z14*z14))
            dist = pscal/normu
            pscal = dist/norm4
            if (abs(pscal) .gt. errorTole) then
                errorCode = BASE_QUAD_NOPLANE
                distAbso = abs(dist)
                distRela = pscal*100
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module plateGeom_module
