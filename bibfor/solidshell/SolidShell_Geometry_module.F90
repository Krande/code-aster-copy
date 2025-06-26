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
! Module for geometry of cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Geometry_module
! ==================================================================================================
    use SolidShell_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public   :: compJacoMatr, shapeDeriHexa
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! shapeDeriHexa
!
! Compute the decomposed derivatives of shape functions at current point for HEXA
!
! In  XI               : parametric coordinates of current point
! Out dN_dXsi          : derivatives of shapes function by XI
! Out dN_dEta          : derivatives of shapes function by ETA
! Out dN_dZeta         : derivatives of shapes function by DZETA
!
! --------------------------------------------------------------------------------------------------
    subroutine shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in)  :: XI(3)
        real(kind=8), intent(out) :: dN_dXsi(SSH_NBNODEG_MAX)
        real(kind=8), intent(out) :: dN_dEta(SSH_NBNODEG_MAX)
        real(kind=8), intent(out) :: dN_dZeta(SSH_NBNODEG_MAX)
! - Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom
!   ------------------------------------------------------------------------------------------------
!
        dN_dXsi = 0.d0
        dN_dEta = 0.d0
        dN_dZeta = 0.d0
        do iNodeGeom = 1, nbNodeGeom
            dN_dXsi(iNodeGeom) = hexaVectG1(iNodeGeom)+ &
                                 XI(2)*hexaVectH1(iNodeGeom)+ &
                                 XI(3)*hexaVectH3(iNodeGeom)+ &
                                 XI(2)*XI(3)*hexaVectH4(iNodeGeom)
            dN_dEta(iNodeGeom) = hexaVectG2(iNodeGeom)+ &
                                 XI(1)*hexaVectH1(iNodeGeom)+ &
                                 XI(3)*hexaVectH2(iNodeGeom)+ &
                                 XI(1)*XI(3)*hexaVectH4(iNodeGeom)
            dN_dZeta(iNodeGeom) = hexaVectG3(iNodeGeom)+ &
                                  XI(2)*hexaVectH2(iNodeGeom)+ &
                                  XI(1)*hexaVectH3(iNodeGeom)+ &
                                  XI(2)*XI(1)*hexaVectH4(iNodeGeom)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compJacoMatr
!
! Compute jacobian matrix of element at current point
!
! In  elemProp         : general properties of element
! In  geomCurr         : current configuration of element
! In  XI               : parametric coordinates of current point
! Out Jac              : jacobian matrix at current point
! Out JacInv           : inverse of jacobian matrix at current point
! Out det              : determinant of Jacobian  at current point
!
! --------------------------------------------------------------------------------------------------
    subroutine compJacoMatr(elemProp, geomCurr, XI, &
                            Jac_, JacInv_, det_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        real(kind=8), intent(in)  :: geomCurr(SSH_NBDOFG_MAX), XI(3)
        real(kind=8), optional, intent(out) :: Jac_(3, 3), JacInv_(3, 3), det_
! - Local
        integer(kind=8) :: iNodeGeom, nbNodeGeom
       real(kind=8) :: dN_dXsi(SSH_NBNODEG_MAX), dN_dEta(SSH_NBNODEG_MAX), dN_dZeta(SSH_NBNODEG_MAX)
        real(kind=8) :: dX_dXsi, dX_dEta, dX_dZeta
        real(kind=8) :: dZ_dXsi, dZ_dEta, dZ_dZeta
        real(kind=8) :: dY_dXsi, dY_dEta, dY_dZeta
        real(kind=8) :: XCurr(SSH_NBNODEG_MAX), YCurr(SSH_NBNODEG_MAX), ZCurr(SSH_NBNODEG_MAX)
        real(kind=8) :: JacInv(3, 3), Jac(3, 3), det
!   ------------------------------------------------------------------------------------------------
!
        Jac = 0.d0
        JacInv = 0.d0
        det = 0.d0

! - Get parameters
        nbNodeGeom = elemProp%nbNodeGeom
        do iNodeGeom = 1, nbNodeGeom
            XCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+1)
            YCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+2)
            ZCurr(iNodeGeom) = geomCurr(3*(iNodeGeom-1)+3)
        end do

! - Derivatives of shape functions
        if (elemProp%cellType == SSH_CELL_HEXA) then
            call shapeDeriHexa(XI, dN_dXsi, dN_dEta, dN_dZeta)
        else
            ASSERT(ASTER_FALSE)
        end if

! - Compute jacobian
        dX_dXsi = sum(XCurr*dN_dXsi)
        dX_dEta = sum(XCurr*dN_dEta)
        dX_dZeta = sum(XCurr*dN_dZeta)
        dY_dXsi = sum(YCurr*dN_dXsi)
        dY_dEta = sum(YCurr*dN_dEta)
        dY_dZeta = sum(YCurr*dN_dZeta)
        dZ_dXsi = sum(ZCurr*dN_dXsi)
        dZ_dEta = sum(ZCurr*dN_dEta)
        dZ_dZeta = sum(ZCurr*dN_dZeta)
        Jac(1, 1) = dX_dXsi
        Jac(1, 2) = dY_dXsi
        Jac(1, 3) = dZ_dXsi
        Jac(2, 1) = dX_dEta
        Jac(2, 2) = dY_dEta
        Jac(2, 3) = dZ_dEta
        Jac(3, 1) = dX_dZeta
        Jac(3, 2) = dY_dZeta
        Jac(3, 3) = dZ_dZeta
        Jac = transpose(Jac)

! - Inverse jacobian matrix
        call matinv('S', 3, Jac, JacInv, det)

! - Set outputs
        if (present(det_)) then
            det_ = det
        end if
        if (present(JacInv_)) then
            JacInv_ = JacInv
        end if
        if (present(Jac_)) then
            Jac_ = Jac
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Geometry_module
