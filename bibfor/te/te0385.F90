! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine te0385(nomopt, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
!
    implicit none
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/ntfcma.h"
#include "asterfort/runge6.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: nomte, nomopt
! ----------------------------------------------------------------------
!
!    - FONCTION REALISEE:  OPTION : 'HYDR_ELGA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! THERMIQUE NON LINEAIRE
!
! ----------------------------------------------------------------------
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer :: itemps, icomp, imate, itempm, itempp
    integer ::  ifon(6), kp, i
    real(kind=8) :: deltat, err, tpgm, tpgp, funcEF(MAX_BS)
    real(kind=8), pointer :: hydrgm(:) => null(), hydrgp(:) => null()
!
    if (nomopt .ne. "HYDR_ELGA") then
        ASSERT(ASTER_FALSE)
    end if
!
    call jevech('PCOMPOR', 'L', icomp)

    if (zk16(icomp) (1:9) .eq. 'THER_HYDR') then
        call FECell%init()
        call FEBasis%initCell(FECell)
        call FEQuadCell%initCell(FECell, "MASS")

        call jevech('PMATERC', 'L', imate)
        call jevech('PTEMPMR', 'L', itempm)
        call jevech('PTEMPPR', 'L', itempp)
        call jevech('PHYDRMR', 'L', vr=hydrgm)
        call jevech('PHYDRPR', 'E', vr=hydrgp)
        call jevech('PTEMPSR', 'L', itemps)

        deltat = zr(itemps+1)

        call ntfcma(zk16(icomp), zi(imate), ASTER_FALSE, ifon)

        do kp = 1, FEQuadCell%nbQuadPoints
            tpgm = 0.d0
            tpgp = 0.d0
            hydrgp(kp) = 0.d0
            funcEF = FEBasis%func(FEQuadCell%points_param(1:3, kp))
            do i = 1, FEBasis%size
                tpgm = tpgm+zr(itempm+i-1)*funcEF(i)
                tpgp = tpgp+zr(itempp+i-1)*funcEF(i)
            end do

            call runge6(ifon(3), deltat, tpgp, tpgm, hydrgm(kp), &
                        hydrgp(kp), err)
        end do
    end if
end subroutine
