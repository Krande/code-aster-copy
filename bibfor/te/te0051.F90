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
subroutine te0051(option, nomte)
!
    use FE_topo_module
    use plateGeom_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/plateGeom_module.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT
!
! Options: VERI_PLAN
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    integer(kind=8) :: jvGeom, jvCodret, jvIndicr, jvPara
    integer(kind=8) :: nbNode, iNode
    integer(kind=8) :: iadzi, iazk24
    real(kind=8) :: nodeCoor(3, 4)
    integer(kind=8) :: errorCode
    real(kind=8) :: errorTole, distAbso, distRela
    aster_logical :: printError
    character(len=8) :: cellName
!
! --------------------------------------------------------------------------------------------------
!
    call FECell%init()
    nbNode = FECell%nbnodes

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PCHCKPR', 'L', jvPara)

! - Get output fields
    call jevech('PCODRET', 'E', jvCodret)
    call jevech('PINDICR', 'E', jvIndicr)

! - Check planeity
    if (nbNode .eq. 4) then
        do iNode = 1, nbNode
            nodeCoor(1, iNode) = zr(jvGeom+3*(iNode-1)-1+1)
            nodeCoor(2, iNode) = zr(jvGeom+3*(iNode-1)-1+2)
            nodeCoor(3, iNode) = zr(jvGeom+3*(iNode-1)-1+3)
        end do
        errorTole = zr(jvPara-1+1)
        printError = zr(jvPara-1+2) .gt. 0.d0
        call checkPlaneity(nodeCoor, errorTole, errorCode, distAbso, distRela)
        if (errorCode .eq. BASE_QUAD_NOPLANE) then
            zi(jvCodret) = 1
            zr(jvIndicr-1+1) = distAbso
            zr(jvIndicr-1+2) = distRela
            if (printError) then
                call tecael(iadzi, iazk24)
                cellName = zk24(iazk24-1+3) (1:8)
                call utmess("I", "PLATE1_80", sk=cellName, nr=2, valr=[distAbso, distRela])
            end if
        end if
    end if
!
end subroutine
