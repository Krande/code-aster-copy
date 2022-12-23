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

subroutine te0455(nomopt, nomte)
!
    use Behaviour_module, only: behaviourOption
    use HHO_type
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Meca_module
    use HHO_init_module, only: hhoInfoInitCell
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmtstm.h"
#include "asterfort/writeVector.h"
#include "asterfort/writeMatrix.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - rigidity and residual
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer :: cbs, fbs, total_dofs
    integer :: jmatt, icompo, npg
    integer :: codret, jcret, jgrad, jstab, icarcr
    aster_logical :: l_largestrains, lMatr, lVect, lSigm, lVari, matsym
    character(len=4), parameter :: fami = 'RIGI'
    character(len=8) :: typmod(2)
    character(len=16) :: defo_comp, type_comp
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC)
    real(kind=8), dimension(MSIZE_CELL_MAT, MSIZE_TDOFS_VEC)   :: gradfull
    real(kind=8), dimension(MSIZE_TDOFS_VEC, MSIZE_TDOFS_VEC)  :: lhs, stab
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(cbs <= MSIZE_CELL_VEC)
    ASSERT(fbs <= MSIZE_FACE_VEC)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)
!
    if (nomopt /= "RIGI_MECA_TANG" .and. &
        nomopt /= "RIGI_MECA_ELAS" .and. &
        nomopt /= "RIGI_MECA" .and. &
        nomopt /= "FULL_MECA" .and. &
        nomopt /= "RAPH_MECA") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature for the rigidity
!
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Type of finite element
!
    select case (hhoCell%ndim)
    case (3)
        typmod(1) = '3D'
    case (2)
        if (lteatt('AXIS', 'OUI')) then
            ASSERT(ASTER_FALSE)
            typmod(1) = 'AXIS'
        else if (lteatt('C_PLAN', 'OUI')) then
            ASSERT(ASTER_FALSE)
            typmod(1) = 'C_PLAN'
        else if (lteatt('D_PLAN', 'OUI')) then
            typmod(1) = 'D_PLAN'
        else
            ASSERT(ASTER_FALSE)
        end if
    case default
        ASSERT(ASTER_FALSE)
    end select
    typmod(2) = 'HHO'
!
    if (nomopt .ne. "RIGI_MECA") then
!
! --- Get input fields
!
        call jevech('PCOMPOR', 'L', icompo)
!
! --- Properties of behaviour
!
        defo_comp = zk16(icompo-1+DEFO)
        type_comp = zk16(icompo-1+INCRELAS)
!
! --- Large strains ?
!
        l_largestrains = isLargeStrain(defo_comp)
!
        call behaviourOption(nomopt, zk16(icompo), lMatr, lVect, lVari, lSigm, codret)
    else
        l_largestrains = ASTER_FALSE
        lMatr = ASTER_TRUE
        lSigm = ASTER_FALSE
        lVect = ASTER_FALSE
        codret = 0
        icompo = 1
    end if
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
        call jevech('PCHHOGT', 'L', jgrad)
        call jevech('PCHHOST', 'L', jstab)
!
        call hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, zr(jgrad), zr(jstab), &
                                  gradfull, stab)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, gradfull, stab)
    end if
!
! --- Compute local contribution
!
    call hhoLocalContribMeca(hhoCell, hhoData, hhoQuadCellRigi, gradfull, stab, &
                                & fami, typmod, zk16(icompo), nomopt, &
                                & l_largestrains, lhs, rhs, codret)
!
! --- Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = codret
    end if
!
! --- Save rhs
!
    if (lVect) then
        call hhoRenumMecaVec(hhoCell, hhoData, rhs)
        call writeVector('PVECTUR', total_dofs, rhs)
    end if
!
! --- Save of lhs
!
    if (lMatr) then
        if (nomopt .ne. "RIGI_MECA") then
            call jevech('PCARCRI', 'L', icarcr)
            call nmtstm(zr(icarcr), jmatt, matsym)
        else
            matsym = ASTER_TRUE
        end if
        call hhoRenumMecaMat(hhoCell, hhoData, lhs)
!
        if (matsym) then
            call writeMatrix('PMATUUR', total_dofs, total_dofs, ASTER_TRUE, lhs)
        else
            call writeMatrix('PMATUNS', total_dofs, total_dofs, ASTER_FALSE, lhs)
        end if
    end if
!
end subroutine
