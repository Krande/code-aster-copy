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
subroutine te0494(nomopt, nomte)
!
    use HHO_type
    use HHO_size_module, only: hhoTherDofs
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_basis_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/writeVector.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/HHO_basis_module.h"
#include "blas/dcopy.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO - Thermics
!  Option: AFFE_CHAR_CINE_R
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomte, nomopt
!
! -- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_basis_face) :: hhoBasisFace
    real(kind=8) :: basis(6*MAX_FACE_COEF+MAX_CELL_COEF)
    integer(kind=8) :: dec, iFace, size
    blas_int :: b_n
    blas_int, parameter :: b_incx = 1, b_incy = 1
!
    ASSERT(nomopt .eq. 'HHO_PRECALC_BS')
!
! --- Retrieve HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    dec = 1
    do iFace = 1, hhoCell%nbfaces
        call hhoBasisFace%initialize(hhoCell%faces(iFace))
        size = maxval(hhoBasisFace%coeff_shift)-1
        b_n = to_blas_int(size)
        call dcopy(b_n, hhoBasisFace%coeff_mono, b_incx, basis(dec), b_incy)
        dec = dec+size
    end do
!
    call hhoBasisCell%initialize(hhoCell)
    size = maxval(hhoBasisCell%coeff_shift)-1
    b_n = to_blas_int(size)
    call dcopy(b_n, hhoBasisCell%coeff_mono, b_incx, basis(dec), b_incy)
    dec = dec+size
!
! -- Save - the name is not PCHHOBS because reading this field in basis
!
    call writeVector('PCHHOBO', dec-1, basis)
!
end subroutine
