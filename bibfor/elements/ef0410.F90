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
subroutine ef0410(nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdrepe.h"
#include "asterfort/vdxefgeElno.h"
!
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgt = 10
    integer(kind=8) :: jvNbsp, jvEfge, jvGeom, lzi, nb2
    integer(kind=8) :: nbLayer
    real(kind=8) :: efgeElno(8, 9)
    real(kind=8) :: matevn(2, 2, npgt), matevg(2, 2, npgt)
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PEFFORR', 'E', jvEfge)

! - Get objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb2 = zi(lzi-1+2)

! - Get properties of shell
    call jevech('PNBSP_I', 'L', jvNbsp)
    nbLayer = zi(jvNbsp)

! - Compute
    call vdxefgeElno(nomte, zr(jvGeom), &
                     nbLayer, efgeElno)

! - Compute matrix to change base
    call vdrepe(nomte, matevn, matevg)

! - From local to global
    call vdefro(nb2, matevn, efgeElno, zr(jvEfge))
!
end subroutine
