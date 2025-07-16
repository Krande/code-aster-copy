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

subroutine crprno(nume_equaz, base, meshz, gran_namez, nb_equa)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jenonu.h"
#include "asterfort/nume_equa_crsd.h"
#include "asterfort/jexnom.h"
!
!
    character(len=*), intent(in) :: nume_equaz
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: gran_namez
    character(len=*), intent(in) :: meshz
    integer(kind=8), intent(in) :: nb_equa
!
! --------------------------------------------------------------------------------------------------
!
! Create NUME_EQUA only on mesh
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_equa   : name of NUME_EQUA
! In  base        : JEVEUX base to create NUME_EQUA
! In  nb_equa     : number of equations
! In  gran_name   : name of GRANDEUR
! In  mesh        : name of mesh
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_ligr_mesh
    character(len=24) :: lili
!
! --------------------------------------------------------------------------------------------------
!
    call nume_equa_crsd(nume_equaz, base, nb_equa, meshz=meshz, &
                        gran_namez=gran_namez)

    lili = nume_equaz(1:19)//'.LILI'
    call jenonu(jexnom(lili, '&MAILLA'), i_ligr_mesh)
    ASSERT(i_ligr_mesh .eq. 1)

end subroutine
