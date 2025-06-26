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

subroutine vtcreb_wrap(field_nodez, base, type_scalz, nume_equaz)
!
    implicit none
!
#include "asterfort/vtcreb.h"
#include "asterfort/dismoi.h"
!
!
    character(len=*), intent(in) :: field_nodez
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: type_scalz
    character(len=*), intent(in) :: nume_equaz
!
! --------------------------------------------------------------------------------------------------
!
! Field utility
!
! Create NODE field
!
! --------------------------------------------------------------------------------------------------
!
! In  field_node    : name of field
! In  base          : JEVEUX base to create field
! In  type_scal     : type of GRANDEUR (real or complex)
! With numbering:
!   In  nume_ddl    : name of numbering
!
!   Out nb_equa_out : number of equations
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mesh
    integer(kind=8) :: idx_gd, nb_equa
!
    call dismoi('NUM_GD_SI', nume_equaz, 'NUME_EQUA', repi=idx_gd)
    call dismoi('NB_EQUA', nume_equaz, 'NUME_EQUA', repi=nb_equa)
    call dismoi('NOM_MAILLA', nume_equaz, 'NUME_EQUA', repk=mesh)

    call vtcreb(field_nodez, base, type_scalz, meshz=mesh, nume_equaz=nume_equaz, &
                idx_gdz=idx_gd, nb_equa_inz=nb_equa)
!
end subroutine
