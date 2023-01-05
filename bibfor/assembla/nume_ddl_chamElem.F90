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

subroutine nume_ddl_chamElem(nume, ligrel, modeloc)
    !
    implicit none
    !
#include "asterfort/numero.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
    !
    !
    character(len=*), intent(in) :: nume
    character(len=*), intent(in) :: ligrel
    character(len=*), intent(in) :: modeloc
    !
    ! ----------------------------------------------------------------------------------------------
    !
    ! Factor
    !
    ! Numbering - Create Nmbreing for a CHAM_ELNO
    !
    ! ----------------------------------------------------------------------------------------------
    !
    ! In  nume_ddl       : name of nume_ddl object
    ! In  list_ligrel      : list of ligrel
    ! In  modloc           : name of mode local
    !
    ! ----------------------------------------------------------------------------------------------
    !
    character(len=14) :: nume_ddl
    character(len=24), pointer :: v_ligrel(:) => null()
    integer :: nb_grel
    !
    ! ----------------------------------------------------------------------------------------------
    !
    nume_ddl = nume
    !
    call jeveuo(ligrel, 'L', vk24=v_ligrel)
    call jelira(ligrel, 'LONUTI', nb_grel)
    !
    ! ----- CALCUL DE LA NUMEROTATION PROPREMENT DITE :
    !
    call numero(nume_ddl, 'GG', modelocz=modeloc, &
                nb_ligrel=nb_grel, list_ligrel=v_ligrel)
    !
    !
end subroutine
