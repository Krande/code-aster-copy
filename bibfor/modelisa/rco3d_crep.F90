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
! aslint: disable=W1501
!
subroutine rco3d_crep(ligrel, noma, chmlrac, cara_elem, epai)
    !
    implicit none
    !
#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"

    character(len=19), intent(in) :: ligrel, chmlrac
    character(len=8), intent(in) :: noma
    character(len=8), intent(in) :: cara_elem
    real(kind=8), intent(in) :: epai



    integer, parameter :: nceld1 = 4
    integer, parameter :: nceld2 = 4
    integer, parameter :: nceld3 = 4
    character(len=24) :: chmlrac_celv
    integer :: jv_chmlrac_celv, nb_grel, i_grel
    integer :: decal, i_liel, nb_liel, vale_indx
    character(len=24) :: chmlrac_celd
    integer, pointer :: v_chmlrac_celd(:) => null()
    integer, pointer :: v_ligrel_liel(:) => null()



    chmlrac_celd = chmlrac//'.CELD'
    chmlrac_celv = chmlrac//'.CELV'
    call jeveuo(chmlrac_celd, 'L', vi=v_chmlrac_celd)
    call jeveuo(chmlrac_celv, 'E', jv_chmlrac_celv)
    nb_grel = v_chmlrac_celd(2)

    write(*,*)  "nb grels    ", nb_grel

    do i_grel = 1, nb_grel
        decal = v_chmlrac_celd(nceld1+i_grel)
        nb_liel = v_chmlrac_celd(decal+1)
        call jeveuo(jexnum(ligrel//'.LIEL', i_grel), 'L', vi=v_ligrel_liel)

        do i_liel = 1, nb_liel
            vale_indx = jv_chmlrac_celv-1+v_chmlrac_celd(decal+nceld2+nceld3*(i_liel-1)+4)
            zr(vale_indx-1+1) = epai
            zr(vale_indx-1+2) = 0.0d0
            zr(vale_indx-1+3) = 0.0d0
            zr(vale_indx-1+4) = 0.0d0
            zr(vale_indx-1+5) = 0.0d0
            zr(vale_indx-1+6) = 0.0d0
        end do

        write(*,*) "############ nbliel  ", nb_liel
    end do



end subroutine
