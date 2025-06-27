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
!
#include "asterf_types.h"
!
interface
    subroutine mmapma(mesh, ds_contact, model_ndim, i_zone,&
                      lexfro, typint, aliase, posmae, node_mast_nume,&
                      nnomae, elem_mast_indx, elem_mast_nume, ksipr1, ksipr2,&
                      tau1m, tau2m, iptm, iptc, norm,&
                      nommam)
        use NonLin_Datastructure_type
        character(len=8) :: mesh
        character(len=8) :: aliase
        type(NL_DS_Contact), intent(in) :: ds_contact
        real(kind=8) :: ksipr1, ksipr2
        integer(kind=8) :: model_ndim
        integer(kind=8) :: posmae, node_mast_nume
        integer(kind=8) :: elem_mast_indx, elem_mast_nume, nnomae
        integer(kind=8) :: i_zone, iptm, iptc
        integer(kind=8) :: typint
        real(kind=8) :: tau1m(3), tau2m(3), norm(3)
        character(len=8) :: nommam
        aster_logical :: lexfro
    end subroutine mmapma
end interface
