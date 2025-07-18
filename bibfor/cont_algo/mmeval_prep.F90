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

subroutine mmeval_prep(mesh, time_curr, model_ndim, ds_contact, &
                       i_zone, &
                       ksipc1, ksipc2, ksipr1, ksipr2, &
                       tau1, tau2, &
                       elem_slav_indx, elem_slav_nbno, &
                       elem_slav_type, elem_slav_coor, &
                       elem_mast_nume, &
                       lagr_cont_node, &
                       norm, &
                       gap, gap_user, lagr_cont_poin, &
                       poin_slav_coor, poin_proj_coor)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/cfdist.h"
#include "asterfort/mmvalp.h"
#include "asterfort/mmvalp_scal.h"
#include "asterfort/mcopco.h"
#include "asterfort/mmnorm.h"
#include "asterfort/utmess.h"
#include "asterfort/mmnewj.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
! aslint: disable=W1504
!
    character(len=8), intent(in) :: mesh
    real(kind=8), intent(in) :: time_curr
    integer(kind=8), intent(in) :: model_ndim
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: i_zone
    real(kind=8), intent(in) :: ksipc1
    real(kind=8), intent(in) :: ksipc2
    real(kind=8), intent(in) :: ksipr1
    real(kind=8), intent(in) :: ksipr2
    real(kind=8), intent(in) :: tau1(3)
    real(kind=8), intent(in) :: tau2(3)
    integer(kind=8), intent(in) :: elem_slav_nbno
    integer(kind=8), intent(in) :: elem_slav_indx
    character(len=8), intent(in) :: elem_slav_type
    real(kind=8), intent(in) :: elem_slav_coor(27)
    integer(kind=8), intent(in) :: elem_mast_nume
    real(kind=8), intent(in) :: lagr_cont_node(9)
    real(kind=8), intent(out) :: norm(3)
    real(kind=8), intent(out) :: gap
    real(kind=8), intent(out) :: gap_user
    real(kind=8), intent(out) :: lagr_cont_poin
    real(kind=8), intent(out), optional :: poin_slav_coor(3), poin_proj_coor(3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Compute gap and contact pressure
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  time_curr        : current time
! In  model_ndim       : size of model
! In  ds_contact       : datastructure for contact management
! In  i_zone           : index of contact zone
! In  ksipc1           : first parametric coordinate of contact point in slave element
! In  ksipc2           : second parametric coordinate of contact point in slave element
! In  ksipr1           : first parametric coordinate of projection of contact point in master elem.
! In  ksipr2           : second parametric coordinate of projection of contact point in master elem.
! In  tau1             : first tangent vector for local basis
! In  tau1             : second tangent vector for local basis
! In  elem_slav_indx   : index of slave element (in contact datastructure)
! In  elem_slav_nbno   : number of nodes of slave element
! In  elem_slav_type   : type of slave element
! In  elem_slav_coor   : coordinates of slave element
! In  elem_mast_nume   : index of master element (in mesh datastructure)
! In  lagr_cont_node   : value of contact lagrangian at (slave) nodes
! Out norm             : normal vector for local basis
! Out gap              : contact gap
! Out gap_user         : contact gap defined by user
! Out lagr_cont_poin   : value of contact lagrangian at contact point
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: noor
    character(len=19) :: newgeo
    character(len=8) :: elem_mast_name
!
! --------------------------------------------------------------------------------------------------
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! - Coordinates of the contact point
!
    call mmvalp(model_ndim, elem_slav_type, elem_slav_nbno, 3, ksipc1, &
                ksipc2, elem_slav_coor, poin_slav_coor)
!
! - Coordinates of the projection of contact point
!
    call mcopco(mesh, newgeo, model_ndim, elem_mast_nume, ksipr1, &
                ksipr2, poin_proj_coor)
!
! - Local basis on master element
!
    call mmnorm(model_ndim, tau1, tau2, norm, noor)
    if (noor .le. r8prem()) then
        elem_mast_name = int_to_char8(elem_mast_nume)
        call utmess('F', 'CONTACT3_23', sk=elem_mast_name, nr=3, valr=poin_proj_coor)
    end if
!
! - Compute gap
!
    call mmnewj(model_ndim, poin_slav_coor, poin_proj_coor, norm, gap)
!
! - Get user gap
!
    call cfdist(ds_contact, i_zone, elem_slav_indx, poin_slav_coor, time_curr, &
                gap_user)
!
! - Interpolate contact pressure (Lagrange) at point
!
    call mmvalp_scal(model_ndim, elem_slav_type, elem_slav_nbno, ksipc1, ksipc2, &
                     lagr_cont_node, lagr_cont_poin)
!
end subroutine
