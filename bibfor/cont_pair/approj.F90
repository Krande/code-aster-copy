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

subroutine approj(mesh, newgeo, sdcont_defi, node_mast_indx, l_pair_dire, &
                  pair_vect, iter_maxi, epsi_maxi, tole_proj_ext, poin_coor, &
                  elem_mast_mini, proj_stat_mini, ksi1_mini, ksi2_mini, tau1_mini, &
                  tau2_mini, dist_mini, vect_pm_mini, err_appa)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterfort/apchoi.h"
#include "asterfort/apcoma.h"
#include "asterfort/apdist.h"
#include "asterfort/aptypm.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfinvm.h"
#include "asterfort/cfnben.h"
#include "asterfort/cfnumm.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/mmproj.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    character(len=24), intent(in) :: sdcont_defi
    integer(kind=8), intent(in) :: node_mast_indx
    aster_logical, intent(in) :: l_pair_dire
    real(kind=8), intent(in) :: pair_vect(3)
    integer(kind=8), intent(in) :: iter_maxi
    real(kind=8), intent(in) :: epsi_maxi
    real(kind=8), intent(in) :: tole_proj_ext
    real(kind=8), intent(in) :: poin_coor(3)
    real(kind=8), intent(out) :: tau1_mini(3)
    real(kind=8), intent(out) :: tau2_mini(3)
    real(kind=8), intent(out) :: vect_pm_mini(3)
    real(kind=8), intent(out) :: ksi1_mini
    real(kind=8), intent(out) :: ksi2_mini
    real(kind=8), intent(out) :: dist_mini
    integer(kind=8), intent(out) :: proj_stat_mini
    integer(kind=8), intent(out) :: elem_mast_mini
    integer(kind=8), intent(inout) :: err_appa
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Projection of contact point on master element
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  node_mast_indx   : nearest master node index in contact datastructure
! In  l_pair_dire      : .true. if using DIRE_APPA (given direction for pairing)
! In  pair_vect        : value of DIRE_APPA (given direction for pairing)
! In  iter_maxi        : maximum number of Newton iterations
! In  epsi_maxi        : maximum tolerance for Newton algorithm
! In  tole_proj_ext    : tolerance for outside element projection
! In  poin_coor        : coordinate of contact point to project
! Out dist_mini        : distance between point and projection of the point for selection
! Out vect_pm_mini     : vector between point and projection of the point for selection
! Out ksi1_mini        : first parametric coordinate of the projection of the point for selection
! Out ksi2_mini        : second parametric coordinate of the projection of the point for selection
! Out proj_stat_mini   : status of projection
!                            0 - Inside element
!                            1 - Inside element + tole_proj_out
!                            2 - Outside element
! Out elem_indx_mini   : index of element in contact datastructure for selection
! Out tau1_mini        : first tangent vector for local basis for selection
! Out tau2_mini        : second tangent vector for local basis for selection
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: elem_mast_type, elem_mast_name
    integer(kind=8) :: elem_mast_ndim, niverr, elem_mast_nbnode, node_nbelem, ino
    real(kind=8) :: elem_mast_coor(27), vect_pm(3)
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: ksi1, ksi2, dist
    integer(kind=8) :: proj_stat, i_elem_mast, elem_mast_indx, elem_mast_nume
    integer(kind=8) :: jdeciv, mesh_ndim
    aster_logical :: l_poi1, l_cont_disc

!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)

    dist_mini = r8gaem()
    ksi1_mini = r8gaem()
    ksi2_mini = r8gaem()
    tau1_mini(1:3) = 0.d0
    tau2_mini(1:3) = 0.d0
    vect_pm_mini(1:3) = 0.d0
    elem_mast_mini = 0
    proj_stat_mini = -1
!
! - No node excluded
!
    ASSERT(node_mast_indx .ne. 0)
!
! - Dimension of the mesh
!
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=mesh_ndim)
    ASSERT((mesh_ndim == 2) .or. (mesh_ndim == 3))
!
    l_cont_disc = cfdisl(sdcont_defi, 'FORMUL_DISCRETE')
!
! - Number of elements attached to master node
!
    call cfnben(sdcont_defi, node_mast_indx, 'CONINV', node_nbelem, jdeciv)
!
! - Loop on master elements
!
    do i_elem_mast = 1, node_nbelem
!
! ----- Get master elements attached to current master node
!
        call cfinvm(sdcont_defi, jdeciv, i_elem_mast, elem_mast_indx)
!
! ----- Index of master element
!
        call cfnumm(sdcont_defi, elem_mast_indx, elem_mast_nume)
!
! ----- Number of nodes
!
        call cfnben(sdcont_defi, elem_mast_indx, 'CONNEX', elem_mast_nbnode)
!
! ----- Parameters of master element
!
        call aptypm(mesh, elem_mast_nume, elem_mast_ndim, elem_mast_nbnode, elem_mast_type, &
                    elem_mast_name)
!
! ----- Coordinates of master elements
!
        call apcoma(mesh, newgeo, elem_mast_nume, elem_mast_nbnode, elem_mast_coor)
!
        if ((elem_mast_ndim .lt. mesh_ndim) .and. .not. l_cont_disc) then
! ----- The nodes have to be in the plane xOy (problem for shell element)
            do ino = 1, elem_mast_nbnode
                if (abs(elem_mast_coor(3*(ino-1)+3)) > r8prem()) then
                    call utmess('F', 'APPARIEMENT_12', ni=2, vali=[mesh_ndim, elem_mast_ndim])
                end if
            end do
        end if
!
! ----- No POI1 master element
!
        l_poi1 = elem_mast_type .eq. 'PO1'
        if (l_poi1) then
            call utmess('F', 'APPARIEMENT_36', sk=elem_mast_name)
        end if
!
! ----- Projection of node on master element
!
        call mmproj(elem_mast_type, elem_mast_nbnode, elem_mast_ndim, elem_mast_coor, poin_coor, &
                    iter_maxi, epsi_maxi, tole_proj_ext, l_pair_dire, pair_vect, &
                    ksi1, ksi2, tau1, tau2, proj_stat, &
                    niverr)
!
! ----- Management of error
!
        if (niverr .eq. 1) then
            err_appa = 1
            if (niv .ge. 2) then
                call utmess('I', 'APPARIEMENT_13', sk=elem_mast_name, nr=3, valr=poin_coor)
            end if
        end if
!
! ----- Compute distance
!
        call apdist(elem_mast_type, elem_mast_coor, elem_mast_nbnode, ksi1, ksi2, &
                    poin_coor, dist, vect_pm)
!
! ----- Select nearest element
!
        call apchoi(dist, dist_mini, elem_mast_indx, elem_mast_mini, tau1, &
                    tau1_mini, tau2, tau2_mini, ksi1, ksi1_mini, &
                    ksi2, ksi2_mini, proj_stat, proj_stat_mini, vect_pm, &
                    vect_pm_mini)
    end do
!
end subroutine
