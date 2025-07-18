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

subroutine aplcno(mesh, newgeo, sdcont_defi, sdappa, err_appa)
!
    implicit none
!
#include "asterfort/jeveuo.h"
#include "asterfort/cfdisi.h"
#include "asterfort/mminfi.h"
#include "asterfort/aptgnn.h"
#include "asterfort/aptgem.h"
#include "asterfort/aptnol.h"
#include "asterfort/apsvnl.h"
#include "asterfort/sdmpic.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    character(len=24), intent(in) :: sdcont_defi
    character(len=19), intent(in) :: sdappa
    integer(kind=8), intent(inout) :: err_appa
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Compute norms at nodes (smoothing)
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  sdappa           : name of pairing datastructure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: model_ndim, nb_cont_zone, nt_node, nt_node_slav, nt_node_mast
    integer(kind=8) :: nb_node_mast, nb_node_slav, nb_elem_mast, nb_elem_slav
    integer(kind=8) :: i_zone
    integer(kind=8) :: jdecnm, jdecmm, jdecne, jdecme
    character(len=4) :: zone_type
    real(kind=8) :: norm_vect(3), epsi_maxi
    integer(kind=8) :: norm_type
    character(len=16), pointer :: valk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    epsi_maxi = 1.d-12
!
! - Get parameters
!
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nt_node_slav = cfdisi(sdcont_defi, 'NTNOE')
    nt_node_mast = cfdisi(sdcont_defi, 'NTNOM')
    nt_node = nt_node_slav+nt_node_mast
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Get parameters on current zone
!
        nb_node_mast = mminfi(sdcont_defi, 'NBNOM', i_zone)
        nb_node_slav = mminfi(sdcont_defi, 'NBNOE', i_zone)
        nb_elem_mast = mminfi(sdcont_defi, 'NBMAM', i_zone)
        nb_elem_slav = mminfi(sdcont_defi, 'NBMAE', i_zone)
        jdecnm = mminfi(sdcont_defi, 'JDECNM', i_zone)
        jdecmm = mminfi(sdcont_defi, 'JDECMM', i_zone)
        jdecne = mminfi(sdcont_defi, 'JDECNE', i_zone)
        jdecme = mminfi(sdcont_defi, 'JDECME', i_zone)
        norm_type = 0
!
! ----- MPI initialisation
!
        call jeveuo(sdappa(1:19)//'.MPIB', 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
        call jeveuo(sdappa(1:19)//'.MPIC', 'E', vk16=valk)
        valk(1) = 'MPI_INCOMPLET'
!
! ----- Compute tangents at each node for each master element
!
        zone_type = 'MAIT'
        call aptgem(sdappa, mesh, newgeo, sdcont_defi, model_ndim, &
                    i_zone, zone_type, epsi_maxi, jdecmm, &
                    nb_elem_mast, err_appa)
!
! ----- Compute tangents at each node for each slave element
!
        zone_type = 'ESCL'
        call aptgem(sdappa, mesh, newgeo, sdcont_defi, model_ndim, &
                    i_zone, zone_type, epsi_maxi, jdecme, &
                    nb_elem_slav, err_appa)
        call sdmpic('SD_APPA_TGEL', sdappa)
!
! ----- Compute tangents at each node by smoothing
!
        call aptgnn(sdappa, mesh, sdcont_defi, model_ndim, jdecnm, &
                    nb_node_mast, norm_type, norm_vect)
        call aptgnn(sdappa, mesh, sdcont_defi, model_ndim, jdecne, &
                    nb_node_slav, norm_type, norm_vect)
        call sdmpic('SD_APPA_TGNO', sdappa)
!
! ----- Compute normals at nodes
!
        call aptnol(sdappa, model_ndim, nt_node)
!
    end do
!
! - Smooth normals at nodes
!
    call apsvnl(sdcont_defi, sdappa, model_ndim, nt_node)
!
end subroutine
