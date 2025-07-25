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

subroutine aprema(sdappa, mesh, sdcont_defi, newgeo, err_appa)
!
    implicit none
!
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterf_types.h"
#include "asterfort/apcopt.h"
#include "asterfort/apinfi.h"
#include "asterfort/aporth.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisr.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfr.h"
#include "asterfort/approj.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sdcont_defi
    character(len=19), intent(in) :: newgeo
    integer(kind=8), intent(inout) :: err_appa

!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Find nearest element from current contact point
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  mesh             : name of mesh
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  newgeo           : name of field for geometry update from initial coordinates of nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_zone, i_poin, i, node_mast_indx
    mpi_int :: i_proc, nb_proc, mpicou
    integer(kind=8) :: nb_poin_mpi, nbr_poin_mpi, idx_start, idx_end
    integer(kind=8) :: nb_cont_zone, model_ndim, nt_poin
    integer(kind=8) :: nb_poin, proj_stat_mini, elem_mast_mini
    real(kind=8) :: poin_coor(3), tau1_mini(3), tau2_mini(3), dist_mini, ksi1_mini, ksi2_mini
    real(kind=8) :: pair_vect(3), tole_proj_ext, epsi_maxi, vect_pm_mini(3)
    integer(kind=8) :: iter_maxi
    aster_logical :: l_pair_dire, l_pair_masl, l_save, one_proc
    integer(kind=8) :: pair_type, pair_enti
    character(len=24) :: sdappa_dist, sdappa_appa
    integer(kind=8), pointer :: v_sdappa_appa(:) => null()
    real(kind=8), pointer :: v_sdappa_dist(:) => null()
    character(len=24) :: sdappa_tau1, sdappa_tau2, sdappa_proj
    real(kind=8), pointer :: v_sdappa_tau1(:) => null()
    real(kind=8), pointer :: v_sdappa_tau2(:) => null()
    real(kind=8), pointer :: v_sdappa_proj(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
    one_proc = .false.
    if (niv .ge. 2) then
        write (ifm, *) '<APPARIEMENT> RECH. MAILLE PLUS PROCHE'
    end if
!
! - Acces to pairing datastructure
!
    sdappa_appa = sdappa(1:19)//'.APPA'
    sdappa_dist = sdappa(1:19)//'.DIST'
    sdappa_tau1 = sdappa(1:19)//'.TAU1'
    sdappa_tau2 = sdappa(1:19)//'.TAU2'
    sdappa_proj = sdappa(1:19)//'.PROJ'
    call jeveuo(sdappa_appa, 'E', vi=v_sdappa_appa)
    call jeveuo(sdappa_dist, 'E', vr=v_sdappa_dist)
    call jeveuo(sdappa_tau1, 'E', vr=v_sdappa_tau1)
    call jeveuo(sdappa_tau2, 'E', vr=v_sdappa_tau2)
    call jeveuo(sdappa_proj, 'E', vr=v_sdappa_proj)
!
! - Get parameters
!
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    iter_maxi = cfdisi(sdcont_defi, 'PROJ_NEWT_ITER')
    epsi_maxi = cfdisr(sdcont_defi, 'PROJ_NEWT_RESI')
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nt_poin = cfdisi(sdcont_defi, 'NTPT')
!
! - Loop on contact zones
!
    i_poin = 0
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters on current zone
!
        nb_poin = mminfi(sdcont_defi, 'NBPT', i_zone)
        l_pair_masl = mminfi(sdcont_defi, 'APPARIEMENT', i_zone) .eq. 1
        l_pair_dire = mminfi(sdcont_defi, 'TYPE_APPA', i_zone) .eq. 1
        if (l_pair_dire) then
            pair_vect(1) = mminfr(sdcont_defi, 'TYPE_APPA_DIRX', i_zone)
            pair_vect(2) = mminfr(sdcont_defi, 'TYPE_APPA_DIRY', i_zone)
            pair_vect(3) = mminfr(sdcont_defi, 'TYPE_APPA_DIRZ', i_zone)
        end if
        tole_proj_ext = mminfr(sdcont_defi, 'TOLE_PROJ_EXT', i_zone)
!
! ----- Mpi informations
!
        call asmpi_comm('GET', mpicou)
        call asmpi_info(mpicou, rank=i_proc, size=nb_proc)
        if (one_proc) then
            nb_proc = 1
        end if
        nb_poin_mpi = int(nb_poin/nb_proc)
        nbr_poin_mpi = nb_poin-nb_poin_mpi*nb_proc
        idx_start = 1+(i_proc)*nb_poin_mpi
        idx_end = idx_start+nb_poin_mpi-1+(nbr_poin_mpi*int((i_proc+1)/nb_proc))

!
! ----- Loop on points
!
        do i = idx_start, idx_end
!
! --------- Point to paired ?
!
            call apinfi(sdappa, 'APPARI_TYPE', i_poin+i, pair_type)
            ASSERT(pair_type .ne. 0)
            if (l_pair_masl) then
!
! ------------- Coordinates of point
!
                call apcopt(sdappa, i_poin+i, poin_coor)
!
! ------------- Nearest master node
!
                call apinfi(sdappa, 'APPARI_ENTITE', i_poin+i, pair_enti)
                node_mast_indx = pair_enti
!
! ------------- Projection of contact point on master element
!
                call approj(mesh, newgeo, sdcont_defi, node_mast_indx, l_pair_dire, &
                            pair_vect, iter_maxi, epsi_maxi, tole_proj_ext, poin_coor, &
                            elem_mast_mini, proj_stat_mini, ksi1_mini, ksi2_mini, tau1_mini, &
                            tau2_mini, dist_mini, vect_pm_mini, err_appa)
!
! ------------- Orthogonalization of local basis
!
                call aporth(mesh, sdcont_defi, model_ndim, elem_mast_mini, poin_coor, &
                            tau1_mini, tau2_mini)
                if (pair_type .eq. 1) then
                    if (proj_stat_mini .eq. 2) then
                        pair_type = -3
                    else
                        pair_type = 2
                    end if
                end if
                l_save = .true.
            else
                l_save = .false.
            end if
!
! --------- Save
!
            if (l_save) then
                v_sdappa_appa(4*(i_poin+i-1)+1) = pair_type
                v_sdappa_appa(4*(i_poin+i-1)+2) = elem_mast_mini
                v_sdappa_appa(4*(i_poin+i-1)+3) = i_zone
                v_sdappa_dist(4*(i_poin+i-1)+1) = dist_mini
                v_sdappa_dist(4*(i_poin+i-1)+2) = vect_pm_mini(1)
                v_sdappa_dist(4*(i_poin+i-1)+3) = vect_pm_mini(2)
                v_sdappa_dist(4*(i_poin+i-1)+4) = vect_pm_mini(3)
                v_sdappa_proj(2*(i_poin+i-1)+1) = ksi1_mini
                v_sdappa_proj(2*(i_poin+i-1)+2) = ksi2_mini
                v_sdappa_tau1(3*(i_poin+i-1)+1) = tau1_mini(1)
                v_sdappa_tau1(3*(i_poin+i-1)+2) = tau1_mini(2)
                v_sdappa_tau1(3*(i_poin+i-1)+3) = tau1_mini(3)
                v_sdappa_tau2(3*(i_poin+i-1)+1) = tau2_mini(1)
                v_sdappa_tau2(3*(i_poin+i-1)+2) = tau2_mini(2)
                v_sdappa_tau2(3*(i_poin+i-1)+3) = tau2_mini(3)
            end if

        end do
!
! ----- Next zone
!
        i_poin = i_poin+nb_poin
    end do
!
end subroutine
