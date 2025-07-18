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

subroutine apstos(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/mminfr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/aplcno.h"
#include "asterfort/aplcpg.h"
#include "asterfort/cfdisi.h"
#include "asterfort/jelira.h"
#include "asterfort/jerazo.h"
#include "asterfort/apstoc.h"
#include "asterfort/mmbouc.h"
#include "asterfort/infdbg.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"

! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Pairing - Segment to segment
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: newgeo, sdappa
    character(len=8) :: knuzo
    integer(kind=8) :: nb_elem_mast, nb_elem_slav
    character(len=24) :: sdappa_mast, sdappa_slav
    integer(kind=8) :: i_zone, nb_pair_zone, nb_cont_zone, nt_patch, zmeth
    aster_logical :: l_smooth
    real(kind=8) :: pair_tole
    character(len=24) :: pair_method
    integer(kind=8), pointer :: list_pair_zone(:) => null()
    integer(kind=8), pointer :: list_nbptit_zone(:) => null()
    real(kind=8), pointer :: list_ptitsl_zone(:) => null()
    real(kind=8), pointer :: list_ptitma_zone(:) => null()
    real(kind=8), pointer :: list_ptgama_zone(:) => null()
    character(len=24) :: sdappa_gapi
    integer(kind=8), pointer :: v_sdappa_mast(:) => null()
    integer(kind=8), pointer :: v_sdappa_slav(:) => null()
    integer(kind=8), pointer :: v_sdcont_methco(:) => null()
    character(len=24) :: sdcont_methco
    mpi_int :: i_proc, nb_proc, mpicou
    integer(kind=8) :: nb_elem_mpi, nbr_elem_mpi, idx_start, idx_end
    integer(kind=8), pointer :: v_appa_slav_mpi(:) => null()
    integer(kind=8) ::nb_el_slav_mpi, err_appa
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<Pairing> Segment-to-segment pairing'
    end if
!
! - Initializations
!
    nb_pair_zone = 0
!
! - Get parameters
!
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    l_smooth = cfdisl(ds_contact%sdcont_defi, 'LISSAGE')
    nt_patch = ds_contact%nt_patch
    zmeth = cfmmvd('ZMETH')
!
! - Updated geometry
!
    newgeo = ds_contact%sdcont_solv(1:14)//'.NEWG'
!
! - Access to pairing datastructures
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
    sdappa_gapi = sdappa(1:19)//'.GAPI'
    call jerazo(sdappa_gapi, nt_patch, 1)
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Get parameters
!
        pair_tole = mminfr(ds_contact%sdcont_defi, 'RESI_APPA', i_zone)
        sdcont_methco = ds_contact%sdcont_defi(1:16)//'.METHCO'
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 2) then
            pair_method = 'ROBUSTE'
        else if (v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 3) then
            pair_method = 'RAPIDE'
        end if
!
! ----- Generate name of objects
!
        ASSERT(i_zone .le. 100)
        call codent(i_zone-1, 'G', knuzo)
        sdappa_mast = sdappa(1:19)//'.MS'//knuzo(1:2)
        sdappa_slav = sdappa(1:19)//'.EC'//knuzo(1:2)
!
! ----- Get objects
!
        call jelira(sdappa_mast, 'LONMAX', nb_elem_mast)
        call jelira(sdappa_slav, 'LONMAX', nb_elem_slav)
        call jeveuo(sdappa_mast, 'L', vi=v_sdappa_mast)
        call jeveuo(sdappa_slav, 'L', vi=v_sdappa_slav)
!
! ----- Pairing
!
!
! ----- MPI initialisation
!
        call asmpi_comm('GET', mpicou)
        call asmpi_info(mpicou, rank=i_proc, size=nb_proc)
        nb_elem_mpi = int(nb_elem_slav/nb_proc)
        nbr_elem_mpi = nb_elem_slav-nb_elem_mpi*nb_proc
        idx_start = 1+(i_proc)*nb_elem_mpi
        idx_end = idx_start+nb_elem_mpi-1+nbr_elem_mpi*int((i_proc+1)/nb_proc)
        !write(*,*)"Proc : ", i_proc, "idx_start", idx_start, "idx_end", idx_end
        nb_el_slav_mpi = idx_end-idx_start+1
        AS_ALLOCATE(vi=v_appa_slav_mpi, size=nb_el_slav_mpi)
        v_appa_slav_mpi(:) = v_sdappa_slav(idx_start:idx_end)
        !write(*,*)"I_PROC = ", i_proc
        !write(*,*)"NB_ELEM_MPI = ",nb_el_slav_mpi, "LIST_ELEM_SLAV_MPI = ",v_appa_slav_mpi(:)
        !write(*,*)"NB_ELEM_SLAV = ", nb_elem_slav, "LIST_ELEM_SLAV = ",v_sdappa_slav(:)
        call aplcpg(mesh, newgeo, sdappa, i_zone, pair_tole, &
                    nb_elem_mast, v_sdappa_mast, nb_el_slav_mpi, v_appa_slav_mpi, &
                    nb_pair_zone, list_pair_zone, list_nbptit_zone, list_ptitsl_zone, &
                    list_ptitma_zone, list_ptgama_zone, int(i_proc), int(nb_proc), pair_method)
        AS_DEALLOCATE(vi=v_appa_slav_mpi)
    end do
!
! - Save pairing information in sdappa data structure
!
    !write(*,*)"Debut apstoc"
    call apstoc(ds_contact, nb_pair_zone, list_pair_zone, list_nbptit_zone, list_ptitsl_zone, &
                list_ptitma_zone, list_ptgama_zone)
    !write(*,*)"Fin apstoc"
!
! - Compute smooth normals at nodes
!

    if (l_smooth) then
        err_appa = 0
        call aplcno(mesh, newgeo, ds_contact%sdcont_defi, sdappa, err_appa)
        if (err_appa .eq. 1) then
            call mmbouc(ds_contact, 'Geom', 'Set_Error')
        else
            call mmbouc(ds_contact, 'Geom', 'Set_NoError')
        end if
    end if
!
end subroutine
