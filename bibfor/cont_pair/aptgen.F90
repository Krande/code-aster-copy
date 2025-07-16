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

subroutine aptgen(sdappa, mesh, sdcont_defi, newgeo, err_appa)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisr.h"
#include "asterfort/aptgem.h"
#include "asterfort/mminfi.h"
#include "asterfort/cfcald.h"
#include "asterfort/infdbg.h"
#include "asterfort/jelira.h"
#include "asterfort/jerazo.h"
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
! Compute tangents at each node for each element
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
    integer(kind=8) :: nb_cont_zone, model_ndim
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_zone, length
    integer(kind=8) :: jdecmm, nb_elem_mast
    integer(kind=8) :: jdecme, nb_elem_slav
    character(len=4) :: zone_type
    aster_logical :: apcald
    real(kind=8) :: epsi_maxi
    integer(kind=8) :: iter_maxi
    character(len=24) :: sdappa_tgel
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<APPARIEMENT> ...... TANGENTES SUR LES NOEUDS PAR ELEMENT (ELNO)'
    end if
!
! - Get parameters
!
    epsi_maxi = cfdisr(sdcont_defi, 'PROJ_NEWT_RESI')
    iter_maxi = cfdisi(sdcont_defi, 'PROJ_NEWT_ITER')
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    sdappa_tgel = sdappa(1:19)//'.TGEL'
    length = 0
    call jelira(sdappa_tgel, 'LONT', length)
    call jerazo(sdappa_tgel, length, 1)
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters on current zone - Master
!
        nb_elem_mast = mminfi(sdcont_defi, 'NBMAM', i_zone)
        jdecmm = mminfi(sdcont_defi, 'JDECMM', i_zone)
        zone_type = 'MAIT'
!
! ----- Compute tangents at each node for each element - Master
!
        apcald = cfcald(sdcont_defi, i_zone, 'MAIT')
        if (apcald) then
            call aptgem(sdappa, mesh, newgeo, sdcont_defi, model_ndim, &
                        i_zone, zone_type, epsi_maxi, jdecmm, &
                        nb_elem_mast, err_appa)
        end if
!
! ----- Parameters on current zone - Slave
!
        nb_elem_slav = mminfi(sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(sdcont_defi, 'JDECME', i_zone)
        zone_type = 'ESCL'
!
! ----- Compute tangents at each node for each element - Slave
!
        apcald = cfcald(sdcont_defi, i_zone, 'ESCL')
        if (apcald) then
            call aptgem(sdappa, mesh, newgeo, sdcont_defi, model_ndim, &
                        i_zone, zone_type, epsi_maxi, jdecme, &
                        nb_elem_slav, err_appa)
        end if
    end do
!
end subroutine
