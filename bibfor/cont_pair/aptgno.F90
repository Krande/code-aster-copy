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

subroutine aptgno(sdappa, mesh, sdcont_defi)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfcald.h"
#include "asterfort/jeveuo.h"
#include "asterfort/cfdisi.h"
#include "asterfort/aptgnn.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfr.h"
#include "asterfort/infdbg.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=19), intent(in) :: sdappa
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: sdcont_defi
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing
!
! Compute tangents at each node (average)
!
! --------------------------------------------------------------------------------------------------
!
! In  sdappa           : name of pairing datastructure
! In  mesh             : name of mesh
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_cont_zone, model_ndim
    integer(kind=8) :: i_zone, norm_type
    integer(kind=8) :: jdecnm, nb_node_mast
    integer(kind=8) :: jdecne, nb_node_slav
    aster_logical :: apcald
    real(kind=8) :: norm_vect(3)
    real(kind=8), pointer :: v_sdappa_tgno(:) => null()
    character(len=24) :: sdappa_tgno
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<APPARIEMENT> ...... TANGENTES SUR LES NOEUDS'
    end if
!
! - Get parameters
!
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    sdappa_tgno = sdappa(1:19)//'.TGNO'
    call jeveuo(sdappa_tgno, 'E', vr=v_sdappa_tgno)
    v_sdappa_tgno(:) = 0.d0
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters on current zone - Master
!
        nb_node_mast = mminfi(sdcont_defi, 'NBNOM', i_zone)
        jdecnm = mminfi(sdcont_defi, 'JDECNM', i_zone)
        norm_type = mminfi(sdcont_defi, 'VECT_MAIT', i_zone)
        norm_vect(1:3) = 0.d0
        if (norm_type .ne. 0) then
            norm_vect(1) = mminfr(sdcont_defi, 'VECT_MAIT_DIRX', i_zone)
            norm_vect(2) = mminfr(sdcont_defi, 'VECT_MAIT_DIRY', i_zone)
            norm_vect(3) = mminfr(sdcont_defi, 'VECT_MAIT_DIRZ', i_zone)
        end if
!
! ----- Compute tangents at each node by smoothing - On current zone/Master
!
        apcald = cfcald(sdcont_defi, i_zone, 'MAIT')
        if (apcald) then
            call aptgnn(sdappa, mesh, sdcont_defi, model_ndim, jdecnm, &
                        nb_node_mast, norm_type, norm_vect)
        end if
!
! ----- Parameters on current zone - Slave
!
        nb_node_slav = mminfi(sdcont_defi, 'NBNOE', i_zone)
        jdecne = mminfi(sdcont_defi, 'JDECNE', i_zone)
        norm_type = mminfi(sdcont_defi, 'VECT_ESCL', i_zone)
        norm_vect(1:3) = 0.d0
        if (norm_type .ne. 0) then
            norm_vect(1) = mminfr(sdcont_defi, 'VECT_ESCL_DIRX', i_zone)
            norm_vect(2) = mminfr(sdcont_defi, 'VECT_ESCL_DIRY', i_zone)
            norm_vect(3) = mminfr(sdcont_defi, 'VECT_ESCL_DIRZ', i_zone)
        end if
!
! ----- Compute tangents at each node by smoothing - On current zone/salve
!
        apcald = cfcald(sdcont_defi, i_zone, 'ESCL')
        if (apcald) then
            call aptgnn(sdappa, mesh, sdcont_defi, model_ndim, jdecne, &
                        nb_node_slav, norm_type, norm_vect)
        end if
    end do
!
end subroutine
