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

subroutine mmprel(sdcont, mesh, model, slavElemLigr)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ajellt.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: sdcont, model, mesh
    character(len=19), intent(in) :: slavElemLigr
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Continue method - Create slave elements
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  model            : name of model
! In  mesh             : name of mesh
! In  slavElemLigr     : LIGREL for virtual elements (slave side)
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_frot_zone, l_veri
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    integer, pointer :: v_sdcont_mailco(:) => null()
    character(len=16) :: modeli, phenom
    integer :: jdecme, i_zone
    integer :: nb_cont_zone, model_ndim, nb_cont_elem, nt_elem_slav
    integer :: i_elem_slav, elem_slav_idx, elem_slav_nume, nb_elem_slav
    aster_logical :: l_verif_all
    character(len=24) :: list_elem
    integer, pointer :: v_list_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    list_elem = '&&MMPREL.LISTE_MAILLES'
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenom)
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
!
! - Parameters
!
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nb_cont_elem = cfdisi(sdcont_defi, 'NMACO')
    nt_elem_slav = cfdisi(sdcont_defi, 'NTMAEC')
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    l_verif_all = cfdisl(sdcont_defi, 'ALL_VERIF')
!
! - Add elements
!
    if (.not. l_verif_all) then
!
! ----- Create list of slave elements
!
        call wkvect(list_elem, 'V V I', nt_elem_slav, vi=v_list_elem)
!
! ----- Set list of slave elements
!
        do i_zone = 1, nb_cont_zone
!
! --------- Type of model
!
            l_frot_zone = mminfl(sdcont_defi, 'FROTTEMENT_ZONE', i_zone)
            l_veri = mminfl(sdcont_defi, 'VERIF', i_zone)
            if (model_ndim .eq. 2) then
                if (l_frot_zone) then
                    modeli = 'FRIC_SL_2D'
                else
                    modeli = 'CONT_SL_2D'
                end if
            else if (model_ndim .eq. 3) then
                if (l_frot_zone) then
                    modeli = 'FRIC_SL_3D'
                else
                    modeli = 'CONT_SL_3D'
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
!
! --------- Type of model
!
            if (.not. l_veri) then
                nb_elem_slav = mminfi(sdcont_defi, 'NBMAE', i_zone)
                jdecme = mminfi(sdcont_defi, 'JDECME', i_zone)
                ASSERT(nb_elem_slav .le. nt_elem_slav)
                do i_elem_slav = 1, nb_elem_slav
                    elem_slav_idx = jdecme+i_elem_slav
                    elem_slav_nume = v_sdcont_mailco(elem_slav_idx)
                    v_list_elem(i_elem_slav) = elem_slav_nume
                end do
                call ajellt(slavElemLigr, mesh, nb_elem_slav, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
            end if
        end do
        call jedetr(list_elem)
    end if
!
end subroutine
