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
subroutine mmprel_lac(sdcont, mesh, model, slavElemLigr)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ajellt.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/getvtx.h"
#include "asterfort/jexnum.h"
#include "asterfort/jelira.h"
#include "asterfort/get_patchzi_num.h"
!
    character(len=8), intent(in) :: sdcont, model, mesh
    character(len=19), intent(in) :: slavElemLigr
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! LAC method - Create slave elements
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
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    character(len=24) :: sdcont_ptrdclac
    integer, pointer :: v_sdcont_mailco(:) => null()
    character(len=16) :: modeli, phenom, nmgrma
    integer :: jdecme, i_zone, i_sub_elem
    integer :: nb_cont_zone, model_ndim, nb_cont_elem, nt_elem_slav, nb_dcl_zi
    integer :: i_elem_slav, elem_slav_indx, elem_slav_nume, nb_elem_slav
    integer :: jdecpa, nb_patch, i_patch, patch_type, jngrma, nupatch_zi
    integer :: nt_sub_elem, nb_sub_elem, nb_list_elem, nb_grma, n1b
    character(len=24) :: list_elem
    integer, pointer :: v_list_elem(:) => null()
    integer, pointer :: v_mesh_lpatch(:) => null()
    integer, pointer :: v_mesh_patch(:) => null()
    integer, pointer :: vi_ptrdclac(:) => null()
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
    sdcont_ptrdclac = sdcont_defi(1:16)//'.PTRDCLC'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)

! - Access to mesh (patches)
    call jeveuo(jexnum(mesh//'.PATCH', 1), 'L', vi=v_mesh_lpatch)
    nmgrma = '&&OP0060.NMMA'

! - Parameters
    model_ndim = cfdisi(sdcont_defi, 'NDIM')
    nb_cont_elem = cfdisi(sdcont_defi, 'NMACO')
    nt_elem_slav = cfdisi(sdcont_defi, 'NTMAE')
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')

! - Check compatiblity DECOUPE_LAC<=>DEFI_CONTACT
    call jelira(mesh//'.PTRNOMPAT', 'LONMAX', nb_dcl_zi)
    if (nb_dcl_zi .ne. nb_cont_zone) then
        call utmess('F', 'CONTACT2_17')
    end if

! - Create pointer index DECOUPE_LAC<=>DEFI_CONTACT
    call wkvect(sdcont_ptrdclac, 'G V I', nb_cont_zone, vi=vi_ptrdclac)

! - Create list of slave elements
    nb_list_elem = 4
    call wkvect(list_elem, 'V V I', nb_list_elem, vi=v_list_elem)

! - Set list of slave elements
    nt_sub_elem = 0
    do i_zone = 1, nb_cont_zone
!
! ----- Get current patches
!
        call getvtx('ZONE', 'GROUP_MA_ESCL', iocc=1, nbval=0, nbret=nb_grma)
        ASSERT(nb_grma .eq. -1)
        call wkvect(nmgrma, 'V V K24', -nb_grma, jngrma)
        call getvtx('ZONE', 'GROUP_MA_ESCL', iocc=i_zone, nbval=-nb_grma, &
                    vect=zk24(jngrma), nbret=n1b)
        ASSERT(n1b .eq. -nb_grma)

        call get_patchzi_num(mesh, zk24(jngrma), nupatch_zi)
        vi_ptrdclac(i_zone) = nupatch_zi

        nb_patch = v_mesh_lpatch(2*(nupatch_zi-1)+2)
        jdecpa = v_mesh_lpatch(2*(nupatch_zi-1)+1)

! ----- Acces to slave elements in zone
        nb_elem_slav = mminfi(sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(sdcont_defi, 'JDECME', i_zone)
        nb_sub_elem = 0
        i_elem_slav = 1
        do i_patch = 1, nb_patch
! --------- Get current patch
            call jeveuo(jexnum(mesh//'.PATCH', jdecpa+i_patch-1), 'L', vi=v_mesh_patch)
            patch_type = v_mesh_patch(1)
!
            if (patch_type .eq. 7) then
! ------------- 2D - TRIA3/SEG2
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2DB'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+2
            elseif (patch_type .eq. 12) then
! ------------- 2D - QUAD4/SEG2
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2DT'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2DB'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+3
            elseif (patch_type .eq. 9 .or. patch_type .eq. 14) then
! ------------- 2D - TRIA6/SEG3
! ------------- 2D - QUAD8/SEG3
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_2D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+1
            elseif (patch_type .eq. 18 .or. patch_type .eq. 19) then
! ------------- 3D - TETRA4/TRIA3
! ------------- 3D - TETRA10/TRIA6
                do i_sub_elem = 1, 3
                    elem_slav_indx = jdecme+i_elem_slav+i_sub_elem-1
                    elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                    v_list_elem(i_sub_elem) = elem_slav_nume
                end do
                nb_list_elem = 3
                i_elem_slav = i_elem_slav+3
                modeli = 'CONT_LAC_SL_3D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+3
            elseif (patch_type .eq. 27) then
! ------------- 3D - HEXA27/QUAD9
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_3D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+1
            else if (patch_type .eq. 25 .or. patch_type .eq. 26) then
! ------------- 3D - HEXA8/QUAD4
! ------------- 3D - HEXA20/QUAD8
                do i_sub_elem = 1, 4
                    elem_slav_indx = jdecme+i_elem_slav+i_sub_elem-1
                    elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                    v_list_elem(i_sub_elem) = elem_slav_nume
                end do
                nb_list_elem = 4
                i_elem_slav = i_elem_slav+4
                modeli = 'CONT_LAC_SL_3D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                elem_slav_indx = jdecme+i_elem_slav
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                v_list_elem(1) = elem_slav_nume
                nb_list_elem = 1
                i_elem_slav = i_elem_slav+1
                modeli = 'CONT_LAC_SL_3DB'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+5
            else if (patch_type .eq. 125 .or. patch_type .eq. 126) then
! ------------- 3D - HEXA8/QUAD4 => PYRA5/TRIA3
! ------------- 3D - HEXA20/QUAD8 => PYRA13/TRIA6
                do i_sub_elem = 1, 4
                    elem_slav_indx = jdecme+i_elem_slav+i_sub_elem-1
                    elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                    v_list_elem(i_sub_elem) = elem_slav_nume
                end do
                nb_list_elem = 4
                i_elem_slav = i_elem_slav+4
                modeli = 'CONT_LAC_SL_3D'
                call ajellt(slavElemLigr, mesh, nb_list_elem, list_elem, ' ', &
                            phenom, modeli, 0, ' ')
                nb_sub_elem = nb_sub_elem+4
            else
                ASSERT(.false.)
            end if
        end do
        call jedetr(nmgrma)
        nt_sub_elem = nt_sub_elem+nb_sub_elem
    end do
    if (nt_sub_elem .ne. nt_elem_slav) then
        call utmess('F', 'CONTACT4_5')
    end if
!
    call jedetr(list_elem)
!
end subroutine
