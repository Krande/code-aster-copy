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
subroutine cfmxr0_lac(mesh, ds_contact, ds_measure_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmrvai.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/detrsd.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Measure), optional, intent(inout) :: ds_measure_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Create CONT_ELEM for storing contact results
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ncmp = 5
    integer(kind=8), parameter :: nceld1 = 4
    integer(kind=8), parameter :: nceld2 = 4
    integer(kind=8), parameter :: nceld3 = 4
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nt_patch, nb_grel, nb_liel, nbliac
    integer(kind=8) :: i_liel, i_grel, patch_indx, elem_slav_nume, iret, vale_indx, decal
    real(kind=8) :: lagc, gapi, coef
    integer(kind=8) :: indi_cont
    character(len=19) :: celinr
    character(len=24) :: celinr_celv
    integer(kind=8) :: jv_celinr_celv
    character(len=24) :: celinr_celd
    integer(kind=8), pointer :: v_celinr_celd(:) => null()
    character(len=8)  :: ligrel_link_slav
    integer(kind=8), pointer :: v_ligrel_liel(:) => null()
    character(len=19) :: sdappa
    character(len=24) :: sdappa_coef
    real(kind=8), pointer :: v_sdappa_coef(:) => null()
    character(len=24) :: sdcont_stat
    integer(kind=8), pointer :: v_sdcont_stat(:) => null()
    character(len=24) :: sdcont_lagc
    real(kind=8), pointer :: v_sdcont_lagc(:) => null()
    character(len=24) :: sdappa_gapi
    real(kind=8), pointer :: v_sdappa_gapi(:) => null()
    integer(kind=8), pointer :: v_mesh_comapa(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'CONTACT5_10')
    end if
!
! - Get pairing datastructure
!
    sdappa = ds_contact%sdcont_solv(1:14)//'.APPA'
!
! - Get parameters
!
    nt_patch = ds_contact%nt_patch
!
! - Access to contact objects
!
    sdcont_stat = ds_contact%sdcont_solv(1:14)//'.STAT'
    sdcont_lagc = ds_contact%sdcont_solv(1:14)//'.LAGC'
    call jeveuo(sdcont_stat, 'L', vi=v_sdcont_stat)
    call jeveuo(sdcont_lagc, 'L', vr=v_sdcont_lagc)
    sdappa_coef = sdappa(1:19)//'.COEF'
    sdappa_gapi = sdappa(1:19)//'.GAPI'
    call jeveuo(sdappa_coef, 'L', vr=v_sdappa_coef)
    call jeveuo(sdappa_gapi, 'L', vr=v_sdappa_gapi)
!
! - Access to mesh
!
    call jeveuo(mesh//'.COMAPA', 'L', vi=v_mesh_comapa)
!
! - Get list of elements for slave surface (create in DEFI_CONTACT)
!
    ligrel_link_slav = ds_contact%ligrel_elem_slav
!
! - <CHELEM> for CONT_ELEM field
!
    celinr = ds_contact%field_cont_elem
!
! - Create new CONT_ELEM field
!
    call detrsd('CHAM_ELEM', celinr)
    call alchml(ligrel_link_slav//'.CONT.LIGRE', 'CONT_ELEM', 'CT_ELEM', 'V', celinr, &
                iret, ' ')
    ASSERT(iret .eq. 0)
!
! - Access to input field
!
    celinr_celd = celinr//'.CELD'
    celinr_celv = celinr//'.CELV'
    call jeveuo(celinr_celd, 'L', vi=v_celinr_celd)
    call jeveuo(celinr_celv, 'E', jv_celinr_celv)
    nb_grel = v_celinr_celd(2)
!
! - Fill input field
!
    nbliac = 0
    do i_grel = 1, nb_grel
        decal = v_celinr_celd(nceld1+i_grel)
        nb_liel = v_celinr_celd(decal+1)
        if (v_celinr_celd(decal+3) .ne. 0) then
            ASSERT(v_celinr_celd(decal+3) .eq. ncmp)
            call jeveuo(jexnum(ligrel_link_slav//'.CONT.LIGRE.LIEL', i_grel), &
                        'L', vi=v_ligrel_liel)
            do i_liel = 1, nb_liel
                elem_slav_nume = v_ligrel_liel(i_liel)
                patch_indx = v_mesh_comapa(elem_slav_nume)
                vale_indx = jv_celinr_celv-1+v_celinr_celd(decal+nceld2+nceld3*(i_liel-1)+4)
                if (patch_indx .ne. 0) then
                    lagc = v_sdcont_lagc(patch_indx)
                    gapi = v_sdappa_gapi(patch_indx)
                    coef = v_sdappa_coef(patch_indx)
                    indi_cont = v_sdcont_stat(patch_indx)
                    zr(vale_indx-1+1) = lagc
                    if (.not. isnan(gapi)) then
                        zr(vale_indx-1+2) = gapi
                    else
                        zr(vale_indx-1+2) = -1.d0
                    end if
                    zr(vale_indx-1+3) = indi_cont
                    zr(vale_indx-1+4) = coef
                    zr(vale_indx-1+5) = coef*lagc
                    if (indi_cont .eq. 1) then
                        nbliac = nbliac+1
                    end if
                end if
            end do
        end if
    end do
!
    if (present(ds_measure_)) then
        call nmrvai(ds_measure_, 'Cont_NCont', input_count=nbliac)
    end if
!
    call jedema()
end subroutine
