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

subroutine mmmreg(mesh, ds_contact, v_disp_cumu, nb_dof, &
                  v_slav_slide, v_mast_slide)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/mmelty.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfm.h"
#include "asterfort/mmnonf.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: nb_dof
    real(kind=8), pointer :: v_disp_cumu(:)
    real(kind=8), pointer :: v_slav_slide(:)
    real(kind=8), pointer :: v_mast_slide(:)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Post-treatment
!
! Continue method - Prepare post-treatment fields / compute slides
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
! In  v_disp_cumu      : pointer to displacement increment from beginning of current time
! In  nb_dof           : number of dof by node
! IO  v_slav_slide     : slides on slave side
! IO  v_mast_slide     : slides on master side
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_cont_poin, i_zone, i_elem_slav, i_poin_elem, i_node
    integer(kind=8) :: nb_poin_elem, nb_elem_slav, nb_cont_zone, model_ndim, nt_cont_poin
    integer(kind=8) :: node_slav_nume, elem_slav_nume
    integer(kind=8) :: node_mast_nume, elem_mast_nume
    integer(kind=8) :: elem_slav_nbnode, elem_mast_nbnode
    integer(kind=8) :: elem_slav_indx, jdecme
    real(kind=8) :: ksipc1, ksipc2, ksipr1, ksipr2
    integer(kind=8) :: ztabf
    real(kind=8) :: disp_mast(3), disp_slav(3)
    real(kind=8) :: tau1(3), tau2(3), ff(9)
    character(len=8) :: elem_slav_type, elem_mast_type
    aster_logical :: l_veri
    character(len=24) :: sdcont_tabfin
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    integer(kind=8), pointer :: v_mesh_connex(:) => null()
    integer(kind=8), pointer :: v_mesh_loncum(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Access to mesh
!
    call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', vi=v_mesh_loncum)
    call jeveuo(mesh(1:8)//'.CONNEX', 'L', vi=v_mesh_connex)
!
! - Acces to contact objects
!
    ztabf = cfmmvd('ZTABF')
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    call jeveuo(sdcont_tabfin, 'L', vr=v_sdcont_tabfin)
!
! - Get parameters
!
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    nt_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
!
! - Compute slides
!
    i_cont_poin = 1
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters of zone
!
        l_veri = mminfl(ds_contact%sdcont_defi, 'VERIF', i_zone)
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', i_zone)
!
! ----- No computation: no contact point
!
        if (l_veri) then
            goto 25
        end if
!
! ----- Loop on slave elements
!
        do i_elem_slav = 1, nb_elem_slav
!
! --------- Slave element index in contact datastructure
!
            elem_slav_indx = jdecme+i_elem_slav
!
! --------- Number of integration points on element
!
            call mminfm(elem_slav_indx, ds_contact%sdcont_defi, 'NPTM', nb_poin_elem)
!
! --------- Loop on integration points
!
            do i_poin_elem = 1, nb_poin_elem
!
! ------------- Get parameters
!
                elem_slav_nume = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+2))
                elem_mast_nume = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+3))
                ksipr1 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+6)
                ksipr2 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+7)
                ksipc1 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+4)
                ksipc2 = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+5)
!
! ------------- Get local basis
!
                tau1(1) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+8)
                tau1(2) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+9)
                tau1(3) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+10)
                tau2(1) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+11)
                tau2(2) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+12)
                tau2(3) = v_sdcont_tabfin(ztabf*(i_cont_poin-1)+13)
!
! ------------- Slave node displacement
!
                call mmelty(mesh, elem_slav_nume, elem_slav_type, elem_slav_nbnode)
                call mmnonf(model_ndim, elem_slav_nbnode, elem_slav_type, ksipc1, ksipc2, &
                            ff)
                disp_slav(1:3) = 0.d0
                do i_node = 1, elem_slav_nbnode
                    node_slav_nume = v_mesh_connex(1+v_mesh_loncum(elem_slav_nume)+i_node-2)
                    disp_slav(1) = disp_slav(1)+ &
                                   v_disp_cumu(nb_dof*(node_slav_nume-1)+1)*ff(i_node)
                    disp_slav(2) = disp_slav(2)+ &
                                   v_disp_cumu(nb_dof*(node_slav_nume-1)+2)*ff(i_node)
                    disp_slav(3) = disp_slav(3)+ &
                                   v_disp_cumu(nb_dof*(node_slav_nume-1)+3)*ff(i_node)
                end do
!
! ------------- Master node displacement
!
                call mmelty(mesh, elem_mast_nume, elem_mast_type, elem_mast_nbnode)
                call mmnonf(model_ndim, elem_mast_nbnode, elem_mast_type, ksipr1, ksipr2, &
                            ff)
                disp_mast(1:3) = 0.d0
                do i_node = 1, elem_mast_nbnode
                    node_mast_nume = v_mesh_connex(1+v_mesh_loncum(elem_mast_nume)+i_node-2)
                    disp_mast(1) = disp_mast(1)+ &
                                   v_disp_cumu(nb_dof*(node_mast_nume-1)+1)*ff(i_node)
                    disp_mast(2) = disp_mast(2)+ &
                                   v_disp_cumu(nb_dof*(node_mast_nume-1)+2)*ff(i_node)
                    if (model_ndim .eq. 3) then
                        disp_mast(3) = disp_mast(3)+ &
                                       v_disp_cumu(nb_dof*(node_mast_nume-1)+3)*ff(i_node)
                    end if
                end do
!
! ------------- Compute slides
!
                if (model_ndim .eq. 3) then
                    v_slav_slide(2*(i_cont_poin-1)+1) = disp_slav(1)*tau1(1)+ &
                                                        disp_slav(2)*tau1(2)+ &
                                                        disp_slav(3)*tau1(3)
                    v_slav_slide(2*(i_cont_poin-1)+2) = disp_slav(1)*tau2(1)+ &
                                                        disp_slav(2)*tau2(2)+ &
                                                        disp_slav(3)*tau2(3)
                    v_mast_slide(2*(i_cont_poin-1)+1) = disp_mast(1)*tau1(1)+ &
                                                        disp_mast(2)*tau1(2)+ &
                                                        disp_mast(3)*tau1(3)
                    v_mast_slide(2*(i_cont_poin-1)+2) = disp_mast(1)*tau2(1)+ &
                                                        disp_mast(2)*tau2(2)+ &
                                                        disp_mast(3)*tau2(3)
                else if (model_ndim .eq. 2) then
                    v_slav_slide(i_cont_poin) = disp_slav(1)*tau1(1)+disp_slav(2)*tau1(2)
                    v_mast_slide(i_cont_poin) = disp_mast(1)*tau1(1)+disp_mast(2)*tau1(2)
                else
                    ASSERT(.false.)
                end if
!
! ------------- Next point
!
                i_cont_poin = i_cont_poin+1
            end do
        end do
25      continue
    end do
!
    call jedema()
end subroutine
