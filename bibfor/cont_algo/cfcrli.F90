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

subroutine cfcrli(mesh, nume_dof, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfecrd.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnomm.h"
#include "asterfort/cfverd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/posddl.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Create datastructures for linear relations
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_cont_node, model_ndim, nt_cont_poin
    character(len=24) :: sdcont_apcoef, sdcont_apcofr
    integer(kind=8) :: jv_sdcont_apcoef, jv_sdcont_apcofr
    character(len=24) :: sdcont_ddlco, sdcont_nbddl
    integer(kind=8), pointer :: v_sdcont_ddlco(:) => null()
    integer(kind=8), pointer :: v_sdcont_nbddl(:) => null()
    character(len=24) :: sdcont_apddl, sdcont_approj
    integer(kind=8) :: jv_sdcont_apddl, jv_sdcont_approj
    character(len=24) :: sdcont_coco
    integer(kind=8) :: jv_sdcont_coco
    character(len=24) :: sdcont_appoin, sdcont_numlia
    integer(kind=8) :: jv_sdcont_appoin, jv_sdcont_numlia
    character(len=19) :: sdcont_liac, sdcont_liot
    integer(kind=8) :: jv_sdcont_liac, jv_sdcont_liot
    integer(kind=8) :: dof_nume, node_nume
    integer(kind=8) :: i_node, i_dof, node_indx
    character(len=8) :: node_name
    integer(kind=8) :: nb_equa, nt_slav_nmax, nb_dof
    aster_logical :: l_frot
    integer(kind=8) :: zcoco
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... CREATION DE LA SD POUR LES LIAISONS LINEAIRES'
    end if
!
! - Number of equations
!
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=nb_equa)
!
! - Get contact parameters
!
    zcoco = cfmmvd('ZCOCO')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
    nb_cont_node = cfdisi(ds_contact%sdcont_defi, 'NNOCO')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nt_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
    nt_slav_nmax = nt_cont_poin
!
! - Check model dimension
!
    if (model_ndim .eq. 2) then
        call cfverd(mesh, nume_dof, ds_contact%sdcont_defi)
    end if
!
! - Datastructure for number of dof by node
!
    sdcont_nbddl = ds_contact%sdcont_solv(1:14)//'.NBDDL'
    call wkvect(sdcont_nbddl, 'V V I', nb_cont_node+1, vi=v_sdcont_nbddl)
    v_sdcont_nbddl(1) = 0
    nb_dof = 0
    do i_node = 1, nb_cont_node
        nb_dof = nb_dof+model_ndim
        v_sdcont_nbddl(i_node+1) = nb_dof
    end do
!
! - Datastructure for index of dof
!
    sdcont_ddlco = ds_contact%sdcont_solv(1:14)//'.DDLCO'
    call wkvect(sdcont_ddlco, 'V V I', nb_dof, vi=v_sdcont_ddlco)
    do i_node = 1, nb_cont_node
        node_indx = i_node
        call cfnomm(mesh, ds_contact%sdcont_defi, 'NOEU', node_indx, node_name)
        i_dof = v_sdcont_nbddl(i_node)+1
        call posddl('NUME_DDL', nume_dof, node_name, 'DX', node_nume, dof_nume)
        ASSERT(dof_nume .gt. 0)
        v_sdcont_ddlco(i_dof) = dof_nume
        call posddl('NUME_DDL', nume_dof, node_name, 'DY', node_nume, dof_nume)
        ASSERT(dof_nume .gt. 0)
        v_sdcont_ddlco(i_dof+1) = dof_nume
        if (model_ndim .eq. 3) then
            call posddl('NUME_DDL', nume_dof, node_name, 'DZ', node_nume, dof_nume)
            ASSERT(dof_nume .gt. 0)
            v_sdcont_ddlco(i_dof+2) = dof_nume
        end if
    end do
!
! - Pointer for dof
!
    sdcont_apddl = ds_contact%sdcont_solv(1:14)//'.APDDL'
    call wkvect(sdcont_apddl, 'V V I', 30*nt_cont_poin, jv_sdcont_apddl)
!
! - Number of dof by relation
!
    sdcont_appoin = ds_contact%sdcont_solv(1:14)//'.APPOIN'
    call wkvect(sdcont_appoin, 'V V I', nt_cont_poin+1, jv_sdcont_appoin)
!
! - Number of point and slave node for relation
!
    sdcont_numlia = ds_contact%sdcont_solv(1:14)//'.NUMLIA'
    call wkvect(sdcont_numlia, 'V V I', 4*nt_cont_poin, jv_sdcont_numlia)
!
! - Management of relations
!
    sdcont_coco = ds_contact%sdcont_solv(1:14)//'.COCO'
    call wkvect(sdcont_coco, 'V V I', zcoco, jv_sdcont_coco)
    call cfecrd(ds_contact%sdcont_solv, 'NDIM', model_ndim)
    call cfecrd(ds_contact%sdcont_solv, 'NEQ', nb_equa)
    call cfecrd(ds_contact%sdcont_solv, 'NESMAX', nt_slav_nmax)
!
! - PIVOT NULS
!
    sdcont_liot = ds_contact%sdcont_solv(1:14)//'.LIOT'
    call wkvect(sdcont_liot, 'V V I', 4*nt_cont_poin+4, jv_sdcont_liot)
!
! - Activated relations
!
    sdcont_liac = ds_contact%sdcont_solv(1:14)//'.LIAC'
    call wkvect(sdcont_liac, 'V V I', 3*nt_cont_poin+1, jv_sdcont_liac)
!
! - Coefficient for relations
!
    sdcont_apcoef = ds_contact%sdcont_solv(1:14)//'.APCOEF'
    call wkvect(sdcont_apcoef, 'V V R', 30*nt_cont_poin, jv_sdcont_apcoef)
    if (l_frot) then
        sdcont_apcofr = ds_contact%sdcont_solv(1:14)//'.APCOFR'
        call wkvect(sdcont_apcofr, 'V V R', 60*nt_cont_poin, jv_sdcont_apcofr)
    end if
!
! - Projection coordinates
!
    sdcont_approj = ds_contact%sdcont_solv(1:14)//'.APPROJ'
    call wkvect(sdcont_approj, 'V V R', 3*nt_cont_poin, jv_sdcont_approj)
!
end subroutine
