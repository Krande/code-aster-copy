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
subroutine nonlinDSContactInit(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/armin.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/matdis.h"
#include "asterfort/utmess.h"
#include "asterfort/lac_rela.h"
#include "asterfort/wkvect.h"
#include "asterfort/jeveuo.h"
#include "asterfort/isParallelMesh.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Contact management
!
! Initializations for contact management
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: cont_form
    character(len=8) :: sdcont
    character(len=24) :: iden_rela
    aster_logical :: l_cont, l_unil
    aster_logical :: l_form_disc, l_form_cont, l_form_lac
    aster_logical :: l_all_verif, l_iden_rela
    aster_logical :: l_unil_pena, l_inte_node
    integer(kind=8) :: nt_patch
    integer(kind=8) :: i_exist
    character(len=8), pointer :: v_load_type(:) => null()
    character(len=24) :: sdcont_paraci
    integer(kind=8), pointer :: v_sdcont_paraci(:) => null()
    character(len=3) :: matd
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
!
! - Initializations
!
    l_cont = ASTER_FALSE
    l_unil = ASTER_FALSE
    l_form_disc = ASTER_FALSE
    l_form_cont = ASTER_FALSE
    l_form_lac = ASTER_FALSE
    l_iden_rela = ASTER_FALSE
    iden_rela = '&&CFMXR0.IDEN_RELA'
!
    if (ds_contact%l_contact) then
! ----- Print
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_3')
        end if

! ----- Forbiden for a ParallelMesh
        ASSERT(.not. isParallelMesh(mesh))

! ----- Datastructure from DEFI_CONTACT
        sdcont = ds_contact%sdcont
        sdcont_paraci = sdcont(1:8)//'.PARACI'
        ds_contact%sdcont_defi = sdcont(1:8)//'.CONTACT'
        ds_contact%sdunil_defi = sdcont(1:8)//'.UNILATE'
        call jeveuo(sdcont_paraci, 'E', vi=v_sdcont_paraci)

! ----- Contact formulation
        cont_form = cfdisi(ds_contact%sdcont_defi, 'FORMULATION')
        ASSERT(cont_form .ge. 1 .and. cont_form .le. 5)
        l_form_disc = cont_form .eq. 1
        l_form_cont = cont_form .eq. 2
        l_unil = cont_form .eq. 4
        l_form_lac = cont_form .eq. 5
        l_cont = cont_form .ne. 4
        l_all_verif = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
        l_inte_node = cfdisl(ds_contact%sdcont_defi, 'ALL_INTEG_NOEUD')
        l_iden_rela = ASTER_FALSE

! ----- Fields for CONT_NOEU
        if (l_form_cont .or. l_form_disc) then
            ds_contact%field_cont_node = '&&CFMXR0.CNOINR'
            ds_contact%fields_cont_node = '&&CFMXR0.CNSINR'
            ds_contact%field_cont_perc = '&&CFMXR0.CNSPER'
            if (l_inte_node) then
                ds_contact%l_cont_node = ASTER_TRUE
            end if
        end if

! ----- Fields for CONT_ELEM
        if (l_form_lac) then
            ds_contact%field_cont_elem = '&&CFMXR0.CEOINR'
            ds_contact%fields_cont_elem = '&&CFMXR0.CESINR'
            ds_contact%l_cont_elem = ASTER_TRUE
        end if

! ----- Special for contact stabilization with elastic matrix
        if (l_form_lac .or. l_form_cont) then
            ds_contact%sContStab = v_sdcont_paraci(31)
            if (ds_contact%sContStab .gt. 0) then
                ds_contact%lContStab = ASTER_TRUE
            end if
        end if

! ----- Special for discrete contact
        if (l_form_disc) then
            ds_contact%nume_dof_frot = '&&CFMXSD.NUMDF'
            call jeexin(sdcont(1:8)//'.CHME.LIGRE.LGRF', i_exist)
            ds_contact%l_dof_rela = i_exist .gt. 0
            if (i_exist .gt. 0) then
                ds_contact%ligrel_dof_rela = sdcont
            end if
        end if

! ----- Special for continue contact
        if (l_form_cont) then
            call matdis(matd)
            if (matd .eq. 'OUI') then
                call utmess('F', 'MECANONLINE_6')
            end if
            ds_contact%field_input = ds_contact%sdcont_solv(1:14)//'.CHML'
            ds_contact%l_elem_slav = ASTER_TRUE
            ds_contact%ligrel_elem_slav = sdcont
            ds_contact%l_elem_cont = ASTER_TRUE
            ds_contact%ligrel_elem_cont = '&&LIGRCF.CONT.LIGRE'
            call wkvect(ds_contact%ligrel_elem_cont(1:8)//'.TYPE', 'V V K8', 1, vk8=v_load_type)
            v_load_type(1) = 'ME'
            ds_contact%it_cycl_maxi = 6
        end if

! ----- Special for LAC contact
        if (l_form_lac) then
            ds_contact%field_input = ds_contact%sdcont_solv(1:14)//'.CHML'
            ds_contact%l_elem_slav = ASTER_TRUE
            ds_contact%ligrel_elem_slav = sdcont
            ds_contact%l_elem_cont = ASTER_TRUE
            ds_contact%ligrel_elem_cont = '&&LIGRCF.CONT.LIGRE'
            call wkvect(ds_contact%ligrel_elem_cont(1:8)//'.TYPE', 'V V K8', 1, vk8=v_load_type)
            v_load_type(1) = 'ME'
            call jelira(mesh//'.PATCH', 'NUTIOC', nt_patch)
            nt_patch = nt_patch-1
            ds_contact%nt_patch = nt_patch
            call lac_rela(mesh, ds_contact, iden_rela, l_iden_rela)
            ds_contact%arete_min = armin(mesh)
        end if

! ----- Identity relation
        ds_contact%l_iden_rela = l_iden_rela
        if (l_iden_rela) then
            ds_contact%iden_rela = iden_rela
        else
            ds_contact%iden_rela = ' '
        end if

! ----- Flag for (re) numbering
        if (l_form_cont) then
            if (l_all_verif) then
                ds_contact%l_renumber = ASTER_FALSE
            else
                ds_contact%l_renumber = ASTER_TRUE
            end if
        end if
        if (l_form_lac) then
            ds_contact%l_renumber = ASTER_TRUE
        end if

! ----- Flag for pairing
        if (l_form_disc) then
            ds_contact%l_pair = ASTER_TRUE
            ds_contact%l_first_geom = ASTER_TRUE
        end if

! ----- Special for UNIL contact
        if (l_unil) then
            l_unil_pena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
            if (l_unil_pena) then
                ds_contact%nume_dof_unil = '&&NMASUN.NUME'
            end if
        end if

! ----- Save parameters
        ds_contact%l_meca_cont = l_cont
        ds_contact%l_meca_unil = l_unil
        ds_contact%l_form_cont = l_form_cont
        ds_contact%l_form_disc = l_form_disc
        ds_contact%l_form_lac = l_form_lac
    end if
!
end subroutine
