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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cfcrsd(mesh, nume_dof, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/cfcrje.h"
#include "asterfort/cfcrli.h"
#include "asterfort/cfcrma.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnum.h"
#include "asterfort/vtcreb.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Create datastructures for DISCRETE methods
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  nume_dof         : name of numbering object (NUME_DDL)
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ztacf
    integer(kind=8) :: model_ndim, nt_cont_poin, nb_cont_node, nb_equa
    integer(kind=8) :: nbliai, ii
    character(len=19) :: sdcont_mu, sdcont_atmu, sdcont_afmu
    integer(kind=8) :: jv_sdcont_mu, jv_sdcont_atmu, jv_sdcont_afmu
    character(len=19) :: sdcont_copo
    real(kind=8), pointer :: v_sdcont_copo(:) => null()
    character(len=19) :: sdcont_del0, sdcont_ddel, sdcont_delc, sdcont_dep0, sdcont_depc
    character(len=19) :: sdcont_cm1a, sdcont_enat, sdcont_fro1, sdcont_fro2
    character(len=19) :: sdcont_cin0
    character(len=19) :: sdcont_sgdm, sdcont_sgdp, sdcont_sgpm, sdcont_sgpp
    integer(kind=8) :: jv_sdcont_sgdm, jv_sdcont_sgdp, jv_sdcont_sgpm, jv_sdcont_sgpp
    character(len=19) :: sdcont_dire, sdcont_mum, sdcont_secm
    integer(kind=8) :: jv_sdcont_dire, jv_sdcont_mum
    character(len=19) :: sdcont_pcrs, sdcont_pcdr, sdcont_pcuu
    integer(kind=8) :: jv_sdcont_pcrs, jv_sdcont_pcdr, jv_sdcont_pcuu
    character(len=19) :: sdcont_svm0, sdcont_svmu
    integer(kind=8) :: jv_sdcont_svm0, jv_sdcont_svmu
    integer(kind=8) :: nbcm1a, nbenat, nbfro1, nbfro2
    character(len=24) :: sdcont_rea1, sdcont_rea2
    character(len=24) :: sdcont_tacfin, sdcont_tangco
    integer(kind=8) :: jv_sdcont_tacfin, jv_sdcont_tangco
    aster_logical :: l_frot, l_pena_cont, l_pena_frot, l_matr_cont, l_gcp, l_pre_cond
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'CONTACT5_4')
    end if
!
! - Get contact parameters
!
    ztacf = cfmmvd('ZTACF')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nt_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
    nb_cont_node = cfdisi(ds_contact%sdcont_defi, 'NNOCO')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
    l_pena_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_PENA')
    l_pena_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_PENA')
    l_matr_cont = cfdisl(ds_contact%sdcont_defi, 'MATR_CONT')
    l_gcp = cfdisl(ds_contact%sdcont_defi, 'CONT_GCP')
    l_pre_cond = cfdisl(ds_contact%sdcont_defi, 'PRE_COND_DIRICHLET')
!
! - For geometric loop
!
    sdcont_rea1 = ds_contact%sdcont_solv(1:14)//'.REA1'
    sdcont_rea2 = ds_contact%sdcont_solv(1:14)//'.REA2'
    call vtcreb(sdcont_rea1, 'V', 'R', nume_ddlz=nume_dof, &
                nb_equa_outz=nb_equa)
    call vtcreb(sdcont_rea2, 'V', 'R', nume_ddlz=nume_dof)
!
! - Parameters for "PENALISATION" and "LAGRANGIEN"
!
    sdcont_tacfin = ds_contact%sdcont_solv(1:14)//'.TACFIN'
    call wkvect(sdcont_tacfin, 'V V R', nt_cont_poin*ztacf, jv_sdcont_tacfin)
!
! - Tangents
!
    sdcont_tangco = ds_contact%sdcont_solv(1:14)//'.TANGCO'
    call wkvect(sdcont_tangco, 'V V R', 6*nt_cont_poin, jv_sdcont_tangco)
!
! - Datastructures for gaps
!
    call cfcrje(ds_contact)
!
! - Datastructures for linear relations
!
    call cfcrli(mesh, nume_dof, ds_contact)
!
! - Lagrange multipliers
!
    sdcont_mu = ds_contact%sdcont_solv(1:14)//'.MU'
    call wkvect(sdcont_mu, 'V V R', 4*nt_cont_poin, jv_sdcont_mu)
!
! - Pseudo-penalization value for 3D lagangrian
!
    sdcont_copo = ds_contact%sdcont_solv(1:14)//'.COPO'
    call wkvect(sdcont_copo, 'V V R', 1, vr=v_sdcont_copo)
    v_sdcont_copo(1) = r8vide()
!
! - Contact forces
!
    sdcont_atmu = ds_contact%sdcont_solv(1:14)//'.ATMU'
    call wkvect(sdcont_atmu, 'V V R', nb_equa, jv_sdcont_atmu)
!
! - Friction forces
!
    if (l_frot) then
        sdcont_afmu = ds_contact%sdcont_solv(1:14)//'.AFMU'
        call wkvect(sdcont_afmu, 'V V R', nb_equa, jv_sdcont_afmu)
    end if
!
! - Contact forces for contact penalization methods
!
    if (l_pena_cont .and. (.not. l_frot)) then
        sdcont_afmu = ds_contact%sdcont_solv(1:14)//'.AFMU'
        call wkvect(sdcont_afmu, 'V V R', nb_equa, jv_sdcont_afmu)
    end if
!
! - Solution increment without contact
!
    sdcont_del0 = ds_contact%sdcont_solv(1:14)//'.DEL0'
    call vtcreb(sdcont_del0, 'V', 'R', nume_ddlz=nume_dof)
!
! - Solution increment for contact iteration
!
    sdcont_ddel = ds_contact%sdcont_solv(1:14)//'.DDEL'
    call vtcreb(sdcont_ddel, 'V', 'R', nume_ddlz=nume_dof)
!
! - Solution increment with contact
!
    sdcont_delc = ds_contact%sdcont_solv(1:14)//'.DELC'
    call vtcreb(sdcont_delc, 'V', 'R', nume_ddlz=nume_dof)
!
! - Solution increment since beginning of time step without contact
!
    if (l_frot) then
        sdcont_dep0 = ds_contact%sdcont_solv(1:14)//'.DEP0'
        call vtcreb(sdcont_dep0, 'V', 'R', nume_ddlz=nume_dof)
    end if
!
! - Solution increment since beginning of time step with contact
!
    if (l_frot) then
        sdcont_depc = ds_contact%sdcont_solv(1:14)//'.DEPC'
        call vtcreb(sdcont_depc, 'V', 'R', nume_ddlz=nume_dof)
    end if
!
! - Void kinematic load
!
    sdcont_cin0 = ds_contact%sdcont_solv(1:14)//'.CIN0'
    call vtcreb(sdcont_cin0, 'V', 'R', nume_ddlz=nume_dof)
!
! - Fields for GCP
!
    if (l_gcp) then
        sdcont_sgdm = ds_contact%sdcont_solv(1:14)//'.SGDM'
        sdcont_sgdp = ds_contact%sdcont_solv(1:14)//'.SGDP'
        sdcont_dire = ds_contact%sdcont_solv(1:14)//'.DIRE'
        sdcont_sgpm = ds_contact%sdcont_solv(1:14)//'.SGPM'
        sdcont_sgpp = ds_contact%sdcont_solv(1:14)//'.SGPP'
        sdcont_mum = ds_contact%sdcont_solv(1:14)//'.MUM'
        sdcont_secm = ds_contact%sdcont_solv(1:14)//'.SECM'
        call vtcreb(sdcont_secm, 'V', 'R', nume_ddlz=nume_dof)
        call wkvect(sdcont_sgdm, 'V V R', nt_cont_poin, jv_sdcont_sgdm)
        call wkvect(sdcont_sgdp, 'V V R', nt_cont_poin, jv_sdcont_sgdp)
        call wkvect(sdcont_sgpm, 'V V R', nt_cont_poin, jv_sdcont_sgpm)
        call wkvect(sdcont_sgpp, 'V V R', nt_cont_poin, jv_sdcont_sgpp)
        call wkvect(sdcont_dire, 'V V R', nt_cont_poin, jv_sdcont_dire)
        call wkvect(sdcont_mum, 'V V R', nt_cont_poin, jv_sdcont_mum)
        if (l_pre_cond) then
            sdcont_pcrs = ds_contact%sdcont_solv(1:14)//'.PCRS'
            call wkvect(sdcont_pcrs, 'V V R', nt_cont_poin, jv_sdcont_pcrs)
            sdcont_pcdr = ds_contact%sdcont_solv(1:14)//'.PCDR'
            call wkvect(sdcont_pcdr, 'V V R', nt_cont_poin, jv_sdcont_pcdr)
            sdcont_pcuu = ds_contact%sdcont_solv(1:14)//'.PCUU'
            call wkvect(sdcont_pcuu, 'V V R', nb_equa, jv_sdcont_pcuu)
        end if
    end if
!
! - State saving for GCP
!
    if (l_gcp) then
        sdcont_svm0 = ds_contact%sdcont_solv(1:14)//'.SVM0'
        call wkvect(sdcont_svm0, 'V V R', nb_cont_node, jv_sdcont_svm0)
        sdcont_svmu = ds_contact%sdcont_solv(1:14)//'.SVMU'
        call wkvect(sdcont_svmu, 'V V R', nb_cont_node, jv_sdcont_svmu)
    end if
!
! - Contact matrix parameters: sizes
!
    if (.not. l_gcp) then
        nbliai = nt_cont_poin
        nbenat = 0
        nbcm1a = 0
        nbfro1 = 0
        nbfro2 = 0
        if (l_pena_cont) then
            nbenat = nbliai
        else if ((.not. l_frot) .or. l_pena_frot) then
            nbcm1a = nbliai
        else
            nbcm1a = model_ndim*nbliai
        end if
        nbfro1 = (model_ndim-1)*nbliai
        nbfro2 = nbliai
    end if
!
! - Contact matrix parameters: sizes
!
    if (.not. l_gcp) then
        if (l_pena_cont) then
            sdcont_enat = ds_contact%sdcont_solv(1:14)//'.ENAT'
            call jecrec(sdcont_enat, 'V V R', 'NU', 'DISPERSE', 'CONSTANT', nbenat)
            call jeecra(sdcont_enat, 'LONMAX', ival=30)
            do ii = 1, nbenat
                call jecroc(jexnum(sdcont_enat, ii))
            end do
        else
            sdcont_cm1a = ds_contact%sdcont_solv(1:14)//'.CM1A'
            call jecrec(sdcont_cm1a, 'V V R', 'NU', 'DISPERSE', 'CONSTANT', nbcm1a)
            call jeecra(sdcont_cm1a, 'LONMAX', ival=nb_equa)
            do ii = 1, nbcm1a
                call jecroc(jexnum(sdcont_cm1a, ii))
            end do
        end if
        if (l_frot) then
            sdcont_fro1 = ds_contact%sdcont_solv(1:14)//'.FRO1'
            sdcont_fro2 = ds_contact%sdcont_solv(1:14)//'.FRO2'
            call jecrec(sdcont_fro1, 'V V R', 'NU', 'DISPERSE', 'CONSTANT', nbfro1)
            call jecrec(sdcont_fro2, 'V V R', 'NU', 'DISPERSE', 'CONSTANT', nbfro2)
            call jeecra(sdcont_fro1, 'LONMAX', ival=30)
            call jeecra(sdcont_fro2, 'LONMAX', ival=30)
            do ii = 1, nbfro1
                call jecroc(jexnum(sdcont_fro1, ii))
            end do
            do ii = 1, nbfro2
                call jecroc(jexnum(sdcont_fro2, ii))
            end do
        end if
    end if
!
! - Matrix A.C-1.AT
!
    if (l_matr_cont) then
        call cfcrma(nbcm1a, mesh, ds_contact%sdcont_solv)
    end if
!
! - Forces to solve
!
    call vtcreb(ds_contact%cnctdc, 'V', 'R', nume_ddlz=nume_dof)
    ds_contact%l_cnctdc = ASTER_TRUE
    if (l_frot .or. l_pena_cont) then
        call vtcreb(ds_contact%cnctdf, 'V', 'R', nume_ddlz=nume_dof)
        ds_contact%l_cnctdf = ASTER_TRUE
    end if
!
    call jedema()
end subroutine
