! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine carc_read(behaviourPrepCrit, model_)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getexm.h"
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lcsymm.h"
#include "asterc/lctest.h"
#include "asterfort/assert.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/dismoi.h"
#include "asterfort/exicp.h"
#include "asterfort/getBehaviourAlgo.h"
#include "asterfort/getBehaviourPara.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/getExternalBehaviourPntr.h"
#include "asterfort/getExternalStateVariable.h"
#include "asterfort/getExternalStrainModel.h"
#include "asterfort/getHHOPara.h"
#include "asterfort/getTHMPara.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
type(Behaviour_PrepCrit), intent(inout) :: behaviourPrepCrit
character(len=8), intent(in), optional :: model_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Read from command file
!
! --------------------------------------------------------------------------------------------------
!
! IO  behaviourPrepCrit: datastructure to prepare parameters for constitutive laws
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=16) :: answer
    integer :: iFactorKeyword, iret, nbFactorKeyword
    character(len=16) :: type_matr_tang, method, post_iter, post_incr
    real(kind=8) :: parm_theta, vale_pert_rela
    real(kind=8) :: resi_deborst_max
    real(kind=8) :: resi_radi_rela
    real(kind=8) :: parm_theta_thm, parm_alpha_thm
    integer :: type_matr_t, iter_inte_pas, iter_deborst_max
    integer :: ipostiter, ipostincr, iveriborne
    character(len=8) :: mesh
    character(len=16) :: rela_code_py, defo_code_py, meca_code_py
    character(len=16) :: veri_borne
    character(len=16) :: kit_comp(4)
    character(len=16) :: defo_comp,  rela_comp
    character(len=16) :: thmc_comp, hydr_comp, ther_comp, meca_comp
    aster_logical :: l_kit_thm, l_kit_ddi, l_exist_thm
    aster_logical :: l_kit
    aster_logical :: plane_stress, l_mfront_proto, l_mfront_offi
    character(len=24), parameter :: list_elem_affe = '&&CARCREAD.LIST'
    aster_logical :: l_affe_all, l_matr_unsymm, l_comp_external
    integer :: nb_elem_affe
    integer :: cptr_nbvarext, cptr_namevarext, cptr_fct_ldc
    integer :: cptr_nameprop, cptr_nbprop
    integer :: variExteCode(2), exte_strain
    character(len=16) :: texte(3)
    integer, pointer :: modelCell(:) => null()
    character(len=16) :: algo_inte
    real(kind=8) :: algo_inte_r, iter_inte_maxi, resi_inte_rela
    type(Behaviour_ParaExte) :: paraExte
!
! --------------------------------------------------------------------------------------------------
!
    iFactorKeyword=0
    iret=0
    type_matr_tang=' '
    method=' '
    post_iter=' '
    post_incr=' '
    parm_theta=0.d0
    vale_pert_rela=0.d0
    resi_deborst_max=0.d0
    resi_radi_rela=0.d0
    parm_theta_thm=0.d0
    parm_alpha_thm=0.d0
    type_matr_t=0
    iter_inte_pas=0
    iter_deborst_max=0
    ipostiter=0
    ipostincr=0
    iveriborne=0
    mesh = ' '
    veri_borne=' '
    kit_comp(1:4) = (/'VIDE','VIDE','VIDE','VIDE'/)
    defo_comp=' '
    rela_comp=' '
    thmc_comp=' '
    hydr_comp=' '
    ther_comp=' '
    meca_comp=' '
    l_kit_thm=ASTER_FALSE
    l_kit_ddi = ASTER_FALSE
    l_exist_thm = ASTER_FALSE
    l_kit = ASTER_FALSE
    texte(:)=(/ ' ',' ',' '/)
    nbFactorKeyword = behaviourPrepCrit%nb_comp
    mesh           = ' '
    l_exist_thm    = ASTER_FALSE

! - Pointer to the list of cells in model
    mesh = ' '
    if ( present(model_) ) then
        call jeveuo(model_//'.MAILLE', 'L', vi = modelCell)
        call dismoi('NOM_MAILLA', model_, 'MODELE', repk = mesh)
    endif

! - Read informations
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get RELATION from command file
        call compGetRelation(factorKeyword, iFactorKeyword, rela_comp)

! ----- Get DEFORMATION from command file
        call getvtx(factorKeyword, 'DEFORMATION', iocc = iFactorKeyword, scal = defo_comp)

! ----- Detection of specific cases
        call comp_meca_l(rela_comp, 'KIT'         , l_kit)
        call comp_meca_l(rela_comp, 'KIT_THM'     , l_kit_thm)
        call comp_meca_l(rela_comp, 'KIT_DDI'     , l_kit_ddi)
        call comp_meca_l(rela_comp, 'MFRONT_OFFI' , l_mfront_offi)
        call comp_meca_l(rela_comp, 'MFRONT_PROTO', l_mfront_proto)
        if (l_kit_thm) then
            l_exist_thm = ASTER_TRUE
        endif

! ----- For KIT
        if (l_kit) then
            call comp_meca_rkit(factorKeyword, iFactorKeyword, rela_comp, kit_comp)
        endif

! ----- Get mechanical part of behaviour
        call compGetMecaPart(rela_comp, kit_comp, meca_comp)

! ----- Coding comportment (Python)
        call lccree(1, rela_comp, rela_code_py)
        call lccree(1, defo_comp, defo_code_py)
        call lccree(1, meca_comp, meca_code_py)

! ----- Symmetric or not ?
        l_matr_unsymm = ASTER_FALSE
        call lcsymm(rela_code_py, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        call lcsymm(meca_code_py, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        call lcsymm(defo_code_py, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        if (l_mfront_proto) then
            call getvtx(factorKeyword, 'SYME_MATR_TANG', iocc = iFactorKeyword,&
                        scal = answer, nbret = iret)
            if (iret .ne. 0) then
                l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'NON'
            endif
        endif

! ----- Get ITER_INTE_PAS
        call getvis(factorKeyword, 'ITER_INTE_PAS', iocc = iFactorKeyword,&
                    scal=iter_inte_pas, nbret=iret)
        if (iret .eq. 0) then
            iter_inte_pas = 0
        endif

! ----- Get ITER_CPLAN_MAXI/RESI_CPLAN_MAXI/RESI_CPLAN_RELA (Deborst method)
        resi_deborst_max = 1.d-6
        iter_deborst_max = 1
        call getvis(factorKeyword, 'ITER_CPLAN_MAXI', iocc = iFactorKeyword,&
                    scal = iter_deborst_max)
        call getvr8(factorKeyword, 'RESI_CPLAN_MAXI', iocc = iFactorKeyword,&
                    scal = resi_deborst_max, nbret = iret)
        if (iret .ne. 0) then
            resi_deborst_max = -resi_deborst_max
        else
            call getvr8(factorKeyword, 'RESI_CPLAN_RELA', iocc = iFactorKeyword,&
                        scal = resi_deborst_max)
        endif

! ----- Get TYPE_MATR_TANG/VALE_PERT_RELA
        vale_pert_rela = 0.d0
        type_matr_t    = 0
        type_matr_tang = ' '
        call getvtx(factorKeyword, 'TYPE_MATR_TANG', iocc = iFactorKeyword,&
                    scal = type_matr_tang, nbret = iret)
        if (iret .eq. 0) then
            type_matr_t = 0
        else
            if (type_matr_tang .eq. 'PERTURBATION') then
                type_matr_t = 1
                call getvr8(factorKeyword, 'VALE_PERT_RELA', iocc = iFactorKeyword,&
                            scal = vale_pert_rela)
            else if (type_matr_tang .eq. 'VERIFICATION') then
                type_matr_t = 2
                call getvr8(factorKeyword, 'VALE_PERT_RELA', iocc = iFactorKeyword,&
                            scal = vale_pert_rela)
            else
                ASSERT(.false.)
            endif
            call lctest(rela_code_py, 'TYPE_MATR_TANG', type_matr_tang, iret)
            if (iret .eq. 0) then
                texte(1) = type_matr_tang
                texte(2) = rela_comp
                call utmess('F', 'COMPOR1_46', nk = 2, valk = texte)
            endif
        endif

! ----- Get PARM_THETA (for viscous laws)
        parm_theta = 1.d0
        call getvr8(factorKeyword, 'PARM_THETA', iocc = iFactorKeyword, scal = parm_theta)

! ----- Get RESI_RADI_RELA
        if (type_matr_t .eq. 0) then
            call getvr8(factorKeyword, 'RESI_RADI_RELA', iocc = iFactorKeyword,&
                        scal = resi_radi_rela, nbret = iret)
            if (iret .eq. 0) then
                resi_radi_rela = -10.d0
            endif
        endif

! ----- Get POST_ITER
        ipostiter = 0
        if (getexm(factorKeyword,'POST_ITER') .eq. 1) then
            post_iter = ' '
            if (type_matr_t .eq. 0) then
                call getvtx(factorKeyword, 'POST_ITER', iocc = iFactorKeyword,&
                            scal = post_iter, nbret = iret)
                if (iret .eq. 1) then
                    if (post_iter .eq. 'CRIT_RUPT') then
                        ipostiter = 1
                    endif
                endif
            endif
        endif

! ----- Get POST_INCR
        ipostincr = 0
        if (getexm(factorKeyword,'POST_INCR') .eq. 1) then
            post_incr = ' '
            call getvtx(factorKeyword, 'POST_INCR', iocc = iFactorKeyword,&
                        scal = post_incr, nbret = iret)
            if (iret .eq. 1) then
               if (post_incr .eq. 'REST_ECRO') then
                    ipostincr = 1
               endif
            endif
        endif

! ----- Get VERI_BORNE
        iveriborne = 0
        if (getexm(factorKeyword,'VERI_BORNE') .eq. 1) then
            call getvtx(factorKeyword, 'VERI_BORNE', iocc = iFactorKeyword,&
                         scal = veri_borne, nbret = iret )
            if (iret .eq. 0) then
                iveriborne = 2
            else
                if ( veri_borne .eq. 'ARRET' ) then
                    iveriborne = 2
                elseif ( veri_borne .eq. 'MESSAGE' ) then
                    iveriborne = 1
                else
                    iveriborne = 0
                endif
            endif
        endif

! ----- Get parameters for external programs (MFRONT/UMAT)
        call getExternalBehaviourPara(mesh, modelCell, rela_comp, kit_comp,&
                                      l_comp_external, paraExte,&
                                      factorKeyword, iFactorKeyword)

! ----- Get list of elements where comportment is defined
        plane_stress = ASTER_FALSE
        if (present(model_)) then
            call comp_read_mesh(mesh, factorKeyword, iFactorKeyword,&
                                list_elem_affe, l_affe_all , nb_elem_affe)
            plane_stress = exicp(model_, l_affe_all, list_elem_affe, nb_elem_affe)
        endif

! ----- Get ALGO_INTE
        algo_inte_r = 0.d0
        call getBehaviourAlgo(plane_stress, rela_comp   ,&
                              rela_code_py, meca_code_py,&
                              factorKeyword , iFactorKeyword      ,&
                              algo_inte   , algo_inte_r)

! ----- Get function pointers for external programs (MFRONT/UMAT)
        cptr_fct_ldc    = 0
        cptr_nbvarext   = 0
        cptr_namevarext = 0
        cptr_nbprop     = 0
        cptr_nameprop   = 0
        if (l_comp_external) then
            call getExternalBehaviourPntr(paraExte,&
                                          cptr_fct_ldc ,&
                                          cptr_nbvarext, cptr_namevarext,&
                                          cptr_nbprop  , cptr_nameprop)
        endif

! ----- Get RESI_INTE_RELA/ITER_INTE_MAXI
        resi_inte_rela = 0.d0
        iter_inte_maxi = 0
        call getBehaviourPara(l_mfront_offi , l_mfront_proto, l_kit_thm,&
                              factorKeyword   , iFactorKeyword        , algo_inte,&
                              iter_inte_maxi, resi_inte_rela)

! ----- Get external state variables
        call getExternalStateVariable(rela_comp, rela_code_py ,&
                                      l_mfront_offi, l_mfront_proto ,&
                                      cptr_nbvarext, cptr_namevarext,&
                                      variExteCode)

! ----- Get model of strains for external programs (MFRONT)
        exte_strain = 0
        call getExternalStrainModel(l_mfront_offi, l_mfront_proto,&
                                    paraExte,&
                                    defo_comp, exte_strain)

! ----- Discard
        call lcdiscard(meca_code_py)
        call lcdiscard(rela_code_py)
        call lcdiscard(defo_code_py)

! ----- Save parameters
        behaviourPrepCrit%v_crit(iFactorKeyword)%rela_comp        = rela_comp
        behaviourPrepCrit%v_crit(iFactorKeyword)%meca_comp        = meca_comp
        behaviourPrepCrit%v_crit(iFactorKeyword)%type_matr_t      = type_matr_t
        behaviourPrepCrit%v_crit(iFactorKeyword)%parm_theta       = parm_theta
        behaviourPrepCrit%v_crit(iFactorKeyword)%iter_inte_pas    = iter_inte_pas
        behaviourPrepCrit%v_crit(iFactorKeyword)%vale_pert_rela   = vale_pert_rela
        behaviourPrepCrit%v_crit(iFactorKeyword)%resi_deborst_max = resi_deborst_max
        behaviourPrepCrit%v_crit(iFactorKeyword)%iter_deborst_max = iter_deborst_max
        behaviourPrepCrit%v_crit(iFactorKeyword)%resi_radi_rela   = resi_radi_rela
        behaviourPrepCrit%v_crit(iFactorKeyword)%ipostiter        = ipostiter
        behaviourPrepCrit%v_crit(iFactorKeyword)%ipostincr        = ipostincr
        behaviourPrepCrit%v_crit(iFactorKeyword)%iveriborne       = iveriborne
        behaviourPrepCrit%v_crit(iFactorKeyword)%l_matr_unsymm    = l_matr_unsymm
        behaviourPrepCrit%v_crit(iFactorKeyword)%algo_inte_r      = algo_inte_r
        behaviourPrepCrit%v_crit(iFactorKeyword)%resi_inte_rela   = resi_inte_rela
        behaviourPrepCrit%v_crit(iFactorKeyword)%iter_inte_maxi   = iter_inte_maxi
        behaviourPrepCrit%v_crit(iFactorKeyword)%cptr_fct_ldc     = cptr_fct_ldc
        behaviourPrepCrit%v_crit(iFactorKeyword)%cptr_nbvarext    = cptr_nbvarext
        behaviourPrepCrit%v_crit(iFactorKeyword)%cptr_namevarext  = cptr_namevarext
        behaviourPrepCrit%v_crit(iFactorKeyword)%cptr_nbprop      = cptr_nbprop
        behaviourPrepCrit%v_crit(iFactorKeyword)%cptr_nameprop    = cptr_nameprop
        behaviourPrepCrit%v_crit(iFactorKeyword)%jvariext1        = variExteCode(1)
        behaviourPrepCrit%v_crit(iFactorKeyword)%jvariext2        = variExteCode(2)
        behaviourPrepCrit%v_crit(iFactorKeyword)%exte_strain      = exte_strain
        behaviourPrepCrit%v_crit(iFactorKeyword)%paraExte         = paraExte
    end do

! - Get SCHEMA_THM parameters
    if (l_exist_thm) then
        call getTHMPara(behaviourPrepCrit)
    endif

! - Get HHO parameters
    call getHHOPara(behaviourPrepCrit)
!
end subroutine
