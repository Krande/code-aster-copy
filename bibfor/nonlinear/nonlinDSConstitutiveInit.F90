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
subroutine nonlinDSConstitutiveInit(modelZ, caraElemZ, ds_constitutive, verbose_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/nmdoco.h"
#include "asterfort/nmcpqu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/comp_info.h"
#include "asterfort/cesvar.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ
    type(NL_DS_Constitutive), intent(inout) :: ds_constitutive
    aster_logical, intent(in), optional :: verbose_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Constitutive laws
!
! Initializations for constitutive laws management datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! IO  ds_constitutive  : datastructure for constitutive laws management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    aster_logical :: lLinear, lDisCtc, verbose
    integer(kind=8) :: nb_affe, i_affe
    character(len=16), pointer :: v_compor_vale(:) => null()
    integer(kind=8), pointer :: v_compor_desc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_2')
    end if
!
    if (present(verbose_)) then
        verbose = verbose_
    else
        verbose = ASTER_TRUE
    end if
!
! - Construct CHAM_ELEM_S
!
    call nmdoco(modelZ, caraElemZ, ds_constitutive%compor)
!
! - Some functionnalities
!
    call nmcpqu(ds_constitutive%compor, 'C_PLAN', 'DEBORST', ds_constitutive%l_deborst)
    call nmcpqu(ds_constitutive%compor, 'RELCOM', 'DIS_CHOC', ds_constitutive%l_dis_choc)
    call nmcpqu(ds_constitutive%compor, 'POSTINCR', 'REST_ECRO', ds_constitutive%lAnnealing)
    call nmcpqu(ds_constitutive%compor, 'RELCOM', 'JOINT_MECA_FROT', ds_constitutive%lJoiFrot)
    call nmcpqu(ds_constitutive%compor, 'RELCOM', 'JOINT_MECA_RUPT', ds_constitutive%lJoiRupt)
    call nmcpqu(ds_constitutive%compor, 'RELCOM', 'JOINT_MECA_ENDO', ds_constitutive%lJoiEndo)

!
! - Look if geometric matrix is included in global tangent matrix
!
    call jeveuo(ds_constitutive%compor(1:19)//'.VALE', 'L', vk16=v_compor_vale)
    call jeveuo(ds_constitutive%compor(1:19)//'.DESC', 'L', vi=v_compor_desc)
    nb_affe = v_compor_desc(3)
    do i_affe = 1, nb_affe
        if ((v_compor_vale(DEFO+COMPOR_SIZE*(i_affe-1)) .eq. 'GROT_GDEP') .or. &
            (v_compor_vale(DEFO+COMPOR_SIZE*(i_affe-1)) .eq. 'SIMO_MIEHE') .or. &
            (v_compor_vale(DEFO+COMPOR_SIZE*(i_affe-1)) .eq. 'GREEN_LAGRANGE') .or. &
            (v_compor_vale(DEFO+COMPOR_SIZE*(i_affe-1)) .eq. 'GDEF_LOG')) then
            ds_constitutive%l_matr_geom = ASTER_TRUE
        end if
    end do
!
! - Check if linear
!
    lLinear = ASTER_TRUE
    do i_affe = 1, nb_affe
        if ((v_compor_vale(RELA_NAME+COMPOR_SIZE*(i_affe-1)) .ne. 'ELAS')) then
            lLinear = ASTER_FALSE
        end if
        if ((v_compor_vale(DEFO+COMPOR_SIZE*(i_affe-1)) .ne. 'PETIT')) then
            lLinear = ASTER_FALSE
        end if
    end do
    ds_constitutive%lLinear = lLinear
!
! - Check if DIS_*
!
    lDisCtc = ASTER_FALSE
    do i_affe = 1, nb_affe
        if ((v_compor_vale(RELA_NAME+COMPOR_SIZE*(i_affe-1)) .eq. 'DIS_CHOC') .or. &
            (v_compor_vale(RELA_NAME+COMPOR_SIZE*(i_affe-1)) .eq. 'DIS_CONTACT') .or. &
            (v_compor_vale(RELA_NAME+COMPOR_SIZE*(i_affe-1)) .eq. 'CHOC_ENDO_PENA') .or. &
            (v_compor_vale(RELA_NAME+COMPOR_SIZE*(i_affe-1)) .eq. 'CHOC_ENDO')) then
            lDisCtc = ASTER_TRUE
        end if
    end do
    ds_constitutive%lDisCtc = lDisCtc
!
! - Print informations about COMPORTEMENT keyword
!
    if (verbose) call comp_info(modelZ, ds_constitutive%compor)
!
! - Preallocation of output stress field for prediction
!
    call cesvar(caraElemZ, ds_constitutive%compor, modelZ(1:8)//'.MODELE', &
                ds_constitutive%sigm_pred)

!
end subroutine
