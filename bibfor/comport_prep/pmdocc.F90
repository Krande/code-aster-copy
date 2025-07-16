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
subroutine pmdocc(compor, nbVari, type_comp, mult_comp)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/comp_meca_cvar.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_info.h"
#include "asterfort/comp_meca_pvar.h"
#include "asterfort/comp_meca_read.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/imvari.h"
#include "asterfort/utmess.h"
#include "asterfort/setBehaviourTypeValue.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(out) :: compor(COMPOR_SIZE)
    integer(kind=8), intent(out) :: nbVari
    character(len=16), intent(out) :: type_comp, mult_comp
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviours (mechanics/SIMU_POINT_MAT)
!
! Get list of parameters for constitutive law
!
! --------------------------------------------------------------------------------------------------
!
! Out compor           : list of parameters for constitutive law
! Out nbVari           : number of internal variables
! Out type_comp        : type of comportment (INCR/ELAS)
! Out mult_comp        : multi-comportment (DEFI_COMPOR for PMF)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: comporInfo = '&&PMDORC.LIST_VARI'
    integer(kind=8) :: nbocc1, nbocc2, nbocc3
    character(len=16) :: rela_comp
    aster_logical :: l_etat_init, l_kit_thm
    type(BehaviourPrep_MapCompor) :: prepMapCompor
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    type_comp = ' '
    mult_comp = ' '
    compor(1:COMPOR_SIZE) = 'VIDE'

! - Initial state
    call getfac('SIGM_INIT', nbocc1)
    call getfac('EPSI_INIT', nbocc2)
    call getfac('VARI_INIT', nbocc3)
    l_etat_init = (nbocc1+nbocc2+nbocc3) > 0

! - Create datastructure to prepare comportement
    call comp_meca_info(prepMapCompor)
    if (prepMapCompor%nb_comp .eq. 0) then
        call utmess('F', 'COMPOR4_63')
    end if

! - Read informations from command file
    call comp_meca_read(l_etat_init, prepMapCompor)

! - Count internal variables
    call comp_meca_cvar(prepMapCompor)

! - Some properties
    nbVari = prepMapCompor%prepPara(1)%nbVari
    rela_comp = prepMapCompor%prepPara(1)%rela_comp
    type_comp = prepMapCompor%prepPara(1)%type_comp
    mult_comp = prepMapCompor%prepPara(1)%mult_comp

! - Detection of specific cases
    call comp_meca_l(rela_comp, 'KIT_THM', l_kit_thm)
    if (l_kit_thm) then
        call utmess('F', 'COMPOR2_7')
    end if

! - Save informations in the field <COMPOR>
    call setBehaviourTypeValue(prepMapCompor, comporList_=compor(1:COMPOR_SIZE))

! - Prepare informations about internal variables
    call comp_meca_pvar(comporList_=compor, comporInfo=comporInfo)

! - Print informations about internal variables
    call imvari(comporInfo)

! - Cleaning
    deallocate (prepMapCompor%prepPara)
    deallocate (prepMapCompor%prepExte)
!
end subroutine
