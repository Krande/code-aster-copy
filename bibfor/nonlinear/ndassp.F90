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
subroutine ndassp(listFuncActi, ds_contact, ds_system, &
                  sddyna, nlDynaDamping, &
                  hval_veasse, cndonn)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisl.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndasva.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmasdi.h"
#include "asterfort/nmasfi.h"
#include "asterfort/nmasva.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nonlinDSVectCombInit.h"
#include "asterfort/nonlinDSVectCombCompute.h"
#include "asterfort/nonlinDSVectCombAddAny.h"
#include "asterfort/nonlinDSVectCombAddHat.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19), intent(in) :: hval_veasse(*)
    character(len=19), intent(in) :: cndonn
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Evaluate second member for prediction (dynamic)
!
! --------------------------------------------------------------------------------------------------
!
! In  listFuncActi     : list of active functionnalities
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  ds_system        : datastructure for non-linear system management
! In  cndonn           : name of nodal field for "given" forces
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: cnffdo, cndfdo, cnfvdo, cnvady
    aster_logical :: lSuperElement
    real(kind=8) :: coeequ
    type(NL_DS_VectComb) :: ds_vectcomb
    aster_logical :: l_unil_pena
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_20')
    end if

! - Initializations
    call nonlinDSVectCombInit(ds_vectcomb)
    cnffdo = '&&CNCHAR.FFDO'
    cndfdo = '&&CNCHAR.DFDO'
    cnfvdo = '&&CNCHAR.FVDO'
    cnvady = '&&CNCHAR.FVDY'
    lSuperElement = isfonc(listFuncActi, 'MACR_ELEM_STAT')

! - Coefficient for multi-step scheme
    coeequ = ndynre(sddyna, 'COEF_MPAS_EQUI_COUR')

! - Get dead Neumann loads and multi-step dynamic schemes forces
    call nmasfi(listFuncActi, hval_veasse, cnffdo, sddyna)
    call nonlinDSVectCombAddAny(cnffdo, +1.d0, ds_vectcomb)

! - Get Dirichlet loads
    call nmasdi(listFuncActi, hval_veasse, cndfdo)
    call nonlinDSVectCombAddAny(cndfdo, +1.d0, ds_vectcomb)

! - Get undead Neumann loads and multi-step dynamic schemes forces
    call nmasva(listFuncActi, hval_veasse, cnfvdo, sddyna)
    call nonlinDSVectCombAddAny(cnfvdo, +1.d0, ds_vectcomb)

! - Get undead Neumann loads for dynamic
    call ndasva(sddyna, nlDynaDamping, hval_veasse, cnvady)
    call nonlinDSVectCombAddAny(cnvady, coeequ, ds_vectcomb)

! - Add DISCRETE contact force
    if (ds_contact%l_cnctdf) then
        call nonlinDSVectCombAddAny(ds_contact%cnctdf, -1.d0, ds_vectcomb)
    end if

! - Add LIAISON_UNIL penalized force
    if (ds_contact%l_cnunil) then
        l_unil_pena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
        if (l_unil_pena) then
            call nonlinDSVectCombAddAny(ds_contact%cnunil, -1.d0, ds_vectcomb)
        end if
    end if

! - Add CONTINUE contact forces
    if (ds_contact%l_cneltc) then
        call nonlinDSVectCombAddAny(ds_contact%cneltc, -1.d0, ds_vectcomb)
    end if
    if (ds_contact%l_cneltf) then
        call nonlinDSVectCombAddAny(ds_contact%cneltf, -1.d0, ds_vectcomb)
    end if

! - Force from sub-structuring
    if (lSuperElement) then
        call nonlinDSVectCombAddHat(hval_veasse, 'CNSSTR', -1.d0, ds_vectcomb)
    end if
!
! - Get Dirichlet boundary conditions - B.U
!
    call nonlinDSVectCombAddHat(hval_veasse, 'CNBUDI', -1.d0, ds_vectcomb)
!
! - Get force for Dirichlet boundary conditions (dualized) - BT.LAMBDA
!
    call nonlinDSVectCombAddHat(hval_veasse, 'CNDIRI', -1.d0, ds_vectcomb)
!
! - Add internal forces to second member
!
    call nonlinDSVectCombAddAny(ds_system%cnfint, -1.d0, ds_vectcomb)
!
! - Combination
!
    call nonlinDSVectCombCompute(ds_vectcomb, cndonn)
!
! - Debug
!
    if (niv .ge. 2) then
        call nmdebg('VECT', cndonn, 6)
    end if

!
end subroutine
