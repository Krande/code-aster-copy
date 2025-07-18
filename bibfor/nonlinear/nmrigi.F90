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
subroutine nmrigi(modelz, cara_elem, &
                  ds_material, ds_constitutive, &
                  list_func_acti, iter_newt, sddyna, ds_measure, ds_system, &
                  hval_incr, hval_algo, &
                  optioz, ldccvg)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/isfonc.h"
#include "asterfort/merimo.h"
#include "asterfort/nmdep0.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/utmess.h"
#include "asterfort/infdbg.h"
!
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: cara_elem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    integer(kind=8), intent(in) :: list_func_acti(*)
    integer(kind=8), intent(in) :: iter_newt
    character(len=19), intent(in) :: sddyna
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=*), intent(in) :: optioz
    integer(kind=8), intent(out) :: ldccvg
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute elementaries for internal forces/rigidity matrix
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  OPTRIG : OPTION DE CALCUL POUR MERIMO
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  SDDYNA : SD POUR LA DYNAMIQUE
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! IN  ITERAT : NUMERO D'ITERATION
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 2 : ERREUR SUR LA NON VERIF. DE CRITERES PHYSIQUES
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=1) :: base
    character(len=24) :: model
    character(len=16) :: optrig
    aster_logical :: lendo, l_xfem, l_macr_elem
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECA_NON_LINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_29')
    end if
!
! - Initializations
!
    base = 'V'
    optrig = optioz
    model = modelz
    ldccvg = 0
!
! - Active functionnalities
!
    l_xfem = isfonc(list_func_acti, 'XFEM')
    l_macr_elem = isfonc(list_func_acti, 'MACR_ELEM_STAT')
    lendo = isfonc(list_func_acti, 'ENDO_NO')
!
! --- INCREMENT DE DEPLACEMENT NUL EN PREDICTION
!
    if (.not. lendo) then
        if (optrig(1:9) .eq. 'RIGI_MECA') then
            call nmdep0('ON ', hval_algo)
        end if
    end if
!
! - Init timer
!
    call nmtime(ds_measure, 'Init', 'Integrate')
    call nmtime(ds_measure, 'Launch', 'Integrate')
!
! - Computation
!
    call merimo(base, &
                l_xfem, l_macr_elem, &
                model, cara_elem, iter_newt+1, &
                ds_constitutive, ds_material, ds_system, &
                hval_incr, hval_algo, &
                optrig, ldccvg, sddyna)
!
! - End timer
!
    call nmtime(ds_measure, 'Stop', 'Integrate')
    call nmrinc(ds_measure, 'Integrate')
!
! --- REMISE INCREMENT DE DEPLACEMENT
!
    if (optrig(1:9) .eq. 'RIGI_MECA') then
        call nmdep0('OFF', hval_algo)
    end if
!
end subroutine
