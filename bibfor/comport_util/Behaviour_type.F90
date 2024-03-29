! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
module Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour: types for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
! - Behaviour - Integration - External state variables
!
    type Behaviour_ESVA
! ----- Flags
        aster_logical :: l_anel = ASTER_FALSE
        aster_logical :: l_temp = ASTER_FALSE
        aster_logical :: l_sech = ASTER_FALSE
        aster_logical :: l_hydr = ASTER_FALSE
        aster_logical :: l_hygr = ASTER_FALSE
        aster_logical :: l_ptot = ASTER_FALSE
! ----- Anelastic strains
        real(kind=8)  :: anel_prev(6) = 0.d0, anel_curr(6) = 0.d0, anel_incr(6) = 0.d0
! ----- Temperature
        real(kind=8)  :: temp_prev = 0.d0, temp_curr = 0.d0, temp_incr = 0.d0, temp_refe = 0.d0
! ----- Drying
        real(kind=8)  :: sech_prev = 0.d0, sech_curr = 0.d0, sech_incr = 0.d0, sech_refe = 0.d0
! ----- Hydratation
        real(kind=8)  :: hydr_prev = 0.d0, hydr_curr = 0.d0, hydr_incr = 0.d0
! ----- Hygrometry
        real(kind=8)  :: hygr_prev = 0.d0, hygr_curr = 0.d0, hygr_incr = 0.d0
! ----- Total pressure
        real(kind=8)  :: ptot_prev = 0.d0, ptot_curr = 0.d0, ptot_incr = 0.d0
! ----- Non-mechanical strains (external state variable as temperature)
        real(kind=8)  :: epsi_varc(6) = 0.d0, depsi_varc(6) = 0.d0
! ----- Thermic strains
        real(kind=8)  :: epsthm = 0.d0, epsth_anism(3) = 0.d0, epsth_metam = 0.d0
        real(kind=8)  :: epsthp = 0.d0, epsth_anisp(3) = 0.d0, epsth_metap = 0.d0
    end type Behaviour_ESVA
!
! - Behaviour - Integration - Parameters on element
!
    type Behaviour_Elem
! ----- Flags
        aster_logical :: l_eltsize1 = ASTER_FALSE
! ----- Size of element
        real(kind=8)  :: eltsize1 = 0.d0
! ----- Size of element for ENDO_PORO_BETON
        real(kind=8)  :: eltsize2(9) = 0.d0
! ----- Gradient of velocity for *CRISTAL
        real(kind=8)  :: gradvelo(9) = 0.d0
! ----- Coordinates of all Gauss points
        real(kind=8)  :: coor_elga(27, 3) = 0.d0
    end type Behaviour_Elem
!
! - Behaviour - Integration - Parameters on current Gauss point
!
    type Behaviour_Elga
! ----- For *_JOINT_HYME models : kinematic matrix
        real(kind=8) :: rotpg(3*3) = 0.d0
! ----- For CABLE_GAINE elements : tension of the cable
        real(kind=8) :: tenscab = 0.d0
! ----- For CABLE_GAINE elements : curvature of the cable
        real(kind=8) :: curvcab = 0.d0
! ----- For GRAD_VARI models : non-local variables PHI
        real(kind=8) :: nonloc(2) = 0.d0
! ----- For CZM_*_MIX behaviours : Lagrange penalty coefficient
        real(kind=8) :: r = 0.d0
    end type Behaviour_Elga
!
! - Behaviour - Integration - Parameters for external solver (UMAT, MFRONT)
!
    type Behaviour_Exte
! ----- Number of external state variables used in external solver
        integer            :: nb_pred = 0
! ----- Value of external state variables used in external solver
        real(kind=8)       :: predef(EXTE_ESVA_NBMAXI) = 0.d0
! ----- Incremental value of external state variables used in external solver
        real(kind=8)       :: dpred(EXTE_ESVA_NBMAXI) = 0.d0
    end type Behaviour_Exte
!
! - Behaviour - Integration
!
    type Behaviour_Integ
! ----- Current time
        real(kind=8)          :: time_curr = 0.d0
! ----- Coded integer to detect external state variables
        integer               :: tabcod(60) = 0
! ----- Parameters on Gauss point
        type(Behaviour_Elga)  :: elga
! ----- Parameters on element
        type(Behaviour_Elem)  :: elem
! ----- Parameters for external solvers (MFRONT, UMAT)
        type(Behaviour_Exte)  :: exte
! ----- Parameters for external state variables
        type(Behaviour_ESVA)  :: esva
! ----- Flag when GEOM external state variable is present
        aster_logical         :: l_varext_geom = ASTER_FALSE
! ----- Flag if MFront is used
        aster_logical         :: l_mfront = ASTER_FALSE
! ----- Flag if UMAT is used
        aster_logical         :: l_umat = ASTER_FALSE
    end type Behaviour_Integ
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour: types for preparation of behaviour (maps COMPOR and CARCRI)
!
! --------------------------------------------------------------------------------------------------

! - Behaviour - Preparation - Parameters for external behaviours
    type Behaviour_ParaExte
! ----- Flag for UMAT law
        aster_logical      :: l_umat = ASTER_FALSE
! ----- Flag for non-official MFront law
        aster_logical      :: l_mfront_proto = ASTER_FALSE
! ----- Flag for official MFront law
        aster_logical      :: l_mfront_offi = ASTER_FALSE
! ----- Type of behaviour: 0 (internal integration), 1 (MFront official),
!       2 (MFront proto), 4 (UMAT)
        integer            :: extern_type = 0
! ----- Address to MGISBehaviour object as hexadecimal
        character(len=16)  :: extern_addr = ' '
! ----- Address to UMAT function
        integer            :: extern_ptr = 0
! ----- Name of subroutine for external UMAT law
        character(len=255) :: subr_name = ' '
! ----- Name of library for external UMAT law
        character(len=255) :: libr_name = ' '
! ----- Model for MFront law
        integer            :: model_mfront = MFRONT_MODEL_UNSET
! ----- Number of dimension for MFront law
        integer            :: model_dim = 0
! ----- Number of internal variables for UMAT
        integer            :: nbVariUMAT = 0
! ----- Identifier for strains model
        integer            :: strain_model = MFRONT_STRAIN_UNSET
    end type Behaviour_ParaExte

! - Behaviour - Preparation - Parameters for behaviour
    type Behaviour_Para
! ----- Keyword RELATION
        character(len=16) :: rela_comp = ' '
! ----- Keyword DEFORMATION
        character(len=16) :: defo_comp = ' '
! ----- Keyword COMP_INCR/COMP_ELAS
        character(len=16) :: type_comp = ' '
! ----- Keyword DEBORST
        character(len=16) :: type_cpla = ' '
! ----- Keyword KIT
        character(len=16) :: kit_comp(4) = ' '
! ----- Keyword COMPOR
        character(len=16) :: mult_comp = ' '
! ----- Keyword POST_ITER
        character(len=16) :: post_iter = ' '
! ----- Type of strain transmitted to the behaviour law : 'OLD', 'MECANIQUE' or 'TOTALE'
        character(len=16) :: defo_ldc = ' '
! ----- Index of law
        integer           :: numeLaw = 0
! ----- Total number of internal state variables
        integer           :: nbVari = 0
! ----- Number of internal state variables for kit
        integer           :: nbVariKit(4) = 0
! ----- Index of law for kit
        integer           :: numeLawKit(4) = 0
! ----- Keyword RIGI_GEOM
        character(len=16) :: rigi_geom = ' '
! ----- Keyword REGU_VISC
        character(len=16) :: regu_visc = ' '
! ----- Mechanical part of behaviour
        character(len=16) :: meca_comp = ' '
! ----- Keyword POST_INCR
        character(len=16) :: post_incr = ' '
    end type Behaviour_Para

! - Behaviour - Preparation - Criteria for behaviour
    type Behaviour_Crit
! ----- Keyword RELATION
        character(len=16) :: rela_comp = ' '
! ----- Mechanical part of behaviour
        character(len=16) :: meca_comp = ' '
! ----- Parameters for external behaviours
        type(Behaviour_ParaExte)  :: paraExte
! ----- Criteria
        integer                   :: type_matr_t = 0
        real(kind=8)              :: parm_theta = 0.d0
        integer                   :: iter_inte_pas = 0
        real(kind=8)              :: vale_pert_rela = 0.d0
        real(kind=8)              :: resi_deborst_max = 0.d0
        integer                   :: iter_deborst_max = 0
        real(kind=8)              :: resi_radi_rela = 0.d0
        integer                   :: ipostiter = 0
        integer                   :: iveriborne = 0
        aster_logical             :: l_matr_unsymm = ASTER_FALSE
        real(kind=8)              :: algo_inte_r = 0.d0
        real(kind=8), pointer     :: resi_inte => null()
        integer, pointer          :: iter_inte_maxi => null()
        integer                   :: extern_ptr = 0
        integer                   :: extern_type = 0
        integer                   :: exte_strain = 0
        integer                   :: jvariext1 = 0
        integer                   :: jvariext2 = 0
    end type Behaviour_Crit

! - Behaviour - Preparation - Map for criteria of behaviours (CARCRI)
    type Behaviour_PrepCrit
! ----- Number of factor keywords
        integer                           :: nb_comp = 0
! ----- Parameters for THM scheme
        real(kind=8)                      :: parm_alpha_thm = 0.d0
        real(kind=8)                      :: parm_theta_thm = 0.d0
! ----- List of criteria (by keyword COMPORTEMENT)
        type(Behaviour_Crit), pointer     :: v_crit(:)
    end type Behaviour_PrepCrit

! - Behaviour - Preparation - Map for parameters of behaviours (COMPOR)
    type Behaviour_PrepPara
! ----- Number of factor keywords
        integer :: nb_comp = 0
! ----- List of parameters
        type(Behaviour_Para), pointer :: v_para(:) => null()
! ----- List of parameters for external behaviours
        type(Behaviour_ParaExte), pointer  :: v_paraExte(:) => null()
! ----- Flag for total strain model cases (at least one behaviour)
        aster_logical :: lTotalStrain = ASTER_FALSE
! ----- Flag for debug
        aster_logical :: lDebug = ASTER_FALSE
    end type Behaviour_PrepPara
!
end module
