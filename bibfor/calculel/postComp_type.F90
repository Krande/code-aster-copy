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
! ==================================================================================================
!
! Types to manage the calculation of a post-processing option
!
! ==================================================================================================
!
module postComp_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - Parameters of results
! ==================================================================================================
    type POST_COMP_RESU
! ----- Input result
        character(len=8) :: resultIn = " "
! ----- Output result
        character(len=8) :: resultOut = " "
! ----- Type of results
        character(len=16) :: resultType = " "
! ----- Type of MODE_MECA
        character(len=16) :: modeType = " "
! ----- Number of storing indexes to compute
        integer(kind=8) :: nbStore = 0
! ----- List of storing indexes to compute
        integer(kind=8), pointer :: listStore(:) => null()
! ----- JEVEUX name of datastructure for list of storing indexes
        character(len=19) :: listStoreJv = " "
! ----- Flag for transient result
        aster_logical :: lTransient = ASTER_FALSE
! ----- Flag for thermal result
        aster_logical :: lTher = ASTER_FALSE
! ----- Flag for mode result
        aster_logical :: lMode = ASTER_FALSE
    end type POST_COMP_RESU
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - Input fields
! ==================================================================================================
    type POST_COMP_FIELDS
! ----- Number of equations
        integer(kind=8) :: nbEqua
! ----- Zero field
        character(len=24) :: vectZero = "&&CCFNRN.ZERO"
! ----- Input fields: DEPL, VITE, ACCE, SIEF_ELGA, STRX_ELGA
        character(len=24) :: disp = " "
        aster_logical :: hasVite = ASTER_FALSE
        character(len=24) :: vite = " "
        aster_logical :: hasAcce = ASTER_FALSE
        character(len=24) :: acce = " "
        character(len=24) :: sigm = " "
        character(len=24) :: sigmPrev = " "
        character(len=24) :: strx = " "
    end type POST_COMP_FIELDS
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - General parameters
! ==================================================================================================
    type POST_COMP_PARA
! ----- Model, elementary characteristics and materiel fields
        character(len=8) :: model = " "
        character(len=8) :: caraElem = " "
        character(len=8) :: materField = " "
        character(len=24) :: materCode = " "
! ----- External state variables
        character(len=19) :: chvarc = "&&CCFNRN.CHVARC"
! ----- Map for non-linear behaviour
        character(len=24) :: compor = " "
! ----- Parameters for loads
        aster_logical :: lLoadsHarmo = ASTER_FALSE
        aster_logical :: lLoadsFromUser = ASTER_FALSE
        aster_logical :: lLoadsAreDefined = ASTER_FALSE
        aster_logical :: hasPiloLoads = ASTER_FALSE
        character(len=24) :: listLoad = " "
! ----- Flag: model contains STRX fields
        aster_logical :: lElemStrx = ASTER_FALSE
! ----- Flag: model contains STRX fields that can be compute
        aster_logical :: lElemStrxComp = ASTER_FALSE
! ----- Flag: where to compute option
        aster_logical :: lReduComp = ASTER_FALSE
        character(len=24) :: calcLigrel = " "
! ----- Reference numbering
        character(len=24) :: numeDofRefe = " "
! ----- Adress of mass matrix
        integer(kind=8) :: jvMassMatr = 0
! ----- Fourier harmonic
        integer(kind=8) :: nh = -1
! ----- Value for frequency access
        aster_logical :: hasFreq = ASTER_FALSE
        real(kind=8) :: freq = -1.d0
! ----- Value for time access
        aster_logical :: hasTime = ASTER_FALSE
        real(kind=8) :: time = 0.d0
        real(kind=8) :: timePrev = 0.d0
! ----- Value for pulsation
        aster_logical :: hasOmega = ASTER_FALSE
        real(kind=8) :: omega = -1.d0
! ----- Value for "PILOTAGE"
        aster_logical :: hasEta = ASTER_FALSE
        real(kind=8) :: eta = 0.d0
    end type POST_COMP_PARA
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - Parameters for nodal fields
! ==================================================================================================
    type POST_COMP_NODA
! ----- Compute FORC_NODA option
        aster_logical :: lForcNoda = ASTER_FALSE
! ----- Compute FORC_EXTE option
        aster_logical :: lFextNoda = ASTER_FALSE
! ----- Compute M_GAMMA option
        aster_logical :: lMGamma = ASTER_FALSE
! ----- Compute REAC_NODA
        aster_logical :: lReacNoda = ASTER_FALSE
! ----- Coefficients for FEXT_NODA
        real(kind=8) :: coefFextR = 0.d0
        complex(kind=8) :: coefFextC = dcmplx(0.d0, 0.d0)
! ----- Coefficients for M_GAMMA
        real(kind=8) :: coefMGamR = 0.d0
        complex(kind=8) :: coefMGamC = dcmplx(0.d0, 0.d0)
    end type POST_COMP_NODA
! ==================================================================================================
! Type: Management of computation of option for post-processing
! ==================================================================================================
    type POST_COMP
! ----- Parameters for result management
        type(POST_COMP_RESU) :: postCompResu
! ----- General parameters to compute option
        type(POST_COMP_PARA) :: postCompPara
! ----- Input fields
        type(POST_COMP_FIELDS) :: postCompFields
! ----- Parameters for nodal quantities of post-treatment
        type(POST_COMP_NODA) :: postCompNoda
! ----- Scalar type for output option (Real or Complex)
        character(len=1) :: scalType = " "
        aster_logical :: lCplx = ASTER_FALSE
! ----- Need list of loads
        aster_logical :: lNeedLoads = ASTER_FALSE
! ----- Compute nodal quantities for post-treatment
        aster_logical :: lPostNoda = ASTER_FALSE
    end type POST_COMP
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - Management of restricted zone
! ==================================================================================================
    type POST_COMP_REST
! ----- Number of FED
        integer(kind=8) :: nbLigrMaxi = 0
        integer(kind=8) :: nbLigr = 0
! ----- List of FED
        character(len=24), pointer :: listFED(:) => null()
! ----- List of models
        character(len=8), pointer :: listModel(:) => null()
! ----- List of memory base
        character(len=8), pointer :: listJvBase(:) => null()
    end type POST_COMP_REST
! ==================================================================================================
! Type: Manage the calculation of a post-processing option - Special for POUX beams
! ==================================================================================================
    type POST_COMP_POUX
! ----- POUX beams in model
        aster_logical :: lPoux = ASTER_FALSE
! ----- Index of distributed load in list of loads for POUX
        integer(kind=8) :: loadIndx = 0
! ----- Name of option to compute mass matrix
        character(len=24) :: optionMass = " "
! ----- Name of vector for M.Gamma
        character(len=24) :: chdynr = "&&MECALM.M.GAMMA"
    end type POST_COMP_POUX
!===================================================================================================
    public :: POST_COMP, POST_COMP_RESU, POST_COMP_PARA, POST_COMP_FIELDS, POST_COMP_NODA
    public :: POST_COMP_REST, POST_COMP_POUX
!
end module postComp_type
