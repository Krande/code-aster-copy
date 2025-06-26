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
module Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Metallurgy_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------
!

! - Metallurgy - Parameters for operator
    type META_ParaOperator
! ----- Name of output datastructure
        character(len=8) :: resultName = " "
! ----- For storing indexes
        character(len=19) :: listStoreJv = '&&OP0194.LISTSTORE'
        integer(kind=8) :: nbStore = 0
        integer(kind=8), pointer :: listStore(:) => null()
! ----- List of options to compute
        character(len=19) :: listOptionsJv = '&&OP0194.LES_OPTION'
        integer(kind=8) :: nbOption = 0
        character(len=16), pointer :: listOption(:) => null()
! ----- Main parameters
        character(len=24) :: metaLigrel = "&&OP0194.METALIGREL"
        character(len=24) :: modelLigrel = " "
        character(len=8) :: model = " "
        character(len=8) :: materialField = " "
        character(len=24) :: materialCoding = " "
        character(len=24) :: comporMetaTemper = "&&OP0194.COMPORTEMPER"
        character(len=24) :: comporMeta = "&&OP0194.COMPOR"
        aster_logical :: hasTRC = ASTER_FALSE
        character(len=24) :: TRCField = "&&SMEVOL.ADRESSES"
! ----- Flag for tempering
        aster_logical :: hasTemper = ASTER_FALSE
    end type META_ParaOperator

! - Metallurgy - Parameters for behaviour
    type META_ParaBehaviour
! ----- Keyword RELATION (steel, zirc, etc.)
        character(len=16) :: metaType = ' '
! ----- Keyword LOI_META
        character(len=16) :: metaLaw = ' '
! ----- Total number of internal state variables
        integer(kind=8) :: nbVari = 0
! ----- Number of phases
        integer(kind=8) :: nbPhase = 0
! ----- Index of behaviour
        integer(kind=8) :: numeComp = 0
    end type META_ParaBehaviour

! - Metallurgy - Preparation - Map for parameters of behaviours (COMPOR_META)
    type META_PrepBehaviour
! ----- Factor keyword to read
        character(len=16) :: factorKeyword = " "
! ----- Number of factor keywords
        integer(kind=8) :: nbFactorKeyword = 0
! ----- List of parameters
        type(META_ParaBehaviour), pointer :: paraBehaviour(:) => null()
! ----- Flag for tempering
        aster_logical :: hasTemper = ASTER_FALSE
    end type META_PrepBehaviour

! - Metallurgy - Parameters for austenite phase
    type META_AusteniteParameters
        real(kind=8) :: lambda0 = 0.d0
        real(kind=8) :: qsr_k = 0.d0
        real(kind=8) :: d10 = 0.d0
        real(kind=8) :: wsr_k = 0.d0
    end type META_AusteniteParameters
    type META_TRCAusteniteGrain
        real(kind=8) :: dref = 0.d0
        real(kind=8) :: a = 0.d0
    end type META_TRCAusteniteGrain

! - Metallurgy - Parameters for martensite phase
    type META_TRCMartensiteLaw
        real(kind=8) :: austeniteMin = 0.d0
        real(kind=8) :: akm = 0.d0, bkm = 0.d0
        real(kind=8) :: lowerSpeed = 0.d0
    end type META_TRCMartensiteLaw

! - Metallurgy - Parameters for tempering
    type META_TemperingParameters
        real(kind=8) :: bainite_b = 0.d0
        real(kind=8) :: bainite_n = 0.d0
        real(kind=8) :: martensite_b = 0.d0
        real(kind=8) :: martensite_n = 0.d0
        real(kind=8) :: temp = 0.d0
        real(kind=8) :: tempHold = 0.d0
    end type META_TemperingParameters

! - Metallurgy - Parameters for TRC curves
    type META_TRCParameters
        integer(kind=8) :: jv_ftrc = 0, jv_trc = 0
        integer(kind=8) :: iadexp = 0, iadtrc = 0
        integer(kind=8) :: nbHist = 0
        type(META_TRCMartensiteLaw) :: martensiteLaw
        type(META_TRCAusteniteGrain) :: austeniteGrain
    end type META_TRCParameters

! - Metallurgy - Parameters for steel
    type META_SteelParameters
        aster_logical :: lNodeDebug = ASTER_FALSE
        real(kind=8) :: ar3 = 0.d0
        real(kind=8) :: alpha = 0.d0
        real(kind=8) :: ms0 = 0.d0
! ----- Quasi-static temperature at which austenite transformation begins on heating.
        real(kind=8) :: ac1 = 0.d0
! ----- Quasi-static temperature at end of austenite transformation
        real(kind=8) :: ac3 = 0.d0
        real(kind=8) :: taux_1 = 0.d0
        real(kind=8) :: taux_3 = 0.d0
        aster_logical :: l_grain_size = ASTER_FALSE
        type(META_AusteniteParameters) :: austenite
        type(META_TemperingParameters) :: temper
        type(META_TRCParameters) :: trc
    end type META_SteelParameters

! - Metallurgy - Parameters for zircaloy
    type META_ZircParameters
        real(kind=8) :: tdeq = 0.d0
        real(kind=8) :: k = 0.d0
        real(kind=8) :: n = 0.d0
        real(kind=8) :: t1c = 0.d0
        real(kind=8) :: t2c = 0.d0
        real(kind=8) :: ac = 0.d0
        real(kind=8) :: m = 0.d0
        real(kind=8) :: qsrk = 0.d0
        real(kind=8) :: t1r = 0.d0
        real(kind=8) :: t2r = 0.d0
        real(kind=8) :: ar = 0.d0
        real(kind=8) :: br = 0.d0
    end type META_ZircParameters

! - Metallurgy - Parameters for material
    type META_MaterialParameters
        type(META_SteelParameters) :: steel
        type(META_ZircParameters) :: zirc
    end type META_MaterialParameters

! - Metallurgy - Parameters for hardness
    type META_HardnessParameters
        character(len=16) :: metaType = ' '
        real(kind=8) :: hardSteel(PRSTEEL_NB) = 0.d0
    end type META_HardnessParameters
!
end module
