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
!
module Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------
!

! - Type: parameters for behaviour
    type META_Parameters
! ----- Keyword RELATION (steel, zirc, etc.)
        character(len=16) :: metaType = ' '
! ----- Keyword LOI_META
        character(len=16) :: metaLaw = ' '
! ----- Total number of internal state variables
        integer :: nbVari = 0
! ----- Number of phases
        integer :: nbPhase = 0
! ----- Index of behaviour
        integer :: numeComp = 0
    end type META_Parameters

! - Metallurgy - Preparation - Map for parameters of behaviours (COMPOR_META)
    type META_PrepPara
! ----- Number of factor keywords
        integer :: nb_comp = 0
! ----- List of parameters
        type(META_Parameters), pointer :: para(:) => null()
    end type META_PrepPara

! - Parameters for austenite phase
    type META_AusteniteParameters
        real(kind=8) :: lambda0 = 0.d0
        real(kind=8) :: qsr_k = 0.d0
        real(kind=8) :: d10 = 0.d0
        real(kind=8) :: wsr_k = 0.d0
    end type META_AusteniteParameters
!
    type META_TRCMartensiteLaw
        real(kind=8) :: austeniteMin = 0.d0
        real(kind=8) :: akm = 0.d0, bkm = 0.d0
        real(kind=8) :: lowerSpeed = 0.d0
    end type META_TRCMartensiteLaw
!
    type META_TRCAusteniteGrain
        real(kind=8) :: dref = 0.d0
        real(kind=8) :: a = 0.d0
    end type META_TRCAusteniteGrain
!
    type META_TRCParameters
        integer :: jv_ftrc = 0, jv_trc = 0
        integer :: iadexp = 0, iadtrc = 0
        integer :: nbHist = 0
        type(META_TRCMartensiteLaw) :: martensiteLaw
        type(META_TRCAusteniteGrain) :: austeniteGrain
    end type META_TRCParameters
!
    type META_SteelParameters
        real(kind=8) :: ar3 = 0.d0
        real(kind=8) :: alpha = 0.d0
        real(kind=8) :: ms0 = 0.d0
! Quasi-static temperature at which austenite transformation begins on heating.
        real(kind=8) :: ac1 = 0.d0
        ! Quasi-static temperature at end of austenite transformation
        real(kind=8) :: ac3 = 0.d0
        real(kind=8) :: taux_1 = 0.d0
        real(kind=8) :: taux_3 = 0.d0
        aster_logical :: l_grain_size = ASTER_FALSE
        type(META_AusteniteParameters) :: austenite
        type(META_TRCParameters) :: trc
    end type META_SteelParameters
!
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
!
    type META_MaterialParameters
        type(META_SteelParameters) :: steel
        type(META_ZircParameters) :: zirc
    end type META_MaterialParameters
!
end module
