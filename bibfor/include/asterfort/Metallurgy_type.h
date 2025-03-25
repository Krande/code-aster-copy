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
! Metallurgy data structure
! -------------------------------------------------------------------------
!
#define META_NONE        0
#define META_STEEL       1
#define META_ZIRC        2

! - Phases for steel: total number (4 cold, 1 hot + 1 sum of cold)
#define NBPHASESTEEL     6

! - Phases for standard steel
#define PSTEEL_NB        5
#define PFERRITE         1
#define PPERLITE         2
#define PBAINITE         3
#define PMARTENS         4
#define PAUSTENITE       5
#define PSUMCOLD         6
! For next ones: add total number of phases to access in internal state variable vector
#define SIZE_GRAIN       1
#define STEEL_TEMP       2
#define TEMP_MARTENSITE  3
! Number of  internal states variables required for initial state
#define PVARIINIT        7

! For steel, law : waeckel
#define NBVARIWAECKEL    9

! For steel, law : jma
#define NBVARIJMA        12

! - Phases for zircaloy
#define PZIRC_NB         3
#define PALPHA1          1
#define PALPHA2          2
#define PBETA            3
! For next ones: add total number of phases to access in internal state variable vector
#define ZIRC_TEMP        1
#define TIME_TRAN        2

! - Phases for steel with tempering: total number (6 cold, 1 hot + 1 sum of cold)
#define NBPHASESTEELR     8

! - Phases for steel with tempering
#define PRSTEEL_NB        7
#define PRFERRITE         1
#define PRPERLITE         2
#define PRBAINITEB        3
#define PRBAINITER        4
#define PRMARTENSB        5
#define PRMARTENSR        6
#define PRAUSTENITE       7
#define PRSUMCOLD         8
! Same as standard steel:
! #define SIZE_GRAIN        1
! #define STEEL_TEMP        2
! #define TEMP_MARTENSITE   3
#define THER_CYCL         4
! Number of  internal states variables required for initial state
#define PRVARIINIT        9

! - Index for internal variables
#define IDX_I_EPSEQ      7
#define IDX_I_IPLAS      8
#define IDX_C_IPLAS      43

! - Kinetic
#define COOLING          0
#define HEATING          1

!
! --------------------------------------------------------------------------------------------------
!
! The field <COMPOR_META>
!   List of strings (K16) for each cell
!   Definition of behaviour (relation, strain model, number of internal state vari, etc.)
!
! --------------------------------------------------------------------------------------------------
! Size
#define COMPORMETA_SIZE  5

! - Slot in COMPOR_META map
#define ZMETATYPE        1
#define ZNBVARI          2
#define ZMETALAW         3
#define ZNUMECOMP        4
#define ZNBPHASE         5

! --------------------------------------------------------------------------------------------------
!
! Parameters for integration of metallurgical laws in CALC_META
!
! --------------------------------------------------------------------------------------------------
! - Bound to cut
#define STEEL_MIN_CUT 1.d-3
! - Maximum temperature step
#define STEEL_MAX_TEMP_STEP 5.d0
! - Maximum number of sub-steps
#define STEEL_MAX_NB_STEP 10
