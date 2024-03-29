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
! Metallurgy data structure
! -------------------------------------------------------------------------
!
#define META_NONE        0
#define META_STEEL       1
#define META_ZIRC        2

! - Phases for steel
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

! - Phases for zircaloy
#define PZIRC_NB         3
#define PALPHA1          1
#define PALPHA2          2
#define PBETA            3
! For next ones: add total number of phases to access in internal state variable vector
#define ZIRC_TEMP        1
#define TIME_TRAN        2

! - Phases for steel with revenu
#define PRSTEEL_NB        7
#define PRFERRITE         1
#define PRPERLITE         2
#define PRBAINITEB        3
#define PRBAINITER        4
#define PRMARTENSB        5
#define PRMARTENSR        6
#define PRAUSTENITE       7
#define PRSUMCOLD         8
! For next ones: add total number of phases to access in internal state variable vector
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

! - Slot in COMPOR_META map
#define ZMETATYPE        1
#define ZNBVARI          2
#define ZMETALAW         3
#define ZNUMECOMP        4
#define ZNBPHASE         5
