! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! CalcG module : Parameters <-> integer definitions
! -------------------------------------------------------------------------
!
! - CALC_G - General
!
! --- number of option computed
#define NB_MAX_OPT 5
!
! --- list of options
#define OPT_G 1
#define OPT_K 2
#define OPT_G_EPSI 3
#define OPT_KJ 4
#define OPT_KJ_EPSI 5
!
! --- number of option computed
#define NB_MAX_TERM 8
!
! --- list of terms
#define G  1
#define K1 2
#define K2 3
#define K3 4
#define G_IRWIN 5
#define G_EPSI 6
#define KJ 7 
#define KJ_EPSI 8 
!
! --- number of parameters for table
#define NB_MAX_PARA 20
!
! --- List of comportement
#define ELAS    0
#define ELAS_NL 1
#define PLAS    2
