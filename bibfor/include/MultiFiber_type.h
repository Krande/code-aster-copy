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
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
! --------------------------------------------------------------------------------------------------
!
!   Pour un groupe de fibre
!       MULTI_FIBER_NAME    : nom du groupe
!       MULTI_FIBER_MATER   : matériau
!       MULTI_FIBER_RELA    : relation
!       MULTI_FIBER_ALGO    : algo (ANALYTIQUE, DEBORST, 1D, ... )
!       MULTI_FIBER_DEFO    : la déformation (PETIT, ...)
!       MULTI_FIBER_NBFI    : nombre de fibres
!       MULTI_FIBER_NBVARI  : nombre de variable internes sur le groupe
!
! --------------------------------------------------------------------------------------------------
!
! Si le nombre de slots change, il faut également mettre à jour
!       "sd_compor.py"
!
!   Size for CPRK
#define MULTI_FIBER_SIZEK    7
!   Slots for CPRK
#define MULTI_FIBER_NAME     1
#define MULTI_FIBER_MATER    2
#define MULTI_FIBER_RELA     3
#define MULTI_FIBER_ALGO     4
#define MULTI_FIBER_DEFO     5
#define MULTI_FIBER_NBFI     6
#define MULTI_FIBER_NBVARI   7
!
!   Size for CPRI
#define MULTI_FIBER_SIZEI    3
!   Slots for CPRI
#define MULTI_FIBER_TYPE     1
#define MULTI_FIBER_NBVARMAX 2
#define MULTI_FIBER_NBGRFIBR 3
!
