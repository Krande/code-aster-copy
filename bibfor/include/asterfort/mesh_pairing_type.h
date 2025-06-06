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

! ==================================================================================================
!
! Pairing/intersection
!
! ==================================================================================================
! Method of pairing
#define PAIR_FAST        1
#define PAIR_OLD         2
#define PAIR_ROBUST      3

! Maximum number of neighbours of a cell
#define MAX_NB_NEIGH 4

! Maximum number of intersection points
#define MAX_NB_INTE 8

! Maximum number of quadrature points
#define MAX_NB_QUAD 56

! ==================================================================================================
!
! Error code for pairing
!
! ERR_CELL_ORTH: cells are orthognal
! ERR_CELL_OOR : out of range (greater than DIST_RATIO)
! ERR_PAIR_PROJ: error during projection
! ERR_INTE_VOID: intersection is not correct
! ERR_PAIR_SLAV: projected slave nodes are not inside master cell
! ERR_PAIR_MAST: master nodes are not inside projected slave cell
! ERR_CELL_DEGE: cell is degenerated
!
! ==================================================================================================
#define ERR_PAIR_NONE 0
#define ERR_CELL_ORTH 1
#define ERR_CELL_OOR  2
#define ERR_PAIR_PROJ 3
#define ERR_INTE_VOID 4
#define ERR_PAIR_SLAV 5
#define ERR_PAIR_MAST 6
#define ERR_CELL_DEGE 7
