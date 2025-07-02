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
! Type of beam
#define BEAM_TYPE_UNDEF        -1
#define BEAM_TYPE_STRAIGHT     0
#define BEAM_TYPE_ELBOW        1

! Type of section
#define BEAM_SECT_UNDEF        -1
#define BEAM_SECT_PIPE         0

! Discretization
#define BEAM_MAX_NODE         4

! Index of DOF
#define BEAM_NBDOF            6
#define BEAM_DOF_DX           1
#define BEAM_DOF_DY           2
#define BEAM_DOF_DZ           3
#define BEAM_DOF_DRX          4
#define BEAM_DOF_DRY          5
#define BEAM_DOF_DRZ          6
