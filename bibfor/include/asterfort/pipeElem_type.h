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
! Flag for Debug
#define PIPE_DEBUG            ASTER_FALSE

! Type of pipe
#define PIPE_TYPE_UNDEF        -1
#define PIPE_TYPE_STRAIGHT     0
#define PIPE_TYPE_ELBOW        1

! Discretization
#define PIPE_TENS_SIZE        4
#define PIPE_MAX_NODE         4
#define PIPE_MAX_NPG          4
#define PIPE_MAX_LAYERS       10
#define PIPE_MAX_SECTORS      32
#define PIPE_MAX_DOF          156

! Index of DOF (beam DOF)
#define PIPE_NBDOF_BEAM       6
#define PIPE_DOF_DX           1
#define PIPE_DOF_DY           2
#define PIPE_DOF_DZ           3
#define PIPE_DOF_DRX          4
#define PIPE_DOF_DRY          5
#define PIPE_DOF_DRZ          6

! Index of components for DEGE field
#define PIPE_NB_DEGE          6
#define PIPE_DEGE_EPXX        1
#define PIPE_DEGE_GAXY        2
#define PIPE_DEGE_GAXZ        3
#define PIPE_DEGE_GAT         4
#define PIPE_DEGE_KY          5
#define PIPE_DEGE_KZ          6

! Index of components for EFGE field
#define PIPE_NB_EFGE          6

! Index of components for EPSI field
#define PIPE_NB_EPSI          6
#define PIPE_EPSI_EPXX        1
#define PIPE_EPSI_EPYY        2
#define PIPE_EPSI_EPZZ        3
#define PIPE_EPSI_EPXY        4
#define PIPE_EPSI_EPXZ        5
#define PIPE_EPSI_EPYZ        6

! Index of components for SIGM field
#define PIPE_NB_SIGM          6
#define PIPE_SIGM_SIXX        1
#define PIPE_SIGM_SIYY        2
#define PIPE_SIGM_SIZZ        3
#define PIPE_SIGM_SIXY        4
#define PIPE_SIGM_SIXZ        5
#define PIPE_SIGM_SIYZ        6

! Index of components for COOR_ELGA
#define PIPE_NB_COORELGA      4

! Tolerance for MODI_METRIC
#define PIPE_METRIC_LIMIT     0.25
