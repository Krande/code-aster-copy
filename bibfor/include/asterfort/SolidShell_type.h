! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
! Type of cell
#define SSH_CELL_UNDEF     0
#define SSH_CELL_HEXA      1

! Number of integration points
#define SSH_NBPG_MAX       5

! Size of tensors
#define SSH_SIZE_TENS      6

! Size of matrices
#define SSH_SIZE_MATR_HEXA 325

! Size of spatial frame
#define SSH_NDIM           3

! Number of nodes
#define SSH_NBNODEG_MAX    8
#define SSH_NBNODE_MAX     9
#define SSH_NBNODEG_HEXA   8
#define SSH_NBNODE_HEXA    9

! Number of DOF
#define SSH_NBDOFG_MAX     24
#define SSH_NBDOF_MAX      25
#define SSH_NBDOFG_HEXA    24
#define SSH_NBDOF_HEXA     25

! For decomposition of jacobian
#define SSH_JACO_J1(J10,J1ETA,J1ZETA,J1ETAZETA,XI2,XI3) J10+XI2*J1ETA+XI3*J1ZETA+XI2*XI3*J1ETAZETA
#define SSH_JACO_J2(J20,J2XI ,J2ZETA,J2XIZETA ,XI1,XI3) J20+XI1*J2XI +XI3*J2ZETA+XI1*XI3*J2XIZETA
#define SSH_JACO_J3(J30,J3ETA,J3XI  ,J3XIETA  ,XI1,XI2) J30+XI2*J3ETA+XI1*J3XI  +XI1*XI2*J3XIETA

! For DEBUG
#define SSH_DBG_UNIT       6
#define SSH_DBG_ELEM       ASTER_FALSE
#define SSH_DBG_MATE       ASTER_FALSE
#define SSH_DBG_GEOM       ASTER_FALSE
#define SSH_DBG_KINE       ASTER_FALSE
#define SSH_DBG_STAB       ASTER_FALSE
#define SSH_DBG_BEHA       ASTER_FALSE
#define SSH_DBG_STRG(strg) WRITE(SSH_DBG_UNIT,'(A)')  strg
#define SSH_DBG_INTE(inte) WRITE(SSH_DBG_UNIT,'(I6)') inte
