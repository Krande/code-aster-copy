! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation  either version 3 of the License  or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not  see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!

! --------------------------------------------------------------------------------------------------
!
! For external state variables
!
! --------------------------------------------------------------------------------------------------

! Maximum number of differents types of external state variables for MFront/MGIS
#define ESVA_EXTE_MGIS_NBMAXI    11

! --------------------------------------------------------------------------------------------------
!
! For external solver MFront/MGIS
!
! Constants
!
! --------------------------------------------------------------------------------------------------
!
! Type of model and strain model (MGIS)
! keep consistency with bibcxx/Behaviours/MGISBehaviourFort.h
!
#define MGIS_MODEL_UNSET          0
#define MGIS_MODEL_TRIDIMENSIONAL 1
#define MGIS_MODEL_AXISYMMETRICAL 2
#define MGIS_MODEL_PLANESTRESS    3
#define MGIS_MODEL_PLANESTRAIN    4
!
#define MGIS_STRAIN_UNSET         0
#define MGIS_STRAIN_SMALL         1
#define MGIS_STRAIN_F             3
!
#define MGIS_BV_SCALAR            0
#define MGIS_BV_VECTOR            2
#define MGIS_BV_VECTOR_1D         10
#define MGIS_BV_VECTOR_2D         18
#define MGIS_BV_VECTOR_3D         26
#define MGIS_BV_STENSOR           1
#define MGIS_BV_STENSOR_1D        9
#define MGIS_BV_STENSOR_2D        17
#define MGIS_BV_STENSOR_3D        25
#define MGIS_BV_TENSOR            3
#define MGIS_BV_TENSOR_1D         11
#define MGIS_BV_TENSOR_2D         19
#define MGIS_BV_TENSOR_3D         27
#define MGIS_BV_HIGHER_ORDER_TENSOR   4
#define MGIS_BV_ARRAY                 5

! Maximum number of material properties
#define MGIS_MAX_PROPS            197

! Type of operator
#define MGIS_BV_PREDICTION_TANGENT_OPERATOR -3
#define MGIS_BV_PREDICTION_SECANT_OPERATOR -2
#define MGIS_BV_PREDICTION_ELASTIC_OPERATOR -1
#define MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR 0
#define MGIS_BV_INTEGRATION_ELASTIC_OPERATOR 1
#define MGIS_BV_INTEGRATION_SECANT_OPERATOR 2
#define MGIS_BV_INTEGRATION_TANGENT_OPERATOR 3
#define MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR 4
