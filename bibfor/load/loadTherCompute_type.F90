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
module loadTherCompute_type
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
! ==================================================================================================
!
! Global variables - General
!
! Definition of parameters for all AFFE_CHAR_THER loads
!
! ==================================================================================================

! - Keyword to define load in AFFE_CHAR_THER
    character(len=24), parameter :: therLoadKeyword(LOAD_NEUT_NBTYPE) = (/ &
                                    'ECHANGE                 ', 'FLUX_REP_XYZ            ', &
                                    'FLUX_REP_NORM           ', 'SOURCE                  ', &
                                    'ECHANGE_PAROI           ', 'PRE_GRAD_TEMP           ', &
                                    'FLUX_NL                 ', 'SOUR_NL                 ', &
                                    'RAYONNEMENT             ', 'EVOL_CHAR               ', &
                                    'SOURCE_CALCULEE         '/)

! - Name of field to define load in AFFE_CHAR_THER
    character(len=6), parameter :: therLoadField(LOAD_NEUT_NBTYPE) = (/ &
                                   '.COEFH', '.FLURE', &
                                   '.FLUR2', '.SOURE', &
                                   '.HECHP', '.GRAIN', &
                                   '.FLUNL', '.SOUNL', &
                                   '.RAYO ', '.EVOL ', &
                                   '.SOURC'/)

! - Name of option (RHS) for real coefficient
    character(len=16), parameter :: therLoadVectR(LOAD_NEUT_NBTYPE) = (/ &
                                    'CHAR_THER_ECHA_R', 'CHAR_THER_FLUN_R', &
                                    'CHAR_THER_FLUX_R', 'CHAR_THER_SOUR_R', &
                                    'CHAR_THER_PARO_R', 'CHAR_THER_GRAI_R', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'CHAR_EVOL_CHAR  ', &
                                    'CHAR_THER_SOUR_R'/)

! - Name of option (RHS) for function coefficient
    character(len=16), parameter :: therLoadVectF(LOAD_NEUT_NBTYPE) = (/ &
                                    'CHAR_THER_ECHA_F', 'CHAR_THER_FLUN_F', &
                                    'CHAR_THER_FLUX_F', 'CHAR_THER_SOUR_F', &
                                    'CHAR_THER_PARO_F', 'CHAR_THER_GRAI_F', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'NoVector        '/)

! - Name of input parameter field for second member (real coefficient)
    character(len=8), parameter :: therLoadParaR(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHR ', 'PFLUXNR ', &
                                   'PFLUXVR ', 'PSOURCR ', &
                                   'PHECHPR ', 'PGRAINR ', &
                                   'NoInput ', 'NoInput ', &
                                   'NoInput ', 'NoInput ', &
                                   'PSOURCR '/)

! - Name of input parameter field for second member (function coefficient)
    character(len=8), parameter :: therLoadParaF(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHF ', 'PFLUXNF ', &
                                   'PFLUXVF ', 'PSOURCF ', &
                                   'PHECHPF ', 'PGRAINF ', &
                                   'NoInput ', 'NoInput ', &
                                   'NoInput ', 'NoInput ', &
                                   'NoInput '/)

! - Name of option (residual) for real coefficient
    character(len=16), parameter :: therLoadResiR(LOAD_NEUT_NBTYPE) = (/ &
                                    'RESI_THER_COEF_R', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'RESI_THER_PARO_R', 'NoVector        ', &
                                    'RESI_THER_FLUXNL', 'RESI_THER_SOURNL', &
                                    'RESI_THER_RAYO_R', 'RESI_EVOL_CHAR  ', &
                                    'NoVector        '/)

! - Name of option (residual) for function coefficient
    character(len=16), parameter :: therLoadResiF(LOAD_NEUT_NBTYPE) = (/ &
                                    'RESI_THER_COEF_F', 'NoVector        ', &
                                    'NoVector        ', 'NoVector        ', &
                                    'RESI_THER_PARO_F', 'NoVector        ', &
                                    'RESI_THER_FLUXNL', 'RESI_THER_SOURNL', &
                                    'RESI_THER_RAYO_F', 'NoVector        ', &
                                    'NoVector        '/)

! - Name of input parameter field for residual (real coefficient)
    character(len=8), parameter :: therResiParaR(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHR ', 'PFLUXNR ', &
                                   'PFLUXVR ', 'PSOURCR ', &
                                   'PHECHPR ', 'PGRAINR ', &
                                   'PFLUXNL ', 'PSOURNL ', &
                                   'PRAYONR ', 'NoInput ', &
                                   'PSOURCR '/)

! - Name of input parameter field for second member (function coefficient)
    character(len=8), parameter :: therResiParaF(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHF ', 'PFLUXNF ', &
                                   'PFLUXVF ', 'PSOURCF ', &
                                   'PHECHPF ', 'PGRAINF ', &
                                   'PFLUXNL ', 'PSOURNL ', &
                                   'PRAYONF ', 'NoInput ', &
                                   'NoInput '/)

! - Name of option for tangent matrix (real coefficient)
    character(len=16), parameter :: therLoadMatrR(LOAD_NEUT_NBTYPE) = (/ &
                                    'RIGI_THER_ECHA_R', "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    'MTAN_THER_PARO_R', "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    'MTAN_THER_RAYO_R', "NoMatrix        ", &
                                    "NoMatrix        "/)

! - Name of option for tangent matrix (function coefficient)
    character(len=16), parameter :: therLoadMatrF(LOAD_NEUT_NBTYPE) = (/ &
                                    'RIGI_THER_ECHA_F', "NoMatrix        ", &
                                    "NoMatrix        ", "NoMatrix        ", &
                                    'MTAN_THER_PARO_F', "NoMatrix        ", &
                                    'MTAN_THER_FLUXNL', 'MTAN_THER_SOURNL', &
                                    'MTAN_THER_RAYO_F', "NoMatrix        ", &
                                    "NoMatrix        "/)

! - Name of input parameter field for matrix (real coefficient)
    character(len=8), parameter :: therMatrParaR(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHR ', 'NoInput ', &
                                   'NoInput ', 'NoInput ', &
                                   'PHECHPR ', 'NoInput ', &
                                   'NoInput ', 'NoInput ', &
                                   'PRAYONR ', 'NoInput ', &
                                   'NoInput '/)

! - Name of input parameter field for matrix (function coefficient)
    character(len=8), parameter :: therMatrParaF(LOAD_NEUT_NBTYPE) = (/ &
                                   'PCOEFHF ', 'NoInput ', &
                                   'NoInput ', 'NoInput ', &
                                   'PHECHPF ', 'NoInput ', &
                                   'PFLUXNL ', 'PSOURNL ', &
                                   'PRAYONF ', 'NoInput ', &
                                   'NoInput '/)
!
end module loadTherCompute_type
