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
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1501
!
module Rom_Datastructure_type
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------

!
! - Datastructure for results datastructure
!
    type ROM_DS_Result
! ----- Type of result: evol_ther, evol_noli, etc .
        character(len=16)       :: resultType = ' '
! ----- Name of datastructure
        character(len=8)        :: resultName = ' '
! ----- Number of time steps saved in results
        integer(kind=8)                 :: nbStore = 0
! ----- Flag for reduced coordinates table
        aster_logical           :: lTablFromResu = ASTER_FALSE
! ----- Reference model
        character(len=8)        :: modelRefe = ' '
! ----- Reference behaviour
        character(len=24)       :: comporRefe = ' '
    end type ROM_DS_Result
!
! - Datastructure to save reduced coordinates
!
    type ROM_DS_TablReduCoor
! ----- Table in result datastructure
        type(NL_DS_TableIO)    :: tablResu
! ----- Flag if table is given by user
        aster_logical          :: lTablFromUser = ASTER_FALSE
! ----- Name of table when given by user
        character(len=8)       :: tablUserName = ' '
! ----- Values of reduced coordinates
        real(kind=8), pointer  :: coorRedu(:) => null()
    end type ROM_DS_TablReduCoor
!
! - Datastructure to select snapshots
!
    type ROM_DS_Snap
! ----- Number of snapshots
        integer(kind=8)           :: nbSnap = 0
! ----- List of snapshots
        integer(kind=8), pointer  :: listSnap(:) => null()
    end type ROM_DS_Snap
!
! - Parameters for lineic base numbering
!
    type ROM_DS_LineicNumb
! ----- Number of slices
        integer(kind=8)           :: nbSlice = 0
! ----- For each node => which slice ?
        integer(kind=8), pointer  :: numeSlice(:) => null()
! ----- For each node => which IN slice ?
        integer(kind=8), pointer  :: numeSection(:) => null()
! ----- Tolerance for separating nodes
        real(kind=8)      :: toleNode = 1.d-7
! ----- Number of components by node
        integer(kind=8)           :: nbCmp = 0
    end type ROM_DS_LineicNumb
!
! - Parameters for field
!
    type ROM_DS_Field

! ----- Name of field (NOM_CHAM)
        character(len=24)         :: fieldName = ' '

! ----- A real field for reference (to manipulate)
        character(len=24)         :: fieldRefe = ' '

! ----- Type of field (NOEU/ELGA)
        character(len=4)          :: fieldSupp = ' '

! ----- Model
        character(len=8)          :: model = ' '

! ----- Mesh
        character(len=8)          :: mesh = ' '

! ----- Number of equations
        integer(kind=8)                   :: nbEqua = 0

! ----- Components in the field: number and name
        character(len=8), pointer :: listCmpName(:) => null()
        integer(kind=8)                   :: nbCmpName = 0

! ----- For each dof: index of name of components (from listCmpName)
        integer(kind=8), pointer          :: equaCmpName(:) => null()

! ----- Flag if has Lagrange multipliers
        aster_logical             :: lLagr = ASTER_FALSE

! ----- Filter to apply on field
        aster_logical             :: lFilter = ASTER_FALSE

! ----- For each equation: 1 to keep this component
        integer(kind=8), pointer          :: equaFilter(:) => null()

    end type ROM_DS_Field
!
! - Datastructure for empiric base
!
    type ROM_DS_Empi
! ----- Name of empiric base to save
        character(len=8)        :: resultName = ' '
! ----- Datastructure for mode
        type(ROM_DS_Field)      :: mode
! ----- Type of reduced base
        character(len=8)        :: baseType = ' '
! ----- Direction of the linear model
        character(len=8)        :: lineicAxis = ' '
! ----- First section of the linear model
        character(len=24)       :: lineicSect = ' '
! ----- Number of modes in base
        integer(kind=8)                 :: nbMode = 0
! ----- Number of modes max
        integer(kind=8)                 :: nbModeMaxi = 0
! ----- Number of snapshots when created base
        integer(kind=8)                 :: nbSnap = 0
! ----- Datastructure for lineic base numbering
        type(ROM_DS_LineicNumb) :: lineicNume
    end type ROM_DS_Empi
!
! - Datastructure to reconstruct field
!
    type ROM_DS_FieldBuild

! ----- Complete field to write
        type(ROM_DS_Field)         :: fieldDom

! ----- Reduced field to read
        type(ROM_DS_Field)         :: fieldRom

! ----- Base to use
        type(ROM_DS_Empi)          :: base

! ----- Flag when field is solution of a linear system (primal variable)
        aster_logical              :: lLinearSolve = ASTER_FALSE

! ----- Name of operation
        character(len=24)          :: operation = ' '

! ----- RID_Total = RID_Trunc + RID_Interface

! ----- Flag to truncate RID
        aster_logical              :: lRIDTrunc = ASTER_FALSE

! ----- Name of GROUP_NO RID_Interface (when lRIDTrunc = .true.)
        character(len=24)          :: grNodeRIDInterface = ' '

! ----- Number of equations in RID (final: complete or truncated)
        integer(kind=8)                    :: nbEquaRID = 0

! ----- Access to equations in complete RID  (when lRIDTrunc = .true.)
        integer(kind=8)                    :: nbEquaRIDTotal = 0
        integer(kind=8), pointer           :: equaRIDTotal(:) => null()

! ----- Access to equation in truncated RID
        integer(kind=8)                    :: nbEquaRIDTrunc = 0
        integer(kind=8), pointer           :: equaRIDTrunc(:) => null()

! ----- [PHI] matrix on RID (size: nbEqua*nbMode)
        real(kind=8), pointer      :: matrPhi(:) => null()

! ----- [PHI] matrix on RID (size: nbEquaRID*nbMode)
        real(kind=8), pointer      :: matrPhiRID(:) => null()

! ----- Matrix of reduced coordinates for all numbering store (size: nbStore * nbMode)
        real(kind=8), pointer      :: reduMatr(:) => null()

! ----- Reconstructed field (on all domain) for all numbering store (size: nbStore * nbEqua)
        real(kind=8), pointer      :: fieldTransientVale(:) => null()

    end type ROM_DS_FieldBuild
!
! - Parameters for REST_REDUIT_COMPLET operator
!
    type ROM_DS_ParaRRC
! ----- Mesh
        character(len=8)                 :: mesh = ' '

! ----- Input result datastructure (ROM)
        type(ROM_DS_Result)              :: resultRom

! ----- Output result datastructure (DOM)
        type(ROM_DS_Result)              :: resultDom

! ----- Model for reduced model
        character(len=8)                 :: modelRom = ' '

! ----- Model for complete model
        character(len=8)                 :: modelDom = ' '

! ----- Table in result datastructure
        type(ROM_DS_TablReduCoor)        :: tablReduCoor

! ----- List of fields to reconstruct
        integer(kind=8)                          :: nbFieldBuild = 0
        character(len=24), pointer       :: fieldName(:) => null()
        type(ROM_DS_FieldBuild), pointer :: fieldBuild(:) => null()

    end type ROM_DS_ParaRRC
!
! - Parameters for definition of multiparametric reduced problem - Evaluation
!
    type ROM_DS_EvalCoef
        integer(kind=8)                     :: nb_para = 0
        real(kind=8)                :: para_vale(5) = 0.d0
        character(len=16)           :: para_name(5) = ' '
    end type ROM_DS_EvalCoef
!
! - Parameters for definition of multiparametric reduced problem - Variations
!
    type ROM_DS_VariPara
        integer(kind=8)                     :: nb_vale_para = 0
        real(kind=8), pointer       :: para_vale(:) => null()
        character(len=16)           :: para_name = ' '
        real(kind=8)                :: para_init = 0.d0
    end type ROM_DS_VariPara
!
! - Parameters for definition of multiparametric reduced problem - Coefficients
!
    type ROM_DS_MultiCoef
! ----- Coefficient is function
        aster_logical               :: l_func = ASTER_FALSE
! ----- Coefficient is constant
        aster_logical               :: l_cste = ASTER_FALSE
! ----- Coefficient is complex
        aster_logical               :: l_cplx = ASTER_FALSE
! ----- Coefficient is real
        aster_logical               :: l_real = ASTER_FALSE
! ----- Value of coefficient if is complex and constant
        complex(kind=8)             :: coef_cste_cplx = dcmplx(0.d0, 0.d0)
! ----- Value of coefficient if is real and constant
        real(kind=8)                :: coef_cste_real = 0.d0
! ----- Value of coefficient if is function
        character(len=8)            :: func_name = ' '
! ----- Value of coefficient if is complex: need evaluation
        complex(kind=8), pointer    :: coef_cplx(:) => null()
! ----- Value of coefficient if is real: need evaluation
        real(kind=8), pointer       :: coef_real(:) => null()
    end type ROM_DS_MultiCoef
!
! - Parameters for definition of multiparametric reduced problem
!
    type ROM_DS_MultiPara
! ----- Type of system to solve
        character(len=1)                 :: syst_type = ' '
! ----- List of matrices for system
        integer(kind=8)                          :: nb_matr = 0
        character(len=8), pointer        :: matr_name(:) => null()
        character(len=8), pointer        :: matr_type(:) => null()
        type(ROM_DS_MultiCoef), pointer  :: matr_coef(:) => null()
! ----- List of vectors for system
        integer(kind=8)                          :: nb_vect = 0
        character(len=8), pointer        :: vect_name(:) => null()
        character(len=8), pointer        :: vect_type(:) => null()
        type(ROM_DS_MultiCoef), pointer  :: vect_coef(:) => null()
! ----- Products matrix by current mode
        character(len=24), pointer       :: matr_mode_curr(:) => null()
! ----- Products matrix by mode
        character(len=24), pointer       :: prod_matr_mode(:) => null()
! ----- Reduced Vector
        character(len=24), pointer       :: vect_redu(:) => null()
! ----- Reduced matrix
        character(len=24), pointer       :: matr_redu(:) => null()
! ----- Variation of coefficients: number (by mode)
        integer(kind=8)                          :: nb_vari_coef = 0
! ----- Variation of coefficients: type (DIRECT, ALEATOIRE, etc. )
        character(len=24)                :: type_vari_coef = ' '
! ----- Variation of coefficients: by parameter
        integer(kind=8)                          :: nb_vari_para = 0
        type(ROM_DS_VariPara), pointer   :: vari_para(:) => null()
! ----- Evaluation of coefficients
        type(ROM_DS_EvalCoef)            :: evalcoef
! ----- Reference field
        type(ROM_DS_Field)               :: field
    end type ROM_DS_MultiPara
!
! - Parameters to solve systems
!
    type ROM_DS_Solve
        character(len=1)         :: syst_type = ' '
        character(len=1)         :: syst_matr_type = ' '
        character(len=1)         :: syst_2mbr_type = ' '
        character(len=19)        :: syst_matr = ' '
        character(len=19)        :: syst_2mbr = ' '
        character(len=19)        :: syst_solu = ' '
        character(len=19)        :: vect_zero = ' '
        integer(kind=8)                  :: syst_size = 0
    end type ROM_DS_Solve
!
! - Parameters for DEFI_BASE_REDUITE operator (POD)
!
    type ROM_DS_ParaDBR_POD
! ----- Name of result to read (high fidelity)
        character(len=8)           :: resultDomName = ' '

! ----- Result to read (high fidelity)
        type(ROM_DS_Result)        :: resultDom

! ----- Name of field to read (NOM_CHAM)
        character(len=24)          :: fieldName = ' '

! ----- Field to read (high fidelity)
        type(ROM_DS_Field)         :: field
        integer(kind=8)                    :: nbCmpToFilter = 0
        character(len=8), pointer  :: cmpToFilter(:) => null()
        integer(kind=8)                    :: nbVariToFilter = 0
        character(len=16), pointer :: variToFilter(:) => null()

! ----- Type of reduced base
        character(len=8)           :: baseType = ' '

! ----- Direction of the linear model
        character(len=8)           :: lineicAxis = ' '

! ----- First section of the linear model
        character(len=24)          :: lineicSect = ' '

! ----- Tolerance for SVD
        real(kind=8)               :: toleSVD = 0.d0

! ----- Tolerance for incremental POD
        real(kind=8)               :: toleIncr = 0.d0

! ----- Table for reduced coordinates
        type(ROM_DS_TablReduCoor)  :: tablReduCoor

! ----- Maximum number of modes
        integer(kind=8)                    :: nbModeMaxi = 0

! ----- Datastructure for snapshot selection
        type(ROM_DS_Snap)          :: snap

    end type ROM_DS_ParaDBR_POD
!
! - Algorithm Greedy
!
    type ROM_DS_AlgoGreedy
! ----- List of reduced components
        character(len=24)       :: coef_redu = ' '
! ----- For residual
        character(len=1)        :: resi_type = ' '
        character(len=24)       :: resi_vect = ' '
        real(kind=8), pointer   :: resi_norm(:) => null()
        real(kind=8)            :: resi_refe = 0.d0
! ----- To solve complete system
        type(ROM_DS_Solve)      :: solveROM
! ----- To solve reduced system
        type(ROM_DS_Solve)      :: solveDOM
! ----- Index of components FSI transient problem
        integer(kind=8)                 :: nume_pres = 0
        integer(kind=8)                 :: nume_phi = 0
    end type ROM_DS_AlgoGreedy
!
! - Parameters for DEFI_BASE_REDUITE operator (GREEDY)
!
    type ROM_DS_ParaDBR_Greedy
! ----- Datastructure for solver's parameters
        character(len=19)       :: solver = ' '

! ----- Datastructure for multiparametric problem
        type(ROM_DS_MultiPara)  :: multiPara

! ----- Maximum number of modes
        integer(kind=8)                 :: nbModeMaxi = 0

! ----- Flag to orthogonalize the basis
        aster_logical           :: lOrthoBase = ASTER_FALSE

! ----- Flag to stabilize the basis for FSI transient problem
        aster_logical           :: lStabFSI = ASTER_FALSE

! ----- Tolerance for greedy algorithm
        real(kind=8)            :: toleGreedy = 0.d0

! ----- Datastructure for greedy algorithm
        type(ROM_DS_AlgoGreedy) :: algoGreedy

    end type ROM_DS_ParaDBR_Greedy
!
! - Parameters for DEFI_BASE_REDUITE operator (TRUNCATION)
!
    type ROM_DS_ParaDBR_Trunc
! ----- Name of base to truncate
        character(len=8)        :: baseInitName = ' '

! ----- Base to truncate
        type(ROM_DS_Empi)       :: baseinit

! ----- Model for truncation
        character(len=8)        :: modelRom = ' '

! ----- List of equations into RID
        integer(kind=8), pointer        :: equaRom(:) => null()

! ----- Profile of nodal field
        character(len=24)       :: profChnoRom = ' '

! ----- Number of equation for RID
        integer(kind=8)                 :: nbEquaRom = 0

! ----- Index of GRANDEUR
        integer(kind=8)                 :: physNume = 0
    end type ROM_DS_ParaDBR_Trunc
!
! - Parameters for DEFI_BASE_REDUITE operator (ORTHO)
!
    type ROM_DS_ParaDBR_Ortho
! ----- Name of base to orthogonalize
        character(len=8)        :: baseInitName = ' '

! ----- Base to orthogonalize
        type(ROM_DS_Empi)       :: baseinit

! ----- Parameter for KAHAN-PARLETT algorithm
        real(kind=8)            :: alpha = 0.d0

    end type ROM_DS_ParaDBR_Ortho
!
! - Parameters for DEFI_BASE_REDUITE operator
!
    type ROM_DS_ParaDBR
! ----- Type of operation (POD, POD_INCR, GREEDY, ...)
        character(len=16)           :: operation = ' '

! ----- Name of empiric base to save
        type(ROM_DS_Result)         :: resultOut

! ----- Parameters for POD/POD_INCR method
        type(ROM_DS_ParaDBR_POD)    :: paraPod

! ----- Parameters for RB method
        type(ROM_DS_ParaDBR_Greedy) :: paraGreedy

! ----- Parameters for truncation method
        type(ROM_DS_ParaDBR_Trunc)  :: paraTrunc

! ----- Parameters for orthogonalization method
        type(ROM_DS_ParaDBR_Ortho)  :: paraOrtho

! ----- Datastructure for modes
        type(ROM_DS_Empi)           :: base

! ----- If operator is "reuse"
        aster_logical               :: lReuse = ASTER_FALSE
    end type ROM_DS_ParaDBR
!
! - Parameters for DEFI_DOMAINE_REDUIT operator
!
    type ROM_DS_ParaDDR
! ----- Mesh
        character(len=8)  :: mesh = ' '
! ----- Datastructure for empiric modes (primal)
        type(ROM_DS_Empi) :: ds_empi_prim
! ----- Datastructure for empiric modes (dual)
        type(ROM_DS_Empi) :: ds_empi_dual
! ----- Name of group of elements for RID
        character(len=24) :: grelem_rid = ' '
! ----- Number of layers in the construction of RID
        integer(kind=8)           :: nb_layer_rid = 0
! ----- The RID in a restricted domain?
        aster_logical     :: l_rid_maxi = ASTER_FALSE
! ----- List of elements restricted
        integer(kind=8), pointer  :: v_rid_maxi(:) => null()
! ----- Number of elements restricted
        integer(kind=8)           :: nb_rid_maxi = 0
! ----- Name of group of nodes for interface
        character(len=24) :: grnode_int = ' '
! ----- Flag for EF corrector?
        aster_logical     :: l_corr_ef = ASTER_FALSE
! ----- Name of group of nodes for outside area of EF corrector
        character(len=24) :: grnode_sub = ' '
! ----- Number of layer in the construction of outside area
        integer(kind=8)           :: nb_layer_sub = 0
        integer(kind=8)           :: nb_rid_mini = 0
! ----- List of nodes for minimal rid
        integer(kind=8), pointer  :: v_rid_mini(:) => null()
    end type ROM_DS_ParaDDR
!
! - Parameters for non_linear operator
!
    type ROM_DS_AlgoPara
! ----- Empiric modes
        type(ROM_DS_Empi) :: ds_empi
! ----- Pointer to list of equations for interface nodes
        integer(kind=8), pointer  :: v_equa_int(:) => null()
! ----- Pointer to list of equation for internal interface nodes
        integer(kind=8), pointer  :: v_equa_sub(:) => null()
! ----- Flag for reduced model
        aster_logical     :: l_rom = ASTER_FALSE
! ----- Flag for hyper-reduced model
        aster_logical     :: l_hrom = ASTER_FALSE
! ----- Flag for hyper-reduced model with EF correction
        aster_logical     :: l_hrom_corref = ASTER_FALSE
! ----- Phase of computation when EF correction
        character(len=24) :: phase = ' '
! ----- Name of GROUP_NO of interface
        character(len=24) :: grnode_int = ' '
! ----- Name of GROUP_NO of internal interface
        character(len=24)   :: grnode_sub = ' '
! ----- Table in result datastructure
        type(NL_DS_TableIO) :: tablResu
! ----- Object to save reduced coordinates
        character(len=24)   :: gamma = '&&ROM.COOR_REDU'
! ----- Identificator for field
        character(len=24) :: field_iden = ' '
! ----- Penalisation parameter for EF correction
        real(kind=8)      :: vale_pena = 0.d0
! ----- Pseuo error indicator for ROM
        real(kind=8)      :: eref_rom = -1.d0
    end type ROM_DS_AlgoPara
!
end module
