! --------------------------------------------------------------------
! Copyright (C) 2016 - 2025 - EDF R&D - www.code-aster.org
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
! A saddle_point_context object is a container used to manage
! a saddle point linear system
! ( k_mat c_mat^T ) (x_1) = (b_1)
! ( c_mat 0       ) (x_2)   (b_2)
! For simplicity the saddle point linear system is embedded
! in the larger system
! ( k_mat c_mat^T 0 ) (x_1) = (b_1)
! ( c_mat 0       0 ) (x_2)   (b_2)
! ( 0     0       Id) (x_3)   (b_3)
!
! It contains
! - Index Sets necessary to extract data from the global (double Lagrange) Aster system
! - matrix data ( k_mat, c_mat, d_mat )
! - vector workspace ( x_1, x_2, x_3, b_1, b_2, b_3 )
!
module saddle_point_context_type
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use matrasse_module
!
    implicit none
!
    private
#include "asterf.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/conlag.h"
!
!
    type, public :: saddlepoint_ctxt
#ifdef ASTER_HAVE_PETSC
        !
        ! Double Lagrange scaling coefficient (from .conl)
        real(kind=8) :: alpha
        !
        ! Two modes are available : distributed/replicated
        ! distributed -> the full matrix is distributed among the
        ! available processors
        ! replicated -> the full matrix is replicated on all the
        ! processors which participate to the computation
        integer(kind=8) :: data_model
        !
        ! Index Sets section
        ! ==================
        ! Flag : true if Index Sets (IS) have been setup
        aster_logical :: is_setup = .false.
        ! Index set of physical degrees of freedom
        ! in PETSc matrix (of the whole double Lagrange system)
        IS :: is_phys
        ! Number of physical degrees of freedom
        integer(kind=8) :: nphys
        !
        ! Index set of Lagrange 1 multipliers
        ! in PETSc matrix (of the whole double Lagrange system)
        IS :: is_lag1
        ! Number of Lagrange 1 multipliers
        integer(kind=8) :: nlag1
        ! Index set of Lagrange 2 multipliers
        ! in PETSc matrix (of the whole double Lagrange system)
        IS :: is_lag2
        ! Number of Lagrange 2 multipliers
        integer(kind=8) :: nlag2
        !
        ! Matrix Data section
        ! ============
        aster_logical :: data_setup = .false.
        ! Stiffness Matrix
        Mat :: k_mat
        ! Constraint Matrix
        ! c_mat is obtained from is_lag1 x is_phys
        ! whereas d_mat is obtained from is_lag2 x is_phys
        ! They are identical on a single proc, but they have different parallel layout
        Mat :: c_mat, d_mat
        !
        !
        ! Work Space
        ! ==========
        aster_logical :: work_setup = .false.
        Vec :: x1, x2, x3, xtmp
        Vec :: y1, y2, y3
        ! Scatter tools
        ! =============
        aster_logical :: scatter_setup = .false.
        VecScatter :: scatter_to_phys, scatter_to_lag1, scatter_to_lag2
        !
#else
        integer(kind=8) :: idummy
#endif
    end type saddlepoint_ctxt
!
    integer(kind=8), parameter, public :: distributed_data = 0, replicated_data = 1
!
    public :: new_saddle_point_context, free_saddle_point_context, get_num_of_constraints
!
#ifdef ASTER_HAVE_PETSC
!
    PetscErrorCode  :: ierr
!
contains
!
    function new_saddle_point_context(full_matas, data_model, a_mat, ak_mat) result(ctxt)
        !
        ! Dummy arguments
        !
        character(len=19), intent(in)                :: full_matas
        integer(kind=8), intent(in)                          :: data_model
        Mat, intent(in)                              :: a_mat
        Mat, intent(in), optional                    :: ak_mat
        type(saddlepoint_ctxt)                       :: ctxt
        !
        ! Local variables
        !
        real(kind=8) :: un_sur_alpha
        !
        ! Modèle de données en parallèle
        ASSERT((data_model == distributed_data) .or. (data_model == replicated_data))
        ctxt%data_model = data_model
        ! On récupère l'inverse du paramètre alpha
        call conlag(full_matas, un_sur_alpha)
        ctxt%alpha = 1.d0/un_sur_alpha
        !
        call set_is(full_matas, ctxt)
        ASSERT(ctxt%is_setup)
        if (present(ak_mat)) then
            ! Les relations lineaires sont stockees dans
            ! la matrice ak_mat
            call set_matrix_data(a_mat, ak_mat, ctxt)
        else
            call set_matrix_data(a_mat, a_mat, ctxt)
        end if
        ASSERT(ctxt%data_setup)
        call set_workspace(ctxt)
        ASSERT(ctxt%work_setup)
        call set_scatter(a_mat, ctxt)
        ASSERT(ctxt%scatter_setup)
        !
    end function new_saddle_point_context
!
    function get_num_of_constraints(ctxt) result(nlag)
!
!   Dummy arguments
        type(saddlepoint_ctxt), intent(in)         :: ctxt
        integer(kind=8)                                    :: nlag
        ASSERT(ctxt%nlag1 == ctxt%nlag2)
        nlag = ctxt%nlag1
!
    end function get_num_of_constraints
!
! This routine initializes the Index sets section of the context object
!
    subroutine set_is(matas, ctxt)
        !
        ! Dummy arguments
        !
        character(len=19), intent(in)                  :: matas
        type(saddlepoint_ctxt), intent(inout)          :: ctxt
        !
        ! Local variables
        PetscInt :: nphys, nlag1, nlag2
        !
        ASSERT(.not. ctxt%is_setup)
        !
        ! Degrés de liberté physiques :
        ! =============================
        ctxt%is_phys = build_is_for_type_of_dof(matas, physical_dof, &
                                                ctxt%data_model)
        call ISGetSize(ctxt%is_phys, nphys, ierr)
        ASSERT(ierr == 0)
        ctxt%nphys = to_aster_int(nphys)
        !
        ! Lagrange 1 :
        ! ============
        ctxt%is_lag1 = build_is_for_type_of_dof(matas, lagrange1_dof, &
                                                ctxt%data_model)
        call ISGetSize(ctxt%is_lag1, nlag1, ierr)
        ASSERT(ierr == 0)
        ctxt%nlag1 = to_aster_int(nlag1)
        !
        ! Lagrange 2 :
        ! ============
        ctxt%is_lag2 = build_is_for_type_of_dof(matas, lagrange2_dof, &
                                                ctxt%data_model)
        call ISGetSize(ctxt%is_lag2, nlag2, ierr)
        ASSERT(ierr == 0)
        ctxt%nlag2 = to_aster_int(nlag2)
        ! Cette étape est terminée
        ctxt%is_setup = .true.
        !
    end subroutine set_is
!
! This routine initializes the data section of the context object
! The IS section is supposed to be OK
! a_mat contains the (double Lagrange) matrix from which k_mat and
! c_mat are extracted.
! The way the submatrix is extracted depends on the parallel data
! model : distributed_data/replicated_data
!
    subroutine set_matrix_data(a_mat, ak_mat, ctxt)
        !
        ! Dummy arguments
        !
        Mat, intent(in)                                :: a_mat
        Mat, intent(in)                                :: ak_mat
        type(saddlepoint_ctxt), intent(inout)          :: ctxt
        !
        ! Local variables
        !
        Mat, dimension(3) :: submat
        PetscInt :: nsub
        !
        ASSERT((ctxt%data_model == distributed_data) .or. (ctxt%data_model == replicated_data))
        ASSERT(ctxt%is_setup)
        ASSERT(.not. ctxt%data_setup)
        ! Définition du bloc k_mat : nouvelle matrice PETSc contenant
        ! les interactions des ddls "physiques" du modèle
        if (ctxt%data_model == distributed_data) then
            call MatCreateSubMatrix(a_mat, ctxt%is_phys, ctxt%is_phys, &
                                    MAT_INITIAL_MATRIX, ctxt%k_mat, ierr)
            ASSERT(ierr == 0)
        else if (ctxt%data_model == replicated_data) then
            nsub = to_petsc_int(1)
            call MatCreateSubMatrices(a_mat, nsub, &
                                      ctxt%is_phys, ctxt%is_phys, MAT_INITIAL_MATRIX, &
                                      submat(1), ierr)
            ASSERT(ierr == 0)
            call MatConvert(submat(1), MATSAME, MAT_INITIAL_MATRIX, &
                            ctxt%k_mat, ierr)
            ASSERT(ierr == 0)
        end if
        !
        ! Définition du bloc c_mat des contraintes: nouvelle matrice PETSc
        ! contenant les interactions des ddls lagrange1 (lignes) avec les ddls
        ! physiques (colonnes)
        ! Attention ! il faut peut-être utiliser la matrice de rigidité
        !
        if (ctxt%data_model == distributed_data) then
            call MatCreateSubMatrix(ak_mat, ctxt%is_lag1, ctxt%is_phys, &
                                    MAT_INITIAL_MATRIX, ctxt%c_mat, ierr)
            ASSERT(ierr == 0)
        else if (ctxt%data_model == replicated_data) then
            nsub = to_petsc_int(1)
            call MatCreateSubMatrices(ak_mat, nsub, &
                                      ctxt%is_lag1, ctxt%is_phys, MAT_INITIAL_MATRIX, &
                                      submat(2), ierr)
            ASSERT(ierr == 0)
            call MatConvert(submat(2), MATSAME, MAT_INITIAL_MATRIX, &
                            ctxt%c_mat, ierr)
            ASSERT(ierr == 0)
        end if
        !
        if (ctxt%data_model == distributed_data) then
            call MatCreateSubMatrix(ak_mat, ctxt%is_lag2, ctxt%is_phys, &
                                    MAT_INITIAL_MATRIX, ctxt%d_mat, ierr)
            ASSERT(ierr == 0)
        else if (ctxt%data_model == replicated_data) then
            nsub = to_petsc_int(1)
            call MatCreateSubMatrices(ak_mat, nsub, &
                                      ctxt%is_lag2, ctxt%is_phys, MAT_INITIAL_MATRIX, &
                                      submat(3), ierr)
            ASSERT(ierr == 0)
            call MatConvert(submat(3), MATSAME, MAT_INITIAL_MATRIX, &
                            ctxt%d_mat, ierr)
            ASSERT(ierr == 0)
        end if
        !
        if (ctxt%data_model == replicated_data) then
            call MatDestroy(submat(1), ierr)
            ASSERT(ierr == 0)
            call MatDestroy(submat(2), ierr)
            ASSERT(ierr == 0)
            call MatDestroy(submat(3), ierr)
            ASSERT(ierr == 0)
        end if
        ctxt%data_setup = .true.
        !
    end subroutine set_matrix_data
!
! This routine initializes the workspace section of the context object
!
    subroutine set_workspace(ctxt)
        !
        ! Dummy arguments
        !
        type(saddlepoint_ctxt), intent(inout) :: ctxt
        !
        ASSERT(ctxt%data_setup)
        ASSERT(.not. ctxt%work_setup)
        !
        call MatCreateVecs(ctxt%k_mat, ctxt%x1, ctxt%y1, ierr)
        ASSERT(ierr == 0)
        call VecDuplicate(ctxt%x1, ctxt%xtmp, ierr)
        ASSERT(ierr == 0)
        call MatCreateVecs(ctxt%c_mat, PETSC_NULL_VEC, ctxt%x2, ierr)
        ASSERT(ierr == 0)
        call VecDuplicate(ctxt%x2, ctxt%y2, ierr)
        ASSERT(ierr == 0)
        call MatCreateVecs(ctxt%d_mat, PETSC_NULL_VEC, ctxt%x3, ierr)
        ASSERT(ierr == 0)
        call VecDuplicate(ctxt%x3, ctxt%y3, ierr)
        ASSERT(ierr == 0)
        !
        ctxt%work_setup = .true.
        ! On n'a plus besoin de d_mat
        call MatDestroy(ctxt%d_mat, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine set_workspace
!
    subroutine set_scatter(a_mat, ctxt)
        !
        ! Dummy arguments
        !
        Mat, intent(in)                                :: a_mat
        type(saddlepoint_ctxt), intent(inout) :: ctxt
        !
        ! Local Variables
        Vec :: x
        !
        ASSERT(ctxt%is_setup)
        ASSERT(ctxt%work_setup)
        ASSERT(.not. ctxt%scatter_setup)
        !
        call MatCreateVecs(a_mat, x, PETSC_NULL_VEC, ierr)
        ASSERT(ierr == 0)
        call VecScatterCreate(ctxt%x1, PETSC_NULL_IS, x, ctxt%is_phys, ctxt%scatter_to_phys, ierr)
        ASSERT(ierr == 0)
        call VecScatterCreate(ctxt%x2, PETSC_NULL_IS, x, ctxt%is_lag1, ctxt%scatter_to_lag1, ierr)
        ASSERT(ierr == 0)
        call VecScatterCreate(ctxt%x3, PETSC_NULL_IS, x, ctxt%is_lag2, ctxt%scatter_to_lag2, ierr)
        ASSERT(ierr == 0)
        !
        call VecDestroy(x, ierr)
        ASSERT(ierr == 0)
        !
        ctxt%scatter_setup = .true.
        !
    end subroutine set_scatter
!
    subroutine free_saddle_point_context(ctxt)
        !
        ! Dummy argument
        !
        type(saddlepoint_ctxt), intent(inout) :: ctxt
        !
        call ISDestroy(ctxt%is_phys, ierr)
        ASSERT(ierr == 0)
        call ISDestroy(ctxt%is_lag1, ierr)
        ASSERT(ierr == 0)
        call ISDestroy(ctxt%is_lag2, ierr)
        ASSERT(ierr == 0)
        call MatDestroy(ctxt%k_mat, ierr)
        ASSERT(ierr == 0)
        call MatDestroy(ctxt%c_mat, ierr)
        ASSERT(ierr == 0)
        if (ctxt%d_mat /= PETSC_NULL_MAT) then
            call MatDestroy(ctxt%d_mat, ierr)
            ASSERT(ierr == 0)
        end if
        call VecDestroy(ctxt%x1, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%x2, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%x3, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%xtmp, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%y1, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%y2, ierr)
        ASSERT(ierr == 0)
        call VecDestroy(ctxt%y3, ierr)
        ASSERT(ierr == 0)
        call VecScatterDestroy(ctxt%scatter_to_phys, ierr)
        ASSERT(ierr == 0)
        call VecScatterDestroy(ctxt%scatter_to_lag1, ierr)
        ASSERT(ierr == 0)
        call VecScatterDestroy(ctxt%scatter_to_lag2, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine free_saddle_point_context
!
!
    function build_is_for_type_of_dof(matas, type_of_dof, data_model) result(is)
        ! Dummy arguments
        character(len=19)   :: matas
        integer(kind=8), intent(in) :: type_of_dof
        integer(kind=8), intent(in) :: data_model
        IS                  :: is
        ! Local variables
        mpi_int :: mpicomm, rank, nbproc
        integer(kind=8) ::  pass, ii, jerr
        PetscInt :: istart, iend, ndof, my_ndof
        integer(kind=8), dimension(:), pointer :: idof => null()
        PetscInt, dimension(:), allocatable    :: my_idof
        Vec :: vtmp
        !
        ASSERT((data_model == distributed_data) .or. (data_model == replicated_data))
        ! Récupération du communicateur MPI
        call asmpi_comm('GET', mpicomm)
        call asmpi_info(rank=rank, size=nbproc)
        !
        ! Détermination des indices (Fortran) dans la matrice aster
        idof => get_indices_of_dofs(type_of_dof, matas)
        ndof = size(idof)
        ! Le vecteur d'indices idof contient les indices globaux des ddls
        ! physiques. Il a la même taille et contient les mêmes
        ! valeurs sur tous les processeurs participant au calcul.
        !
        ! Passage indices Fortran -> Indices C
        idof(:) = idof(:)-1
        !
        ! On distingue à présent deux types de distributions de données
        ! - distributed_data : on veut construire des matrices (rigidité, masse, contraintes ...)
        !                     de type MPIAIJ, réparties sur tous les processeurs.
        ! - replicated_data : on veut construire des matrices (rigidité, masse, contraintes ...)
        !                     de type MATSEQ, répliquées sur tous les processeurs.
        ! Détermination de istart et iend
        if (data_model == replicated_data) then
            istart = 0
            iend = ndof
        else
            ! Sur plusieurs procs
            ! l'index set retourné
            ! par cette fonction a la même répartition
            ! parallèle qu'un vecteur PETSc de taille ndof
            call VecCreateMPI(mpicomm, PETSC_DECIDE, ndof, vtmp, ierr)
            ASSERT(ierr == 0)
            call VecGetOwnerShipRange(vtmp, istart, iend, ierr)
            ASSERT(ierr == 0)
            call VecDestroy(vtmp, ierr)
            ASSERT(ierr == 0)
        end if
        ! Sur chaque processeur, l'IS  possède les degrés de liberté
        ! listés dans my_idof
        do pass = 1, 2
            my_ndof = 0
            do ii = istart, iend-1
                my_ndof = my_ndof+1
                if (pass == 2) then
                    my_idof(my_ndof) = to_petsc_int(idof(ii+1))
                end if
            end do
            if (pass == 1) then
                allocate (my_idof(my_ndof), stat=jerr)
                ASSERT(jerr == 0)
            end if
        end do
        !
        ! Construction d'un Index Set (IS) PETSc à partir du vecteur d'indices C
        call ISCreateGeneral(mpicomm, my_ndof, my_idof, &
                             PETSC_COPY_VALUES, is, ierr)
        ASSERT(ierr == 0)
        if (associated(idof)) then
            deallocate (idof)
        end if
        !
        if (allocated(my_idof)) then
            deallocate (my_idof)
        end if
        !
    end function build_is_for_type_of_dof
!
#else
! Si on ne compile pas avec PETSc, il faut quand même définir les
! interfaces des routines publiques
contains
!
    function new_saddle_point_context(matas, a_mat) result(ctxt)
        !
        character(len=19), intent(in)                :: matas
        integer(kind=8), intent(in)                          :: a_mat
        type(saddlepoint_ctxt)              :: ctxt
        character(len=1) :: kdummy
        integer(kind=8) :: idummy
        kdummy = matas(1:1)
        idummy = a_mat
        ctxt%idummy = 0
    end function new_saddle_point_context
!
    subroutine free_saddle_point_context(ctxt)
        type(saddlepoint_ctxt), intent(inout) :: ctxt
        ctxt%idummy = 0
    end subroutine free_saddle_point_context
!
    function get_num_of_constraints(ctxt) result(nlag)
        type(saddlepoint_ctxt), intent(in)         :: ctxt
        integer(kind=8)                                    :: nlag
        nlag = ctxt%idummy
    end function get_num_of_constraints
!
#endif
end module saddle_point_context_type
