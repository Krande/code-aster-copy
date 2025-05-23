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

module elg_module
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use elg_context_type
    use elg_data_module
    use petsc_data_module
    use saddle_point_context_type
    use csc_matrix_type
    use matrix_conversion_module
    use constraint_module
!
    implicit none
!
    private
#include "asterf.h"
#include "asterfort/apalmc.h"
#include "asterfort/apmamc.h"
#include "asterfort/apmain.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    public ::  build_elg_context
!
contains

#ifdef ASTER_HAVE_PETSC
!
! The routine builds a elg_context that is a container
! with everything needed for the elimination of
! Lagrange multipliers. The main step is to build a nullbasis for
! the constraints matrix.
!
    subroutine build_elg_context(full_matas)
        !
        ! Dummy arguments
        character(len=19), intent(in) :: full_matas
        !
        ! Local variables
        !
        integer(kind=8) :: iret
        integer(kind=8) :: kptsc, kptscr, ifm, niv, nlag
        aster_logical :: verbose
        character(len=3) :: matd
        character(len=19) :: kbid
        type(saddlepoint_ctxt)    :: sp_ctxt
        type(elim_lagr_ctxt), pointer :: elg_ctxt => null()
        logical :: k_mat_to_free
        character(len=19) :: nomat_save
        PetscErrorCode :: ierr
        PetscReal :: tol, val
        PC :: pc
        Mat :: factor
        PetscInt :: ival, icntl
!
        k_mat_to_free = .false.
        !
        call infniv(ifm, niv)
        verbose = (niv == 2)
        ! Retrouve l'identifiant de l'objet elg_context associe
        ! a la matrice full_matas.
        ! En effet, on a appelé elg_gest_data('NOTE', full_matas, reduced_matas)
        ! dans elg_preres.
        call elg_gest_data('CHERCHE', full_matas, ' ', ' ')
        ! Alias vers cet objet
        elg_ctxt => elg_context(ke)
        call free_elg_context(elg_ctxt, keep_basis=.true.)
        ! On recherche egalement l'identifiant de la matrice PETSc
        ! associée à full_matas
        kptsc = get_mat_id(full_matas)
        !
        ! Pour l'instant on interdit les matrices distribuées
        call dismoi('MATR_DISTRIBUEE', full_matas, 'MATR_ASSE', repk=matd)
        ASSERT(matd == 'NON')
        ! Creation d'un clone PETSc de la matrice (matr_asse) definissant
        ! le systeme aster complet (avec doubles multiplicateurs de Lagrange)
        ! La matrice est deja enregistree et possede l'id kptsc
        ! On prealloue la matrice PETSc correspondante
        if (ap(kptsc) /= PETSC_NULL_MAT) then
            call MatDestroy(ap(kptsc), ierr)
        end if
        call apalmc(kptsc)
        ! On copie les valeurs de la matr_asse dans la matrice PETSc
        call apmamc(kptsc)
        ! et on assemble
        call MatAssemblyBegin(ap(kptsc), MAT_FINAL_ASSEMBLY, ierr)
        ASSERT(ierr == 0)
        call MatAssemblyEnd(ap(kptsc), MAT_FINAL_ASSEMBLY, ierr)
        ASSERT(ierr == 0)
        !
        ! Parfois (si matrice de masse ou amortissement), il faut aller chercher
        ! les relations lineaires dans une autre matrice
        !
        kptscr = kptsc
        if (elg_ctxt%k_matas /= " ") then
            if (verbose) then
                write (6, *) " C'est la matrice ", elg_ctxt%k_matas, &
                    " qui contient les relations linéaires de ", elg_ctxt%full_matas
            end if
            kptscr = get_mat_id(elg_ctxt%k_matas)
            !  S'il n'existe pas encore (ou deja plus) un clone PETSc de cette matrice
            if (kptscr == 0) then
                ! on le crée
                ! il faudra donc le détruire : on met à jour le flag k_mat_to_free
                k_mat_to_free = .true.
                call mat_record(elg_ctxt%k_matas, nosols(kptsc), kptscr)
                ASSERT(kptscr /= 0)
                ! on met à jour (temporairement) le nom de la matrice courante
                ! en effet, il est utilise dans les routines apalmc et apmamc
                ! pour recuperer les valeurs de la matr_asse que l'on est
                ! en train de cloner. On n'a pas besoin de mettre à jour
                ! le nom du nume_ddl: c'est le meme que pour la matrice sur laquelle
                ! on est en train de proceder a l'elimination des lagranges.
                nomat_save = nomat_courant
                nomat_courant = elg_ctxt%k_matas
                ! On prealloue la matrice PETSc correspondante
                call apalmc(kptscr)
                ! On copie les valeurs de la matr_asse dans la matrice PETSc
                call apmamc(kptscr)
                ! et on assemble
                call MatAssemblyBegin(ap(kptscr), MAT_FINAL_ASSEMBLY, ierr)
                ASSERT(ierr == 0)
                call MatAssemblyEnd(ap(kptscr), MAT_FINAL_ASSEMBLY, ierr)
                ASSERT(ierr == 0)
                ! On retablit le nom de la matrice courante
                nomat_courant = nomat_save
            end if
        end if
        !
        ! Construction d'un saddle_point_context
        sp_ctxt = new_saddle_point_context(full_matas, replicated_data, &
                                           ap(kptsc), ap(kptscr))
        !  et on libere le clone PETSc de full_matas
        kbid = repeat(" ", 19)
        call apmain('DETR_MAT', kptsc, [0.d0], kbid, 0_8, iret)
        ! ainsi que celui de k_matas
        if (k_mat_to_free) then
            call apmain('DETR_MAT', kptscr, [0.d0], kbid, 0_8, iret)
        end if
        ! On vérifie que la matrice des contraintes n'est pas vide
        nlag = get_num_of_constraints(sp_ctxt)
        if (nlag == 0) then
            if (elg_ctxt%k_matas /= " ") then
                call utmess('F', 'ELIMLAGR_10', sk=elg_ctxt%k_matas)
            else
                call utmess('F', 'ELIMLAGR_10', sk=full_matas)
            end if
        end if
        ! Remplissage de elg_ctxt (à l'aide des matrices définies dans le
        ! saddle_point_context)
        call MatConvert(sp_ctxt%k_mat, MATSAME, MAT_INITIAL_MATRIX, &
                        elg_ctxt%matb, ierr)
        ASSERT(ierr == 0)
        call MatConvert(sp_ctxt%c_mat, MATSAME, MAT_INITIAL_MATRIX, &
                        elg_ctxt%matc, ierr)
        ASSERT(ierr == 0)
        !
        elg_ctxt%nphys = sp_ctxt%nphys
        elg_ctxt%nlag = sp_ctxt%nlag1
        !
        ! On n'a plus besoin du saddle_point context
        call free_saddle_point_context(sp_ctxt)
        !
        if (elg_ctxt%tfinal == PETSC_NULL_MAT) then
            ! Calcul de la base du noyau (appel SuperLU)
            if (verbose) then
                write (6, *) " -- Construction de la base du noyau "
            end if
            call get_nullbasis(elg_ctxt%matc, elg_ctxt%tfinal)
        end if
        if (elg_ctxt%cct == PETSC_NULL_MAT) then
            !   -- Calcul de CCT = C * transpose(C)
            if (verbose) then
                write (6, *) " -- Construction de CCT "
            end if
            call MatMatTransposeMult(elg_ctxt%matc, elg_ctxt%matc, &
                                     MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, elg_ctxt%cct, ierr)
            ASSERT(ierr == 0)
            if (verbose) then
                write (6, *) " -- Factorisation de CCT "
            end if
        end if
        if (elg_ctxt%ksp == PETSC_NULL_KSP) then
#ifdef ASTER_PETSC_HAVE_MUMPS
            !  Create linear solver context
            call KSPCreate(PETSC_COMM_SELF, elg_ctxt%ksp, ierr)
            ASSERT(ierr == 0)
            !  Set operators. Here the matrix that defines the linear system
            !  also serves as the preconditioning matrix.
            !
            call KSPSetOperators(elg_ctxt%ksp, elg_ctxt%cct, elg_ctxt%cct, ierr)
            ASSERT(ierr == 0)
            tol = 1.d-7
            call KSPSetTolerances(elg_ctxt%ksp, tol, PETSC_DEFAULT_REAL, &
                                  PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
            ASSERT(ierr == 0)
            call KSPSetType(elg_ctxt%ksp, KSPPREONLY, ierr)
            ASSERT(ierr == 0)
            call KSPGetPC(elg_ctxt%ksp, pc, ierr)
            ASSERT(ierr == 0)
            call PCSetType(pc, PCCHOLESKY, ierr)
            ASSERT(ierr == 0)
            call PCFactorSetMatSolverType(pc, MATSOLVERMUMPS, ierr)
            ASSERT(ierr == 0)
            call PCFactorSetUpMatSolverType(pc, ierr)
            ASSERT(ierr == 0)
            call PCFactorGetMatrix(pc, factor, ierr)
            ASSERT(ierr == 0)
            !     sequential ordering
            icntl = 7
            ival = 2
            call MatMumpsSetIcntl(factor, icntl, ival, ierr)
            ASSERT(ierr == 0)
            !     threshhold for row pivot detection
            icntl = 24
            ival = 1
            call MatMumpsSetIcntl(factor, icntl, ival, ierr)
            ASSERT(ierr == 0)
            icntl = 3
            val = 1.d-6
            call MatMumpsSetCntl(factor, icntl, val, ierr)
            ASSERT(ierr == 0)
            call KSPSetUp(elg_ctxt%ksp, ierr)
            ASSERT(ierr == 0)
#else
            ! should not pass here
            call utmess('F', 'FERMETUR_4', sk='MUMPS')
#endif
        end if
        !
        ! Projection T'*(MatB*T)
        !
        if (verbose) then
            write (6, *) " -- Projection du problème sur la base du noyau "
        end if
        call MatPtAP(elg_ctxt%matb, elg_ctxt%tfinal, MAT_INITIAL_MATRIX, 1.d0, &
                     elg_ctxt%kproj, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine build_elg_context
!
! On entry, c_mat is a PETSc MATSEQ matrix containing the
! the constraints imposed on the physical dofs.
! On exit, z_mat is a new PETSc MATSEQ matrix containing
! a basis of the nullspace of c_mat
!
    subroutine get_nullbasis(c_mat, z_mat)
        ! Dummy arguments
        Mat, intent(in)  :: c_mat
        Mat, intent(out) :: z_mat
        ! Local variables
        aster_logical :: verbose, debug
        PetscErrorCode :: ierr
        integer(kind=8) :: ifm, niv
        real(kind=8)::  t0, t1
        type(csc_matrix) :: ct_csc, b_csc, z_csc
        Mat :: cz_mat
        PetscScalar :: c_nrm, cz_nrm, z_nrm
        PetscInt, dimension(:), pointer :: ia => null(), ja => null()
        PetscScalar, dimension(:), pointer :: val => null()
        PetscInt :: nn
        PetscInt :: mm_z, nn_z
!
        call infniv(ifm, niv)
        verbose = (niv == 2)
        debug = .false.
!
        call CPU_time(t0)
!   Conversion de c_mat au format CSR
        call matseq2csr(c_mat, nn, ia, ja, val)
!   Les tableaux obtenus définissent C^T au format CSC
        call define_csc_matrix_from_array("CT", one_based, to_aster_int(nn), ja, &
                                          ia, val, ct_csc)
!   Si la matrice C.T n'est pas de rang maximal (existence de contraintes
!   redondantes ou dépendantes), on commence par se ramener à une matrice
!   de rang maximal B
        call get_columnspace_basis(ct_csc, b_csc)
!   On n'a plus besoin de ct_csc
        call free_csc_matrix(ct_csc)
!   On calcule une base du noyau de B stockée dans Z
        call get_nullbasis_trans(b_csc, z_csc)
!   On change de format et on stocke la base dans une matrice PETSc (Tfinal)
        call csc2matseq(z_csc, z_mat)
!   Vérification
        call MatGetSize(z_mat, mm_z, nn_z, ierr)
        ASSERT(ierr == 0)
        call MatMatMult(c_mat, z_mat, MAT_INITIAL_MATRIX, &
                        PETSC_DEFAULT_REAL, cz_mat, ierr)
        ASSERT(ierr == 0)
        call MatNorm(cz_mat, NORM_FROBENIUS, cz_nrm, ierr)
        ASSERT(ierr == 0)
        call MatNorm(c_mat, NORM_FROBENIUS, c_nrm, ierr)
        ASSERT(ierr == 0)
        call MatNorm(z_mat, NORM_FROBENIUS, z_nrm, ierr)
        ASSERT(ierr == 0)

        call CPU_time(t1)
        if (verbose) then
            call utmess('I', 'ELIMLAGR_13', si=to_aster_int(nn_z), &
                        nr=4_8, valr=(/c_nrm, z_nrm, cz_nrm, t1-t0/))
        end if
        if (debug) then
            print *, "ELG Norme de C                             : ", c_nrm
            print *, "ELG Norme de Z                             : ", z_nrm
            print *, "ELG Norme de CZ                            : ", cz_nrm
            print *, "ELG Norme de CZ/( Norme de C * Norme de Z ): ", cz_nrm/(c_nrm*z_nrm)
        end if
        ASSERT(cz_nrm/(c_nrm*z_nrm) < 1.d-2)
! Libération de la mémoire
        call free_csc_matrix(b_csc)
        call free_csc_matrix(z_csc)
        call MatDestroy(cz_mat, ierr)
        ASSERT(ierr == 0)
!
    end subroutine get_nullbasis
!
#else
    subroutine build_elg_context(full_matas)
        character(len=19), intent(in) :: full_matas
        character(len=19)  :: kbid
        kbid = full_matas
    end subroutine build_elg_context
#endif
end module elg_module
