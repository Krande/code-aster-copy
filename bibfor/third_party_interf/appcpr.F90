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

subroutine appcpr(kptsc)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
! person_in_charge: natacha.bereux at edf.fr
    use aster_petsc_module
    use petsc_data_module
    use augmented_lagrangian_module, only: augmented_lagrangian_apply, &
                                           augmented_lagrangian_setup, augmented_lagrangian_destroy
    use lmp_module, only: lmp_apply_right, lmp_destroy
    use lmp_data_module, only: reac_lmp
    use aster_fieldsplit_module, only: mfield_setup

    implicit none
#include "asterf.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/apmams.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ldsp1.h"
#include "asterfort/ldsp2.h"
#include "asterfort/utmess.h"
#include "asterfort/PCHPDDMDumpAuxiliaryMat.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/isParallelMatrix.h"
!
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  CREATION DU PRECONDITIONNEUR PETSC (INSTANCE NUMERO KPTSC)
!  PHASE DE PRE-TRAITEMENT (PRERES)
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!----------------------------------------------------------------
!
!   VARIABLES LOCALES
    integer(kind=8) :: rang, nbproc
    integer(kind=8) :: dimgeo, dimgeo_b, niremp
    integer(kind=8) :: jnequ, jnequl
    integer(kind=8) :: nloc, neqg, ndprop, ieq, numno, icmp
    integer(kind=8) :: iret, i
    integer(kind=8) :: il, ix, iga_f, igp_f
    integer(kind=8) :: fill, reacpr
    integer(kind=8), dimension(:), pointer :: slvi => null()
    integer(kind=8), dimension(:), pointer :: prddl => null()
    integer(kind=8), dimension(:), pointer :: deeq => null()
    integer(kind=8), dimension(:), pointer :: nulg => null()
    integer(kind=8), dimension(:), pointer :: nlgp => null()
    mpi_int :: mpicomm
!
    character(len=24) :: precon
    character(len=19) :: nomat, nosolv, syme
    character(len=14) :: nonu, factor
    character(len=8) :: nomail, nbproc_str
    character(len=4) :: exilag
    character(len=3) :: matd
    character(len=800) :: myopt
    character(len=24), dimension(:), pointer :: slvk => null()
!
    real(kind=8) :: fillin, val
    real(kind=8), dimension(:), pointer :: coordo => null()
    real(kind=8), dimension(:), pointer :: slvr => null()
    !
    aster_logical :: lmd, l_parallel_matrix, lmp_is_active
!
!----------------------------------------------------------------
!   Variables PETSc
    PetscInt :: low, high, bs, nterm, nsmooth
    PetscErrorCode :: ierr
    PetscReal :: fillp
    PetscScalar :: xx_v(1)
    PetscOffset :: xx_i
    PetscInt, pointer :: nulg_ip(:) => null()

    Mat :: a, auxMat
    Vec :: coords
    KSP :: ksp
    IS  :: auxIS
    PC  :: pc
    mpi_int :: mrank, msize
    MatNullSpace :: sp
!----------------------------------------------------------------
    call jemarq()
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicomm)
!
!  -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
    ksp = kp(kptsc)
!
    call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
    call jeveuo(nosolv//'.SLVR', 'L', vr=slvr)
    call jeveuo(nosolv//'.SLVI', 'L', vi=slvi)
    precon = slvk(2)

    fillin = slvr(3)
    niremp = slvi(4)
    call dismoi('MATR_DISTRIBUEE', nomat, 'MATR_ASSE', repk=matd)
    lmd = matd == 'OUI'
    l_parallel_matrix = isParallelMatrix(nomat)
    ASSERT(.not. (lmd .and. l_parallel_matrix))
!
    lmp_is_active = slvk(6) == 'GMRES_LMP'
    if (lmp_is_active .and. precon /= 'LDLT_SP') then
        call utmess('F', 'PETSC_28')
    end if
    reacpr = slvi(6)
!
    fill = niremp
    fillp = fillin
    bs = 1
!
!   -- RECUPERE DES INFORMATIONS SUR LE MAILLAGE POUR
!   LE CALCUL DES MODES RIGIDES
    call dismoi('NOM_MAILLA', nomat, 'MATR_ASSE', repk=nomail)
    call dismoi('DIM_GEOM_B', nomail, 'MAILLAGE', repi=dimgeo_b)
    call dismoi('DIM_GEOM', nomail, 'MAILLAGE', repi=dimgeo)
!
!  -- RECUPERE LE RANG DU PROCESSUS ET LE NB DE PROCS
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
!  -- CAS PARTICULIER (LDLT_INC/SOR)
!  -- CES PC NE SONT PAS PARALLELISES
!  -- ON UTILISE DONC DES VERSIONS PAR BLOC
!   -- QUE L'ON CREERA AU MOMENT DE LA RESOLUTION (DANS APPCRS)
!  -----------------------------------------------------------
    if ((precon == 'LDLT_INC') .or. (precon == 'SOR')) then
        if (nbproc .gt. 1) then
!        EN PARALLELE, ON NE PREPARE PAS LE PRECONDITIONNEUR
!        TOUT DE SUITE CAR ON NE VEUT PAS ETRE OBLIGE
!        D'APPELER KSPSetUp
            goto 999
        end if
    end if
!
!   -- VERIFICATIONS COMPLEMENTAIRES POUR LES PRE-CONDITIONNEURS MULTIGRILLES
!   -- DEFINITION DU NOYAU UNIQUEMENT EN MODELISATION SOLIDE (2D OU 3D)
!   -------------------------------------------------------------------------
    if ((precon == 'ML') .or. (precon == 'BOOMER') .or. (precon == 'GAMG')) then
!  -- PAS DE LAGRANGE
        call dismoi('EXIS_LAGR', nomat, 'MATR_ASSE', repk=exilag, arret='C', ier=iret)
        if (iret == 0) then
            if (exilag == 'OUI') then
                call utmess('F', 'PETSC_18')
            end if
        end if
!  -- NOMBRE CONSTANT DE DDLS PAR NOEUD
        bs = tblocs(kptsc)
        if (bs .le. 0) then
            call utmess('A', 'PETSC_18')
        else
            bs = abs(bs)
        end if
    end if

    if ((precon == 'ML') .or. (precon == 'BOOMER') .or. (precon == 'GAMG')) then
!  -- CREATION DES MOUVEMENTS DE CORPS RIGIDE --
!      * VECTEUR RECEPTABLE DES COORDONNEES AVEC TAILLE DE BLOC
!
!     dimgeo = 3 signifie que le maillage est 3D et que les noeuds ne sont pas
!     tous dans le plan z=0
!     c'est une condition nécessaire au pré-calcul des modes de corps rigides
        if (dimgeo == 3) then
            call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
            call jeveuo(nonu//'.NUME.DEEQ', 'L', vi=deeq)
            call jeveuo(nomail//'.COORDO    .VALE', 'L', vr=coordo)
            !
            call VecCreate(mpicomm, coords, ierr)
            call VecSetBlockSize(coords, bs, ierr)
            ! la matrice est distribuée
            if (lmd) then
                call jeveuo(nonu//'.NUML.NEQU', 'L', jnequl)
                call jeveuo(nonu//'.NUML.PDDL', 'L', vi=prddl)
                neqg = zi(jnequ)
                nloc = zi(jnequl)
                ! Nb de ddls dont le proc courant est propriétaire (pour PETSc)
                ndprop = 0
                do il = 1, nloc
                    if (prddl(il) == rang) then
                        ndprop = ndprop+1
                    end if
                end do
                !
                call VecSetSizes(coords, to_petsc_int(ndprop), to_petsc_int(neqg), ierr)
                ! la matrice est issue d'un ParallelMesh
            else if (l_parallel_matrix) then
                call jeveuo(nonu//'.NUME.PDDL', 'L', vi=prddl)
                call jeveuo(nonu//'.NUME.NEQU', 'L', jnequl)
                nloc = zi(jnequl)
                neqg = zi(jnequl+1)
                ndprop = 0
                do il = 1, nloc
                    if (prddl(il) == rang) then
                        ndprop = ndprop+1
                    end if
                end do
                call VecSetSizes(coords, to_petsc_int(ndprop), to_petsc_int(neqg), ierr)
                ! la matrice est centralisée
            else
                neqg = zi(jnequ)
                call VecSetSizes(coords, PETSC_DECIDE, to_petsc_int(neqg), ierr)
            end if
            !
            call VecSetType(coords, VECMPI, ierr)
            ! * REMPLISSAGE DU VECTEUR
            !   coords: vecteur PETSc des coordonnées des noeuds du maillage,
            !   dans l'ordre de la numérotation PETSc des équations
            if (lmd) then
                call jeveuo(nonu//'.NUML.NULG', 'L', vi=nulg)
                call jeveuo(nonu//'.NUML.NLGP', 'L', vi=nlgp)
                do il = 1, nloc
                    ! Indice global PETSc (F) correspondant à l'indice local il
                    igp_f = nlgp(il)
                    ! Indice global Aster (F) correspondant à l'indice local il
                    iga_f = nulg(il)
                    ! Noeud auquel est associé le ddl global Aster iga_f
                    numno = deeq((iga_f-1)*2+1)
                    ! Composante (X, Y ou Z) à laquelle est associé
                    ! le ddl global Aster iga_f
                    icmp = deeq((iga_f-1)*2+2)
                    ASSERT((numno .gt. 0) .and. (icmp .gt. 0))
                    ! Valeur de la coordonnée (X,Y ou Z) icmp du noeud numno
                    val = coordo(dimgeo_b*(numno-1)+icmp)
                    ! On met à jour le vecteur PETSc des coordonnées
                    nterm = 1
                    call VecSetValues(coords, nterm, [to_petsc_int(igp_f-1)], [val], &
                                      INSERT_VALUES, ierr)
                    ASSERT(ierr == 0)
                end do
                call VecAssemblyBegin(coords, ierr)
                call VecAssemblyEnd(coords, ierr)
                ASSERT(ierr == 0)
                ! la matrice est issue d'un ParallelMesh
            else if (l_parallel_matrix) then
                call jeveuo(nonu//'.NUME.NULG', 'L', vi=nulg)
                do il = 1, nloc
                    if (prddl(il) == rang) then
                        ! Noeud auquel est associé le ddl Aster ieq
                        numno = deeq((il-1)*2+1)
                        ! Composante (X, Y ou Z) à laquelle est associé
                        ! le ddl Aster ieq
                        icmp = deeq((il-1)*2+2)
                        ASSERT((numno .gt. 0) .and. (icmp .gt. 0))
                        ! Valeur de la coordonnée (X,Y ou Z) icmp du noeud numno
                        val = coordo(dimgeo_b*(numno-1)+icmp)
                        ! On met à jour le vecteur PETSc des coordonnées
                        nterm = 1
                        ! Indice global Aster (F) correspondant à l'indice local il
                        iga_f = nulg(il)
                        call VecSetValues(coords, nterm, [to_petsc_int(iga_f)], [val], &
                                          INSERT_VALUES, ierr)
                        ASSERT(ierr == 0)
                    end if
                end do
                ! la matrice est centralisée
            else
                call VecGetOwnershipRange(coords, low, high, ierr)
                call VecGetArray(coords, xx_v, xx_i, ierr)
                ix = 0
                do ieq = low+1, high
                    ! Noeud auquel est associé le ddl Aster ieq
                    numno = deeq((ieq-1)*2+1)
                    ! Composante (X, Y ou Z) à laquelle est associé
                    ! le ddl Aster ieq
                    icmp = deeq((ieq-1)*2+2)
                    ASSERT((numno .gt. 0) .and. (icmp .gt. 0))
                    ix = ix+1
                    xx_v(xx_i+ix) = coordo(dimgeo*(numno-1)+icmp)
                end do
                !
                call VecRestoreArray(coords, xx_v, xx_i, ierr)
                ASSERT(ierr == 0)
            end if
            !
            !
            !      * CALCUL DES MODES A PARTIR DES COORDONNEES
            if (bs .le. 3) then
                call MatNullSpaceCreateRigidBody(coords, sp, ierr)
                ASSERT(ierr == 0)
                call MatSetNearNullSpace(a, sp, ierr)
                ASSERT(ierr == 0)
                call MatNullSpaceDestroy(sp, ierr)
                ASSERT(ierr == 0)
            end if
            call VecDestroy(coords, ierr)
            ASSERT(ierr == 0)
        end if
! Si dimgeo /= 3 on ne pré-calcule pas les modes de corps rigides
    end if

!
!  -- CHOIX DU PRECONDITIONNEUR :
!  ------------------------------
    call KSPGetPC(ksp, pc, ierr)
    ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    if (precon == 'LDLT_INC') then
        call PCSetType(pc, PCILU, ierr)
        ASSERT(ierr == 0)
        call PCFactorSetLevels(pc, to_petsc_int(fill), ierr)
        ASSERT(ierr == 0)
        call PCFactorSetFill(pc, fillp, ierr)
        ASSERT(ierr == 0)
        call PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL, ierr)
        ASSERT(ierr == 0)
        !-----------------------------------------------------------------------
    else if (precon .eq. 'BLOC_LAGR') then
        call KSPGetPC(ksp, pc, ierr)
        ASSERT(ierr == 0)
        call PCSetType(pc, PCSHELL, ierr)
        ASSERT(ierr == 0)
        call PCShellSetSetUp(pc, augmented_lagrangian_setup, ierr)
        ASSERT(ierr == 0)
        call PCShellSetDestroy(pc, augmented_lagrangian_destroy, ierr)
        ASSERT(ierr == 0)
!     Si LMP, on définit un préconditionneur à gauche et à droite
        if (lmp_is_active) then
            ASSERT(ierr == 0)
            call PCShellSetName(pc, "Symmetric Preconditionner: Left BLOC_LAGR, right LMP", ierr)
            ASSERT(ierr == 0)
            call PCShellSetApplySymmetricLeft(pc, augmented_lagrangian_apply, ierr)
            ASSERT(ierr == 0)
            call PCShellSetApplySymmetricRight(pc, lmp_apply_right, ierr)
            ASSERT(ierr == 0)
            call KSPSetPCSide(ksp, PC_SYMMETRIC, ierr)
            ASSERT(ierr == 0)
        else
            call PCShellSetName(pc, "BLOC_LAGR Preconditionner", ierr)
            ASSERT(ierr == 0)
            call PCShellSetApply(pc, augmented_lagrangian_apply, ierr)
            ASSERT(ierr == 0)
        end if
!-----------------------------------------------------------------------
    else if ((precon == 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
        call PCSetType(pc, PCSHELL, ierr)
        ASSERT(ierr == 0)
!     LDLT_SP/LDLT_DP FAIT APPEL A DEUX ROUTINES EXTERNES
        call PCShellSetSetUp(pc, ldsp1, ierr)
        ASSERT(ierr == 0)
!     Si LMP, on définit un préconditionneur à gauche et à droite
        if (lmp_is_active) then
            ASSERT(ierr == 0)
            call PCShellSetName(pc, "Symmetric Preconditionner: Left LDLT_XP, right LMP", ierr)
            ASSERT(ierr == 0)
            call PCShellSetApplySymmetricLeft(pc, ldsp2, ierr)
            ASSERT(ierr == 0)
            call PCShellSetApplySymmetricRight(pc, lmp_apply_right, ierr)
            ASSERT(ierr == 0)
            call KSPSetPCSide(ksp, PC_SYMMETRIC, ierr)
            ASSERT(ierr == 0)
!        Pour le préconditionneur LDLT_XP, on définit reac_lmp à partir de reac_precond
            reac_lmp = reacpr/2
        else
            call PCShellSetName(pc, "LDLT_XP Preconditionner", ierr)
            call PCShellSetApply(pc, ldsp2, ierr)
            ASSERT(ierr == 0)
        end if
!
        ASSERT(spmat == ' ')
        spmat = nomat
        ASSERT(spsolv == ' ')
        spsolv = nosolv
!-----------------------------------------------------------------------
    else if (precon == 'JACOBI') then
        call PCSetType(pc, PCJACOBI, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else if (precon == 'SOR') then
        call PCSetType(pc, PCSOR, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else if (precon == 'FIELDSPLIT') then
        call PCSetType(pc, PCFIELDSPLIT, ierr)
        ASSERT(ierr == 0)
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
        call mfield_setup(kptsc, nonu)
        call PCSetFromOptions(pc, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else if (precon == 'ML') then
        call PCSetType(pc, PCML, ierr)
        if (ierr .ne. 0) then
            call utmess('F', 'PETSC_19', sk=precon)
        end if
        ASSERT(ierr == 0)
!     CHOIX DE LA RESTRICTION (UNCOUPLED UNIQUEMENT ACTUELLEMENT)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_ml_CoarsenScheme', 'Uncoupled', ierr)
        ASSERT(ierr == 0)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_ml_PrintLevel', '0', ierr)
        ASSERT(ierr == 0)
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
!     APPEL OBLIGATOIRE POUR PRENDRE EN COMPTE LES AJOUTS CI-DESSUS
        call PCSetFromOptions(pc, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else if (precon == 'BOOMER') then
        call PCSetType(pc, PCHYPRE, ierr)
        if (ierr .ne. 0) then
            call utmess('F', 'PETSC_19', sk=precon)
        end if
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_type', 'boomeramg', ierr)
        ASSERT(ierr == 0)
!     CHOIX DE LA RESTRICTION (PMIS UNIQUEMENT ACTUELLEMENT)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
        &   '-pc_hypre_boomeramg_coarsen_type', 'HMIS', ierr)
        ASSERT(ierr == 0)
!     CHOIX DU LISSAGE (SOR UNIQUEMENT POUR LE MOMENT)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
        &   '-pc_hypre_boomeramg_relax_type_all', 'SOR/Jacobi', ierr)
        ASSERT(ierr == 0)
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
        & '-pc_hypre_boomeramg_print_statistics', '0', ierr)
        ASSERT(ierr == 0)
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
!     APPEL OBLIGATOIRE POUR PRENDRE EN COMPTE LES AJOUTS CI-DESSUS
        call PCSetFromOptions(pc, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else if (precon == 'GAMG') then
        call PCSetType(pc, PCGAMG, ierr)
        if (ierr .ne. 0) then
            call utmess('F', 'PETSC_19', sk=precon)
        end if
!  CHOIX DE LA VARIANTE AGGREGATED
!   call PCGAMGSetType(pc, "agg", ierr)
!   ASSERT(ierr == 0)
!  CHOIX DU NOMBRE DE LISSAGES
        nsmooth = 1
        call PCGAMGSetNSmooths(pc, nsmooth, ierr)
        ASSERT(ierr == 0)
!
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_gamg_verbose', '2', ierr)
        ASSERT(ierr == 0)
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
!     APPEL OBLIGATOIRE POUR PRENDRE EN COMPTE LES AJOUTS CI-DESSUS
        call PCSetFromOptions(pc, ierr)
        ASSERT(ierr == 0)

!-----------------------------------------------------------------------
    else if (precon == 'UTILISATEUR') then
!     It should be defined by the user
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
        call PCSetFromOptions(pc, ierr)

!-----------------------------------------------------------------------
    else if (precon == 'HPDDM') then
#ifdef ASTER_PETSC_HAVE_HPDDM
!     HPDDM works only in hpc mode
        if (.not. l_parallel_matrix) call utmess('F', 'PETSC_23')
!     Set HPPDM pc
        call PCSetType(pc, PCHPDDM, ierr)
        ASSERT(ierr == 0)
!     Create the auxiliary matrix usefull for the coarse problem
        call MatCreate(PETSC_COMM_SELF, auxMat, ierr)
        ASSERT(ierr == 0)
!     Set the block size here because we have it easily
        bs = tblocs(kptsc)
        call MatSetBlockSize(auxMat, bs, ierr)
        ASSERT(ierr == 0)
!     Transfer the local code_aster matrix to the PETSc one
        call apmams(nomat, auxMat)
        ASSERT(ierr == 0)
!     Get the local dof in the global numbering
        call jeveuo(nonu//'.NUME.NULG', 'L', vi=nulg)

#if ASTER_PETSC_INT_SIZE == 4
        AS_ALLOCATE(size=int(size(nulg), 8), vi4=nulg_ip)
#else
        AS_ALLOCATE(size=int(size(nulg), 8), vi=nulg_ip)
#endif
        nulg_ip(1:size(nulg)) = (/(to_petsc_int(nulg(i)), i=1, size(nulg))/)
        call ISCreateGeneral(PETSC_COMM_SELF, to_petsc_int(size(nulg)), nulg_ip, &
                             PETSC_COPY_VALUES, auxIS, ierr)
        ASSERT(ierr == 0)

        ! call PCHPDDMDumpAuxiliaryMat(pc, auxIS, auxMat)
!     Set the Neumann matrix
        call PCHPDDMSetAuxiliaryMat(pc, auxIS, auxMat, PETSC_NULL_FUNCTION, &
                                    PETSC_NULL_INTEGER, ierr)
        ASSERT(ierr == 0)
        call ISDestroy(auxIS, ierr)
        call MatDestroy(auxMat, ierr)
        ASSERT(ierr == 0)
!     Set the PC definition
        call dismoi('TYPE_MATRICE', nomat, 'MATR_ASSE', repk=syme)
        factor = 'lu'
        if (syme .eq. 'SYMETRI') factor = 'cholesky'

        call codent(int(nbproc/2, 8), 'D', nbproc_str, 'F')
        myopt = '-prefix_push pc_hpddm_ '// &
                '-prefix_push levels_1_ '// &
                '-pc_type asm '// &
                '-sub_mat_mumps_icntl_14 5000 '// &
                '-sub_mat_mumps_icntl_24 1 '// &
                '-sub_mat_mumps_icntl_25 0 '// &
                '-sub_mat_mumps_cntl_3 1.e-50 '// &
                '-sub_mat_mumps_cntl_5 0. '// &
                '-eps_nev 30 '// &
                '-sub_pc_type '//factor//' '// &
                '-sub_pc_factor_mat_solver_type mumps '// &
                '-st_pc_factor_mat_solver_type mumps '// &
                '-st_share_sub_ksp '// &
                '-prefix_pop '// &
                '-prefix_push coarse_ '// &
                '-pc_factor_mat_solver_type mumps '// &
                '-sub_pc_type '//factor//' '// &
                '-mat_mumps_icntl_14 5000 '// &
                '-mat_mumps_icntl_24 1 '// &
                '-mat_mumps_icntl_25 0 '// &
                '-mat_mumps_cntl_3 1.e-50 '// &
                '-mat_mumps_cntl_5 0. '// &
                '-p '//trim(nbproc_str)//' '// &
                '-prefix_pop '// &
                '-define_subdomains '// &
                '-has_neumann '// &
                '-prefix_pop '
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, myopt, ierr)
        ASSERT(ierr == 0)
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, options(kptsc), ierr)
        ASSERT(ierr .eq. 0)
        call PCSetFromOptions(pc, ierr)
        ASSERT(ierr == 0)
#if ASTER_PETSC_INT_SIZE == 4
        AS_DEALLOCATE(vi4=nulg_ip)
#else
        AS_DEALLOCATE(vi=nulg_ip)
#endif
#else
        call utmess('F', 'GENERIC_1', sk='unsupported feature: HPDDM')
#endif

!-----------------------------------------------------------------------
    else if (precon == 'SANS') then
        call PCSetType(pc, PCNONE, ierr)
        ASSERT(ierr == 0)
!-----------------------------------------------------------------------
    else
        ASSERT(.false.)
    end if
!-----------------------------------------------------------------------
!
!  CREATION EFFECTIVE DU PRECONDITIONNEUR
    call PCSetUp(pc, ierr)
!  ANALYSE DU CODE RETOUR
    if (ierr .ne. 0) then
        if ((precon == 'LDLT_SP') .or. (precon == 'LDLT_DP')) then
!        ERREUR : PCENT_PIVOT PAS SUFFISANT
            call utmess('F', 'PETSC_15')
        else
            call utmess('F', 'PETSC_14')
        end if
    end if
!
999 continue
!
    call jedema()
!
#else
    integer(kind=8) :: idummy
    idummy = kptsc
    idummy = reac_lmp
#endif
!
end subroutine
