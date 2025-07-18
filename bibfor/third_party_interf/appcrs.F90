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

subroutine appcrs(kptsc, lmd)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
! person_in_charge: natacha.bereux at edf.fr
! aslint:disable=
    use aster_petsc_module
    use petsc_data_module
    implicit none

#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterc/r8prem.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/crsvfm.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMatrix.h"
!
    integer(kind=8) :: kptsc
    aster_logical :: lmd
!----------------------------------------------------------------
!
!  CREATION DU PRECONDITIONNEUR PETSC (INSTANCE NUMERO KPTSC)
!  PHASE DE RESOLUTION (RESOUD)
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!----------------------------------------------------------------
!
!     VARIABLES LOCALES
    integer(kind=8) :: rang, nbproc
    integer(kind=8) :: jslvk, jslvr, jslvi, jnequ, jnequl, jprddl, jcoll, nloc
    integer(kind=8) :: niremp, nsmdi, pcpiv, redmpi
    mpi_int :: mpicou
!
    character(len=24) :: precon, usersm, renum
    character(len=19) :: nomat, nosolv
    character(len=14) :: nonu
    character :: prec, rank
    real(kind=8) :: fillin, blreps
    aster_logical :: l_parallel_matrix
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscErrorCode ::  ierr
    PetscInt :: nlocal, first
    PetscInt :: fill, neq, ndprop
    PetscReal :: fillp
    Mat :: a
    KSP :: ksp, kspp, subksp(1)
    PC :: pc, pcp
    mpi_int :: mrank, msize
!----------------------------------------------------------------
    call jemarq()
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicou)
!
!     -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
    ksp = kp(kptsc)
!
    call jeveuo(nosolv//'.SLVK', 'L', jslvk)
    call jeveuo(nosolv//'.SLVR', 'L', jslvr)
    call jeveuo(nosolv//'.SLVI', 'L', jslvi)
    precon = zk24(jslvk-1+2)
    fillin = zr(jslvr-1+3)
    niremp = zi(jslvi-1+4)
!
    fill = to_petsc_int(niremp)
    fillp = fillin
!
!     -- RECUPERE LE RANG DU PROCESSUS ET LE NB DE PROCS
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)

!  Est-ce que la matrice est distribuee dans asterxx?
    l_parallel_matrix = isParallelMatrix(nomat)
!
!     -- CAS PARTICULIER (LDLT_INC/SOR)
!     -- CES PC NE SONT PAS PARALLELISES
!     -- ON UTILISE DONC DES VERSIONS PAR BLOC
!     ----------------------------------------
    if ((precon .eq. 'LDLT_INC') .or. (precon .eq. 'SOR')) then
        if (nbproc .gt. 1) then
            kspp = ksp
            call KSPGetPC(kspp, pcp, ierr)
            ASSERT(ierr .eq. 0)
            call PCSetType(pcp, PCBJACOBI, ierr)
            ASSERT(ierr .eq. 0)
            call KSPSetUp(kspp, ierr)
            ASSERT(ierr .eq. 0)
            call PCBJacobiGetSubKSP(pcp, nlocal, first, (/PETSC_NULL_KSP/), ierr)
            ASSERT(ierr .eq. 0)
            ASSERT(nlocal == 1)
            call PCBJacobiGetSubKSP(pcp, nlocal, first, subksp, ierr)
            ASSERT(ierr .eq. 0)
            ksp = subksp(1)
        else
            goto 999
        end if
    end if
!
!     -- choix du preconditionneur :
!     -------------------------------
    call KSPGetPC(ksp, pc, ierr)
    ASSERT(ierr .eq. 0)
!-----------------------------------------------------------------------
    if (precon .eq. 'LDLT_INC') then
        call PCSetType(pc, PCILU, ierr)
        ASSERT(ierr .eq. 0)
        call PCFactorSetLevels(pc, fill, ierr)
        ASSERT(ierr .eq. 0)
        call PCFactorSetFill(pc, fillp, ierr)
        ASSERT(ierr .eq. 0)
        call PCFactorSetMatOrderingType(pc, MATORDERINGNATURAL, ierr)
        ASSERT(ierr .eq. 0)
!-----------------------------------------------------------------------
    else if ((precon .eq. 'LDLT_SP') .or. (precon .eq. 'LDLT_DP')) then
!       CREATION SOLVEUR BIDON SIMPLE PRECISION/LOW_RANK
        spsomu = zk24(jslvk-1+3) (1:19)
        pcpiv = zi(jslvi-1+7)
        usersm = zk24(jslvk-1+9)
        blreps = zr(jslvr-1+4)
        renum = zk24(jslvk-1+4)
        redmpi = zi(jslvi-1+1)
        if (precon == 'LDLT_SP') then
            prec = 'S'
        else if (precon == 'LDLT_DP') then
            prec = 'D'
        end if
        if (blreps < r8prem()) then
            rank = 'F'
        else
            rank = 'L'
        end if
        call crsvfm(spsomu, nomat, prec, rank, pcpiv, usersm, blreps, renum, redmpi)
!        CREATION DES VECTEURS TEMPORAIRES UTILISES DANS LDLT_SP
        if (lmd .or. l_parallel_matrix) then
            if (lmd) then
                call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
                call jeveuo(nonu//'.NUML.NEQU', 'L', jnequl)
                call jeveuo(nonu//'.NUML.PDDL', 'L', jprddl)
                nloc = zi(jnequl)
                neq = to_petsc_int(zi(jnequ))
            else
                call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
                call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
                nloc = zi(jnequ)
                neq = to_petsc_int(zi(jnequ+1))
            end if
            ndprop = 0
            do jcoll = 0, nloc-1
                if (zi(jprddl+jcoll) .eq. rang) ndprop = ndprop+to_petsc_int(1)
            end do
!
            ASSERT(xlocal == PETSC_NULL_VEC)
            call VecCreateMPI(mpicou, ndprop, neq, xlocal, ierr)
        else
            call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
            neq = to_petsc_int(nsmdi)
            ASSERT(xlocal == PETSC_NULL_VEC)
            call VecCreateMPI(mpicou, PETSC_DECIDE, neq, xlocal, ierr)
        end if
        ASSERT(ierr .eq. 0)
!
        ASSERT(xscatt == PETSC_NULL_VECSCATTER)
        ASSERT(xglobal == PETSC_NULL_VEC)
        call VecCreateSeq(PETSC_COMM_SELF, to_petsc_int(neq), xglobal, ierr)
        ASSERT(ierr == 0)
! On passe PETSC_NULL_VEC en lieu et place de xglobal en entrée de VecScatterCreateToAll
! car ainsi, il n'est pas realloué
        call VecScatterCreateToAll(xlocal, xscatt, PETSC_NULL_VEC, ierr)
        ASSERT(ierr .eq. 0)
!-----------------------------------------------------------------------
    else if (precon .eq. 'SOR') then
        call PCSetType(pc, PCSOR, ierr)
        ASSERT(ierr .eq. 0)
    end if
!-----------------------------------------------------------------------
!
!     CREATION EFFECTIVE DES PRECONDITIONNEURS RETARDES
    if ((precon .eq. 'LDLT_INC') .or. (precon .eq. 'SOR')) then
        call PCSetUp(pc, ierr)
        if (ierr .ne. 0) then
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
    aster_logical :: ldummy
    idummy = kptsc
    ldummy = lmd
#endif
!
end subroutine
