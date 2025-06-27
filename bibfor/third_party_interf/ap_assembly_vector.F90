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

subroutine ap_assembly_vector(chno)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use petsc_data_module
    use saddle_point_module, only: convert_rhs_to_saddle_point

    implicit none

#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/ap_on_off.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/crnustd.h"
!
    character(len=*), intent(in) :: chno
!----------------------------------------------------------------
!
!  Le but de cette routine est de completer un vecteur
!  parallele distribue (assemblage parallele)
!
!----------------------------------------------------------------
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: rang, nbproc, jnequ, numglo, jnulg, nuno, nucmp
    integer(kind=8) :: iloc, nloc, nglo, ndprop, nval, jvale, iret
    integer(kind=8), save :: nstep = 0, step
    integer(kind=8), pointer :: pddl(:) => null()
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()

    mpi_int :: mpicomm
!
    character(len=14) :: numddl

    real(kind=8), dimension(:), pointer :: val => null()
    character(len=16) :: typsd
    character(len=19) :: cn19, nume_equa, nommai
    aster_logical :: dbg
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscErrorCode :: ierr
    mpi_int :: mrank, msize
    Vec :: assembly
    PetscInt :: low, high
    PetscScalar :: xx(1)
    PetscOffset :: xidx
    PetscBool :: done
    PetscInt, dimension(:), pointer :: ig_petsc_c => null()

!----------------------------------------------------------------
    call jemarq()
!
    cn19 = chno
    typsd = '****'
    step = 1
    dbg = .false. .and. nstep == step
    nstep = nstep+1

    call dismoi('NUME_EQUA', cn19, 'CHAM_NO', repk=nume_equa)
    call dismoi('NOM_MAILLA', nume_equa, 'NUME_EQUA', repk=nommai)
    call gettco(nommai(1:8), typsd)
    if (typsd .ne. 'MAILLAGE_P') then
        goto 999
    end if
    call jeveuo(nume_equa//'.NULG', 'L', jnulg)
    call jeveuo(nume_equa//'.PDDL', 'L', vi=pddl)
    call jeveuo(nume_equa//'.NEQU', 'L', jnequ)
    call jeveuo(cn19//'.VALE', 'E', jvale)
    nloc = zi(jnequ)
    nglo = zi(jnequ+1)
!
    if (dbg) then
        print *, "DEBUG IN AS_ASSEMBLY_VECTOR"
        call jeexin(nume_equa//'.NULS', iret)
        if (iret == 0) then
            call crnustd(numddl)
        end if
        call jeveuo(nume_equa//'.NULS', 'L', vi=v_nuls)
        call jeveuo(nume_equa//'.DEEG', 'L', vi=v_deeg)
    end if
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicomm)
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!   Nombre de ddls m'appartenant (pour PETSc)
    ndprop = count(pddl(1:nloc) == rang)
!
    call PetscInitialized(done, ierr)
    ASSERT(ierr .eq. 0)
    if (.not. done) then
        call ap_on_off('ON', " ")
    end if
    call VecCreate(mpicomm, assembly, ierr)
    ASSERT(ierr .eq. 0)
    call VecSetSizes(assembly, to_petsc_int(ndprop), to_petsc_int(nglo), ierr)
    ASSERT(ierr .eq. 0)
    call VecSetType(assembly, VECMPI, ierr)
    ASSERT(ierr .eq. 0)
!
#if ASTER_PETSC_INT_SIZE == 4
    AS_ALLOCATE(vi4=ig_petsc_c, size=nloc)
#else
    AS_ALLOCATE(vi=ig_petsc_c, size=nloc)
#endif
    AS_ALLOCATE(vr=val, size=nloc)
    nval = 0
    do iloc = 1, nloc
        if (dbg) then
            nuno = v_deeg(2*(iloc-1)+1)
            nucmp = v_deeg(2*(iloc-1)+2)
!           numéro noeud global, num comp du noeud, nume eq std, rhs, nume eq glob
            write (701+rang, *) nuno, nucmp, v_nuls(iloc), zr(jvale+iloc-1)

        end if
        !if( pddl(iloc) .eq. rang ) then
        nval = nval+1
        ig_petsc_c(nval) = to_petsc_int(zi(jnulg+iloc-1))
        val(nval) = zr(jvale+iloc-1)
        !endif
    end do
    if (dbg) flush (701+rang)
!
    call VecSetValues(assembly, to_petsc_int(nval), ig_petsc_c, val, ADD_VALUES, ierr)
    ASSERT(ierr .eq. 0)
    call VecAssemblyBegin(assembly, ierr)
    ASSERT(ierr .eq. 0)
    call VecAssemblyEnd(assembly, ierr)
    ASSERT(ierr .eq. 0)
!
    call VecGetOwnershipRange(assembly, low, high, ierr)
    ASSERT(ierr .eq. 0)
!
    call VecGetArray(assembly, xx, xidx, ierr)
    ASSERT(ierr .eq. 0)
    do iloc = 1, nloc
        if (pddl(iloc) .eq. rang) then
            numglo = zi(jnulg-1+iloc)
            zr(jvale-1+iloc) = xx(xidx+numglo-low+1)
        end if
        if (dbg) then
            nuno = v_deeg(2*(iloc-1)+1)
            nucmp = v_deeg(2*(iloc-1)+2)
!           numéro noeud global, num comp du noeud, nume eq std, rhs, nume eq glob
            write (801+rang, *) nuno, nucmp, v_nuls(iloc), zr(jvale-1+iloc)

        end if
    end do
    if (dbg) flush (801+rang)
!
    call VecRestoreArray(assembly, xx, xidx, ierr)
    ASSERT(ierr .eq. 0)
#if ASTER_PETSC_INT_SIZE == 4
    AS_DEALLOCATE(vi4=ig_petsc_c)
#else
    AS_DEALLOCATE(vi=ig_petsc_c)
#endif
    AS_DEALLOCATE(vr=val)
    call VecDestroy(assembly, ierr)
    ASSERT(ierr .eq. 0)
!
999 continue
    call jedema()
!
#else
    character(len=19) :: cn19
!
    cn19 = chno
#endif
!
end subroutine
