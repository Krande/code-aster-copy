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

subroutine apvsmbh(kptsc, rsolu)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
! person_in_charge: nicolas.pignet at edf.fr
! aslint:disable=
    use aster_petsc_module
    use petsc_data_module

    implicit none

#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/jeexin.h"
#include "asterfort/crnustd.h"

    integer(kind=8) :: kptsc
    real(kind=8) :: rsolu(*)
#ifdef ASTER_HAVE_PETSC

!----------------------------------------------------------------
!
!  CREATION ET REMPLISSAGE DU SECOND MEMBRE EN HPC
!
!----------------------------------------------------------------
!
!     VARIABLES LOCALES
    integer(kind=8) :: jnequ, jnugll, jprddl, iret
    integer(kind=8) :: nloc, nglo, ndprop
    integer(kind=8) :: bs, jcoll, jvaleu, iterm, nuno, nucmp, step
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()

    mpi_int :: mpicomm
!
    character(len=14) :: nonu

    aster_logical :: dbg
    integer(kind=8), save :: nstep = 0
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscInt :: ndprop4, nn
    PetscInt, pointer :: v_indic(:) => null()
    PetscErrorCode :: ierr
    mpi_int :: rang
!----------------------------------------------------------------
    call jemarq()
!
! -- DEBUG
    step = -1
    dbg = ASTER_FALSE .and. step == nstep
    nstep = nstep+1
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicomm)
    call asmpi_info(rank=rang)
!
!     -- LECTURE DU COMMUN
    nonu = nonu_courant
    bs = tblocs(kptsc)
!
    call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
    call jeveuo(nonu//'.NUME.NULG', 'L', jnugll)
    call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
    if (dbg) then
        call jeexin(nonu//'.NUME.NULS', iret)
        if (iret == 0) then
            call crnustd(nonu)
        end if
        call jeveuo(nonu//'.NUME.NULS', 'L', vi=v_nuls)
        call jeveuo(nonu//'.NUME.DEEG', 'L', vi=v_deeg)
        print *, "DEBUG IN APVSMBH"
    end if
    nloc = zi(jnequ)
    nglo = zi(jnequ+1)
    ndprop = 0
    do jcoll = 0, nloc-1
        if (zi(jprddl+jcoll) .eq. rang) ndprop = ndprop+1
    end do
    ndprop4 = to_petsc_int(ndprop)
    call VecCreate(mpicomm, b, ierr)
    ASSERT(ierr .eq. 0)
    call VecSetBlockSize(b, to_petsc_int(bs), ierr)
    ASSERT(ierr .eq. 0)
    call VecSetSizes(b, ndprop4, to_petsc_int(nglo), ierr)
    ASSERT(ierr .eq. 0)
    call VecSetType(b, VECMPI, ierr)
    ASSERT(ierr .eq. 0)
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect('&&APMAIN.INDICES', 'V V S', ndprop, vi4=v_indic)
#else
    call wkvect('&&APMAIN.INDICES', 'V V I', ndprop, vi=v_indic)
#endif
    call wkvect('&&APMAIN.VALEURS', 'V V R', ndprop, jvaleu)
    iterm = 0
    do jcoll = 0, nloc-1
        if (zi(jprddl+jcoll) .eq. rang) then
            v_indic(iterm+1) = to_petsc_int(zi(jnugll+jcoll))
            zr(jvaleu+iterm) = rsolu(jcoll+1)
            iterm = iterm+1
!
            if (dbg) then
                nuno = v_deeg(2*(jcoll)+1)
                nucmp = v_deeg(2*(jcoll)+2)
!                    numéro noeud global, num comp du noeud, nume eq std, rhs, nume eq glob
                write (601+rang, *) nuno, nucmp, v_nuls(jcoll+1), rsolu(jcoll+1)
                !,zi(jnugll + jcoll)
            end if
        end if
    end do
    if (dbg) flush (601+rang)
    nn = to_petsc_int(iterm)
    call VecSetValues(b, nn, v_indic(1), zr(jvaleu), INSERT_VALUES, ierr)
    call jedetr('&&APMAIN.INDICES')
    call jedetr('&&APMAIN.VALEURS')
    call VecAssemblyBegin(b, ierr)
    ASSERT(ierr .eq. 0)
    call VecAssemblyEnd(b, ierr)
    ASSERT(ierr .eq. 0)

    call jedema()
!
#else
    integer(kind=8) :: idummy
    real(kind=8) :: rdummy
    idummy = kptsc
    rdummy = rsolu(1)
#endif
!
end subroutine
