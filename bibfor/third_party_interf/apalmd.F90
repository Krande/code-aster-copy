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

subroutine apalmd(kptsc)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
!
! person_in_charge: nicolas.sellenet at edf.fr
    use aster_petsc_module
    use petsc_data_module

    implicit none

#include "jeveux.h"
#include "asterf_debug.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/apbloc.h"
#include "asterfort/asmpi_comm_point.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  CREATION DE LA MATRICE PETSC (INSTANCE NUMERO KPTSC)
!  PREALLOCATION DANS LE CAS MATR_DISTRIBUEE
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: rang, nbproc, jnbjoi, nbjoin, jnequ, jnequl, jnugll
    integer(kind=8) :: nsmdi, nsmhc, nz, jprddl, nloc, nglo, jcoll
    integer(kind=8) :: jsmdi, jsmhc, ndprop, procol
    integer(kind=8) :: k, nzdeb, nzfin, jcolg, iligl, jaux
    integer(kind=8) :: prolig, iligg, iaux, numpro, jjoint, numloc
    integer(kind=8) :: lgenvo, numglo, comple
    mpi_int :: mpicou
!
    PetscInt, pointer :: v_idxd(:) => null()
    PetscInt, pointer :: v_idxo(:) => null()
    PetscInt, pointer :: v_idxdc(:) => null()
    PetscInt, pointer :: v_idxoc(:) => null()
    PetscInt, pointer :: v_valeur(:) => null()
    PetscInt, parameter :: one = to_petsc_int(1)
!
!
    character(len=4) :: chnbjo
    character(len=14) :: nonu
    character(len=16) :: idxo, idxd, idxoc, idxdc, cpysol
    character(len=19) :: nomat, nosolv
    character(len=24) :: nojoin
!
    parameter(idxo='&&APALMD.IDXO___')
    parameter(idxd='&&APALMD.IDXD___')
    parameter(idxoc='&&APALMD.IDXOC__')
    parameter(idxdc='&&APALMD.IDXDC__')
    parameter(cpysol='&&APALMD.COPYSOL')
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscInt :: low, high, unused_nz
    PetscErrorCode ::  ierr
    integer(kind=8) :: neql, neqg, bs
    Vec :: tmp
    mpi_int :: mrank, msize
!----------------------------------------------------------------
    call jemarq()
!
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicou)
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    DEBUG_MPI('apalmd', rang, nbproc)
!
!     -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
!
    call jeveuo(nonu//'.NUML.JOIN', 'L', jnbjoi)
    call jelira(nonu//'.NUML.JOIN', 'LONMAX', nbjoin)
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', jsmdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    nz = zi(jsmdi-1+nsmdi)
!
    call apbloc(kptsc)
    bs = tblocs(kptsc)
    ASSERT(bs .ge. 1)
!
    call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
    call jeveuo(nonu//'.NUML.NEQU', 'L', jnequl)
    call jeveuo(nonu//'.NUML.NLGP', 'L', jnugll)
    call jeveuo(nonu//'.NUML.PDDL', 'L', jprddl)
!
    nloc = zi(jnequl)
    nglo = zi(jnequ)
    neqg = nglo
    neql = nloc
!
    ndprop = 0
!
    do jcoll = 1, nloc
        procol = zi(jprddl-1+jcoll)
        if (procol .eq. rang) ndprop = ndprop+one
    end do
    call VecCreateMPI(mpicou, to_petsc_int(ndprop), to_petsc_int(neqg), tmp, ierr)
    ASSERT(ierr .eq. 0)
!
    call VecGetOwnershipRange(tmp, low, high, ierr)
    ASSERT(ierr .eq. 0)
    call VecDestroy(tmp, ierr)
    ASSERT(ierr .eq. 0)
!
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect(idxd, 'V V S', ndprop, vi4=v_idxd)
    call wkvect(idxo, 'V V S', ndprop, vi4=v_idxo)
    call wkvect(idxdc, 'V V S', nloc, vi4=v_idxdc)
    call wkvect(idxoc, 'V V S', nloc, vi4=v_idxoc)
#else
    call wkvect(idxd, 'V V I', ndprop, vi=v_idxd)
    call wkvect(idxo, 'V V I', ndprop, vi=v_idxo)
    call wkvect(idxdc, 'V V I', nloc, vi=v_idxdc)
    call wkvect(idxoc, 'V V I', nloc, vi=v_idxoc)
#endif
!
    jcolg = zi(jnugll)
    if (zi(jprddl) .eq. rang) then
        v_idxd(jcolg-low) = v_idxd(jcolg-low)+one
    else
        v_idxdc(1) = v_idxdc(1)+one
    end if
!
!     ON COMMENCE PAR NOTER DDL PAR DDL LE NOMBRE DE TERMES POSSEDES
!     ET CEUX QU'IL FAUDRA ENVOYER AUX AUTRES PROCESSEURS
    do jcoll = 2, nloc
        nzdeb = zi(jsmdi+jcoll-2)+1
        nzfin = zi(jsmdi+jcoll-1)
        procol = zi(jprddl+jcoll-1)
        jcolg = zi(jnugll+jcoll-1)
        do k = nzdeb, nzfin
            iligl = zi4(jsmhc-1+k)
            prolig = zi(jprddl-1+iligl)
            iligg = zi(jnugll-1+iligl)
!         SOIT LA COLONNE ET LA LIGNE APPARTIENNENT AU PROC COURANT
!         AUQUEL CAS, ON S'EN PREOCCUPE POUR L'ALLOCATION
            if (procol .eq. rang .and. prolig .eq. rang) then
                v_idxd(iligg-low) = v_idxd(iligg-low)+one
                if (iligg .ne. jcolg) then
                    v_idxd(jcolg-low) = v_idxd(jcolg-low)+one
                end if
!           SOIT ILS N'APPARTIENNENT PAS AU PROC COURANT TOUS LES
!         DEUX, DANS CE CAS ON LES OUBLIE
            else if (procol .ne. rang .and. prolig .ne. rang) then
                if (procol .eq. prolig) then
                    v_idxdc(iligl) = v_idxdc(iligl)+one
                    if (iligg .ne. jcolg) then
                        v_idxdc(jcoll) = v_idxdc(jcoll)+one
                    end if
                else
                    v_idxoc(iligl) = v_idxoc(iligl)+one
                    if (iligg .ne. jcolg) then
                        v_idxoc(jcoll) = v_idxoc(jcoll)+one
                    end if
                end if
!         SOIT L'UN DES DEUX APPARTIENT AU PROC COURANT
!         DANS CE CAS, ON LE COMPTE POUR L'ALLOCATION
!         OU ON PREVIENT L'AUTRE PROC
            else
                if (procol .eq. rang) then
                    v_idxo(jcolg-low) = v_idxo(jcolg-low)+one
                    v_idxoc(iligl) = v_idxoc(iligl)+one
                else
                    v_idxo(iligg-low) = v_idxo(iligg-low)+one
                    v_idxoc(jcoll) = v_idxoc(jcoll)+one
                end if
            end if
        end do
    end do
!
    do iaux = 1, nbjoin
        numpro = zi(jnbjoi+iaux-1)
        if (numpro .ne. -1) then
            call codent(iaux, 'G', chnbjo)
!
            nojoin = nonu//'.NUML.'//chnbjo
            call jeveuo(nojoin, 'L', jjoint)
            call jelira(nojoin, 'LONMAX', lgenvo)
            if (lgenvo .gt. 0) then
#if ASTER_PETSC_INT_SIZE == 4
                call wkvect(cpysol, 'V V S', lgenvo, vi4=v_valeur)
#else
                call wkvect(cpysol, 'V V I', lgenvo, vi=v_valeur)
#endif
!
                if (numpro .gt. rang) then
#if ASTER_PETSC_INT_SIZE == 4
                    call asmpi_comm_point('MPI_RECV', 'I4', numpro, iaux, nbval=lgenvo, &
                                          vi4=v_valeur)
#else
                    call asmpi_comm_point('MPI_RECV', 'I', numpro, iaux, nbval=lgenvo, &
                                          vi=v_valeur)
#endif
                    do jaux = 1, lgenvo
                        numloc = zi(jjoint+jaux-1)
                        numglo = zi(jnugll+numloc-1)
                        v_idxo(numglo-low) = v_idxo(numglo-low)+v_valeur(jaux)
                    end do
!
#if ASTER_PETSC_INT_SIZE == 4
                    call asmpi_comm_point('MPI_RECV', 'I4', numpro, iaux, nbval=lgenvo, &
                                          vi4=v_valeur)
#else
                    call asmpi_comm_point('MPI_RECV', 'I', numpro, iaux, nbval=lgenvo, &
                                          vi=v_valeur)
#endif
                    do jaux = 1, lgenvo
                        numloc = zi(jjoint+jaux-1)
                        numglo = zi(jnugll+numloc-1)
                        v_idxd(numglo-low) = v_idxd(numglo-low)+v_valeur(jaux)
                    end do
                else if (numpro .lt. rang) then
                    do jaux = 1, lgenvo
                        numloc = zi(jjoint+jaux-1)
                        v_valeur(jaux) = v_idxoc(numloc)
                    end do
#if ASTER_PETSC_INT_SIZE == 4
                    call asmpi_comm_point('MPI_SEND', 'I4', numpro, iaux, nbval=lgenvo, &
                                          vi4=v_valeur)
#else
                    call asmpi_comm_point('MPI_SEND', 'I', numpro, iaux, nbval=lgenvo, &
                                          vi=v_valeur)
#endif
!
                    do jaux = 1, lgenvo
                        numloc = zi(jjoint+jaux-1)
                        v_valeur(jaux) = v_idxdc(numloc)
                    end do
#if ASTER_PETSC_INT_SIZE == 4
                    call asmpi_comm_point('MPI_SEND', 'I4', numpro, iaux, nbval=lgenvo, &
                                          vi4=v_valeur)
#else
                    call asmpi_comm_point('MPI_SEND', 'I', numpro, iaux, nbval=lgenvo, &
                                          vi=v_valeur)
#endif
                else
                    ASSERT(.false.)
                end if
                call jedetr(cpysol)
            end if
        end if
    end do
!
    comple = nglo-ndprop
    do iaux = 1, ndprop
        v_idxd(iaux) = min(v_idxd(iaux), to_petsc_int(ndprop))
        v_idxo(iaux) = min(v_idxo(iaux), to_petsc_int(comple))
    end do
!
    call MatCreate(mpicou, ap(kptsc), ierr)
    ASSERT(ierr .eq. 0)
    call MatSetSizes(ap(kptsc), to_petsc_int(ndprop), to_petsc_int(ndprop), &
                     to_petsc_int(neqg), to_petsc_int(neqg), ierr)
    ASSERT(ierr .eq. 0)
    call MatSetType(ap(kptsc), MATMPIAIJ, ierr)
    ASSERT(ierr .eq. 0)
!
    call MatSetBlockSize(ap(kptsc), to_petsc_int(bs), ierr)
    ASSERT(ierr .eq. 0)
!
    unused_nz = -1
    call MatMPIAIJSetPreallocation(ap(kptsc), unused_nz, v_idxd(1:ndprop), &
                                   unused_nz, v_idxo(1:ndprop), ierr)
    ASSERT(ierr .eq. 0)
!
!
    call jedetr(idxd)
    call jedetr(idxo)
    call jedetr(idxdc)
    call jedetr(idxoc)
!
    call jedema()
!
#else
    integer(kind=8) :: idummy
    idummy = kptsc
#endif
!
end subroutine
