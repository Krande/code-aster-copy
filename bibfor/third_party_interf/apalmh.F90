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

! This file is part of code_aster.
!
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine apalmh(kptsc)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
!
! person_in_charge: natacha.bereux at edf.fr
    use aster_petsc_module
    use petsc_data_module
    implicit none

#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterfort/apbloc.h"
#include "asterfort/assert.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/asmpi_info.h"
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  CREATION DE LA MATRICE PETSC (INSTANCE NUMERO KPTSC)
!  PREALLOCATION DANS LE CAS GENERAL
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: nsmdi, nsmhc, ndprop, nz, bs
    integer(kind=8) :: jsmdi, jsmhc, procol, prolig
    integer(kind=8) :: k, nzdeb, nzfin, jdelg, jmdla
    integer(kind=8) :: jnequ, iligl, jcoll, iligg, jcolg
    integer(kind=8) :: rang, nbproc, jprddl, jnugll, nloc, nglo
    integer(kind=8) :: nuno1, nuno2, num_ddl_max, imult, ipos, ibid
    integer(kind=8) :: num_ddl_min, iret, nblag, jdeeq
    integer(kind=8) :: imults
!
    mpi_int :: mrank, msize, mpicou
!
    character(len=19) :: nomat, nosolv
    character(len=16) :: idxo, idxd
    character(len=14) :: nonu
!
    parameter(idxo='&&APALMC.IDXO___')
    parameter(idxd='&&APALMC.IDXD___')
!
    PetscInt, pointer :: v_idxd(:) => null()
    PetscInt, pointer :: v_idxo(:) => null()
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscInt, parameter :: one = to_petsc_int(1)
    PetscErrorCode  :: ierr
    PetscInt :: low, high, neql, neqg, unused_nz
    Mat :: a
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

    call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
    call jeveuo(nonu//'.NUME.NULG', 'L', jnugll)
    call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
    call jeveuo(nonu//'.NUME.DELG', 'L', jdelg)
    call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
    nloc = zi(jnequ)
    nglo = zi(jnequ+1)
    neqg = to_petsc_int(nglo)
    neql = to_petsc_int(nloc)
!
#if ASTER_PETSC_INT_SIZE == 4
! maximum number of equation with short integer - use long int to remove this limit
    ASSERT(neqg <= huge(neqg))
#endif
!
!     -- RECUPERE LE RANG DU PROCESSUS ET LE NB DE PROCS
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
    ndprop = 0
!
    num_ddl_max = 0
    num_ddl_min = ismaem()
    do jcoll = 1, nloc
        procol = zi(jprddl-1+jcoll)
        if (procol .eq. rang) then
            ndprop = ndprop+1
            num_ddl_max = max(num_ddl_max, zi(jnugll+jcoll-1))
            num_ddl_min = min(num_ddl_min, zi(jnugll+jcoll-1))
        end if
    end do

    call MatCreate(mpicou, a, ierr)
    ASSERT(ierr .eq. 0)
    call MatSetSizes(a, to_petsc_int(ndprop), to_petsc_int(ndprop), &
                     to_petsc_int(neqg), to_petsc_int(neqg), &
                     ierr)
    ASSERT(ierr .eq. 0)
!
!   IL FAUT APPELER MATSETBLOCKSIZE *AVANT* MAT*SETPREALLOCATION
    call MatSetBlockSize(a, to_petsc_int(bs), ierr)
    ASSERT(ierr .eq. 0)

    call MatSetType(a, MATMPIAIJ, ierr)
    ASSERT(ierr .eq. 0)
    low = to_petsc_int(num_ddl_min)
    high = to_petsc_int(num_ddl_max+1)
!
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect(idxo, 'V V S', ndprop, vi4=v_idxo)
    call wkvect(idxd, 'V V S', ndprop, vi4=v_idxd)
#else
    call wkvect(idxo, 'V V I', ndprop, vi=v_idxo)
    call wkvect(idxd, 'V V I', ndprop, vi=v_idxd)
#endif
!
    jcolg = zi(jnugll)
    if (zi(jprddl) .eq. rang) then
        v_idxd(jcolg-low+1) = v_idxd(jcolg-low+1)+one
    end if
!
!   On commence par s'occuper du nombre de NZ par ligne
!   dans le bloc diagonal
    do jcoll = 2, nloc
        nzdeb = zi(jsmdi+jcoll-2)+1
        nzfin = zi(jsmdi+jcoll-1)
        procol = zi(jprddl+jcoll-1)
        jcolg = zi(jnugll+jcoll-1)
        nuno2 = 0
        if (zi(jdeeq+(jcoll-1)*2) .gt. 0) then
            nuno2 = 1
        end if
        do k = nzdeb, nzfin
            iligl = zi4(jsmhc+k-1)
            prolig = zi(jprddl+iligl-1)
            iligg = zi(jnugll+iligl-1)
            nuno1 = 0
            if (zi(jdeeq+(iligl-1)*2) .gt. 0) then
                nuno1 = 1
            end if
            if (procol .eq. rang .and. prolig .eq. rang) then
                v_idxd(iligg-low+1) = v_idxd(iligg-low+1)+one
                if (iligg .ne. jcolg) then
                    v_idxd(jcolg-low+1) = v_idxd(jcolg-low+1)+one
                end if
            else if (procol .ne. rang .and. prolig .eq. rang) then
                v_idxo(iligg-low+1) = v_idxo(iligg-low+1)+one
            else if (procol .eq. rang .and. prolig .ne. rang) then
                v_idxo(jcolg-low+1) = v_idxo(jcolg-low+1)+one
            end if
        end do
    end do
!
    call jeexin(nonu//'.NUME.MDLA', iret)
    jmdla = 0
    nblag = 0
    if (iret .ne. 0) then
        call jeveuo(nonu//'.NUME.MDLA', 'L', jmdla)
        call jelira(nonu//'.NUME.MDLA', 'LONMAX', nblag)
        nblag = nblag/3
        do ipos = 0, nblag-1
            iligl = zi(jmdla+ipos*3)
            imult = zi(jmdla+ipos*3+1)
            imults = zi(jmdla+ipos*3+2)-imult
            iligg = zi(jnugll+iligl-1)
!           Le but ici est de rajouter juste le bon nombre de termes
!           On utilise le nombre de fois qu'apparaissent les noeuds de Lagrange
!           dans des mailles tardives (sur tous les procs et sur les autres procs
!           que le proc courant)
!           On suppose qu'un ddl de Lagrange sera connecte aux autres ddl de la
!           même maniere que sur le proc qui les possede
!           C'est pour cette raison qu'on utilise v_idxd(iligg - low +1)
!           divise par le nombre de fois qu'apparait un noeud de Lagrange sur le
!           proc courant
!           Dans le cas des doubles Lagrange, la division (v_idxd(iligg - low +1)/imults)
!           ne tombe pas juste à cause des 2 termes diagonaux +1 et -1 qui sont là pour fixer
!           l'égalité des 2 Lagrange
!           Dans le cas des simples Lagrange, elle tombe juste
!           Le "+1" est là dans le cas de double Lagrange qui vivent chacun sur 2 processus
!           différents (il y a alors un terme de couplage dans le bloc hors-diagonal)
!           (voir issue31132)
            !ibid = (v_idxd(iligg - low +1)/imults)*(imult) + 1
            ibid = imult*6
            v_idxo(iligg-low+1) = v_idxo(iligg-low+1)+to_petsc_int(ibid)
        end do
    end if
    unused_nz = -1
    call MatMPIAIJSetPreallocation(a, unused_nz, v_idxd, &
                                   unused_nz, v_idxo, ierr)
    ASSERT(ierr .eq. 0)

    ap(kptsc) = a

    call jedetr(idxd)
    call jedetr(idxo)

    call jedema()

#else
    integer(kind=8) :: idummy
    idummy = kptsc
#endif

end subroutine
