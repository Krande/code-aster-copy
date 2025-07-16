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

subroutine apmams(matass, auxMat)
    !
#include "asterf_types.h"
#include "asterf_petsc.h"
    !
    !
    use aster_petsc_module
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtmchc.h"
#include "asterfort/wkvect.h"
    !----------------------------------------------------------------
    !
    !  CONSTRUCTION DE LA MATRICE DE RIGIDITE LOCALE DU SOUS-DOMAINE
    !
    !  En entrée : la matrice ASTER complète
    !  En sortie : La matrice PETSc définie dans l'appelant
    !
    !  Rq :
    !  - la matrice construite comporte tous les degrés de liberté du sous-domaine
    !    y compris les ghosts
    !  - dans la terminologie de la decomposition de domaine, la matrice construite
    !    est la matrice de Neumann
    !----------------------------------------------------------------
    !
#ifdef ASTER_HAVE_PETSC
    character(len=19) :: matass
    Mat :: auxMat
    !
    !     VARIABLES LOCALES
    integer(kind=8) :: nsmdi, nsmhc, nz, nvalm, nlong
    integer(kind=8) :: jdval1, jdval2, jvalm, jvalm2
    integer(kind=8) :: k, ilig, nzdeb, nzfin, bs
    integer(kind=8) :: iterm, jterm, neq2, nbd
    integer(kind=8) :: jcol1, jcol2
    !
    character(len=16), parameter :: idxi1 = '&&APMAMS.IDXI1__', idxi2 = '&&APMAMS.IDXI2__'
    character(len=16), parameter :: trans1 = '&&APMAMS.TRANS1_', trans2 = '&&APMAMS.TRANS2_'
    character(len=14) :: nume
    !
    aster_logical :: lmnsy
    !
    real(kind=8) :: valm
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()
    character(len=24), pointer :: refa(:) => null()
    !
    PetscInt, pointer :: v_dxi1(:) => null()
    PetscInt, pointer :: v_dxi2(:) => null()
    !
    !----------------------------------------------------------------
    !     Variables PETSc
    PetscInt :: unused_nz
    PetscInt :: mm, nn, bs4
    PetscErrorCode ::  ierr
    PetscInt, parameter :: ione = 1
    !----------------------------------------------------------------
    call jemarq()
    !
    call dismoi('NOM_NUME_DDL', matass, 'MATR_ASSE', repk=nume)
!
    call jeveuo(nume//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nume//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nume//'.SMOS.SMHC', 'L', vi4=smhc)
    call jelira(nume//'.SMOS.SMHC', 'LONMAX', nsmhc)
!   neq2: nb total de ddls locaux au sous-domaine
    neq2 = nsmdi
!   nz: nombre de termes non-nuls (dans la partie triangulaire superieure
!                                             ou inferieure de la matrice)
    nz = smdi(neq2)
!
!   la matrice est-elle symetrique ?
!   ---------------------------------
    call jelira(matass//'.VALM', 'NMAXOC', nvalm)
    if (nvalm .eq. 1) then
        lmnsy = .false.
    else if (nvalm .eq. 2) then
        lmnsy = .true.
    else
        ASSERT(.false.)
    end if

!    elimination des ddls (affe_char_cine)
!   ---------------------------------
    call jeveuo(matass//'.REFA', 'E', vk24=refa)
    if (refa(3) .eq. 'ELIML') call mtmchc(matass, 'ELIMF')
    ASSERT(refa(3) .ne. 'ELIML')

!   les valeurs de la partie triangulaire superieure de la matrice sont stockees
!   dans valm
!   -----------------------------------------------------------------------------
    call jeveuo(jexnum(matass//'.VALM', 1_8), 'L', jvalm)
    call jelira(jexnum(matass//'.VALM', 1_8), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
!   si la matrice n'est pas symetrique, on a aussi besoin des valeurs de
!   la partie triangulaire inferieure
    if (lmnsy) then
        call jeveuo(jexnum(matass//'.VALM', 2_8), 'L', jvalm2)
        call jelira(jexnum(matass//'.VALM', 2_8), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if

!   Usefull vectors for the transfer
!   *WARNING* the indices begin at 0 (C PETSc's convention)
!   -----------------------------------------------------------------------------
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect(idxi1, 'V V S', neq2, vi4=v_dxi1)
    call wkvect(idxi2, 'V V S', neq2, vi4=v_dxi2)
#else
    call wkvect(idxi1, 'V V I', neq2, vi=v_dxi1)
    call wkvect(idxi2, 'V V I', neq2, vi=v_dxi2)
#endif
    call wkvect(trans1, 'V V R', neq2, jdval1)
    call wkvect(trans2, 'V V R', neq2, jdval2)

!   Define the essential attributes of the auxilliary matrix
!   -----------------------------------------------------------------------------
    call MatSetType(auxMat, MATSEQAIJ, ierr)
    ASSERT(ierr == 0)
    call MatGetBlockSize(auxMat, bs4, ierr)
    ASSERT(ierr == 0)
    bs = bs4
    ASSERT(mod(neq2, bs) .eq. 0)
    call MatSetSizes(auxMat, PETSC_DECIDE, PETSC_DECIDE, to_petsc_int(neq2), &
                     to_petsc_int(neq2), ierr)
    ASSERT(ierr .eq. 0)
    call MatSetUp(auxMat, ierr)
    ASSERT(ierr .eq. 0)

!   Compute nnz
!   -----------------------------------------------------------------------------
    do jcol2 = 0, neq2-1
        jcol1 = jcol2+1
        nbd = 0
        if (jcol1 .eq. 1) then
            nzdeb = 1
        else
            nzdeb = smdi(jcol1-1)+1
        end if
        nzfin = smdi(jcol1)
        do k = nzdeb, nzfin
            ilig = smhc(k)
            nbd = nbd+1
            v_dxi1(ilig) = v_dxi1(ilig)+to_petsc_int(1)
        end do
        v_dxi1(jcol2+1) = v_dxi1(jcol2+1)+to_petsc_int(nbd-1)
    end do

!   Preallocate the matrix
!   -----------------------------------------------------------------------------
    unused_nz = -1
    mm = to_petsc_int(neq2)
    call MatSeqAIJSetPreallocation(auxMat, unused_nz, v_dxi1(1:mm), ierr)
    ASSERT(ierr .eq. 0)

!   Fill the matrix
!   -----------------------------------------------------------------------------
    do jcol2 = 0, neq2-1
        iterm = 0
        jterm = 0
        jcol1 = jcol2+1
        if (jcol1 .eq. 1) then
            nzdeb = 1
        else
            nzdeb = smdi(jcol1-1)+1
        end if
        nzfin = smdi(jcol1)
        do k = nzdeb, nzfin
            !       -- ilig : indice ligne (fortran) du terme courant dans la matrice Aster
            ilig = smhc(k)
            jterm = jterm+1
            !           -- si A n'est pas symetrique, on lit valm2
            if (lmnsy) then
                valm = zr(jvalm2-1+k)
            else
                !           -- si A est symetrique, on lit valm1
                valm = zr(jvalm-1+k)
            end if
            zr(jdval2+jterm-1) = valm
            v_dxi2(jterm) = to_petsc_int(ilig-1)

            iterm = iterm+1
            valm = zr(jvalm-1+k)
            zr(jdval1+iterm-1) = valm
            v_dxi1(iterm) = to_petsc_int(ilig-1)
        end do

!       we remove the extra term computed in the algorithm
        jterm = jterm-1
        mm = to_petsc_int(iterm)
        call MatSetValues(auxMat, mm, v_dxi1(1:mm), ione, [to_petsc_int(jcol2)], &
                          zr(jdval1-1+1:jdval1-1+mm), INSERT_VALUES, ierr)
        ASSERT(ierr .eq. 0)
!
        nn = to_petsc_int(jterm)
        call MatSetValues(auxMat, ione, [to_petsc_int(jcol2)], nn, v_dxi2(1:nn), &
                          zr(jdval2-1+1:jdval2-1+nn), INSERT_VALUES, ierr)
        ASSERT(ierr .eq. 0)
!
    end do
!
!   Assembly the auxiliary matrix
!   -----------------------------------------------------------------------------
    call MatAssemblyBegin(auxMat, MAT_FINAL_ASSEMBLY, ierr)
    ASSERT(ierr == 0)
    call MatAssemblyEnd(auxMat, MAT_FINAL_ASSEMBLY, ierr)
    ASSERT(ierr == 0)

    call jelibe(nume//'.SMOS.SMDI')
    call jelibe(nume//'.SMOS.SMHC')
    call jelibe(jexnum(matass//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(matass//'.VALM', 2_8))

    !   -- menage :
    call jedetr(idxi1)
    call jedetr(idxi2)
    call jedetr(trans1)
    call jedetr(trans2)

    call jedema()

#else
    integer(kind=8) :: idummy
    character(len=19) :: matass
    integer(kind=8) :: auxMat
#endif
    !
end subroutine
