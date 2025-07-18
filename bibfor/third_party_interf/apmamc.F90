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

subroutine apmamc(kptsc)
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
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "asterfort/asmpi_info.h"
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  REMPLISSAGE DE LA MATRICE PETSC (INSTANCE NUMERO KPTSC)
!
!  En entrée : la matrice ASTER complète
!  En sortie : les valeurs de la matrice PETSc sont remplies à
!              partir des valeurs de la matrice ASTER
!
!  Rq :
!  - la matrice PETSc n'a pas de stockage symétrique: que la matrice
!    ASTER soit symétrique ou non, la matrice PETSc est stockée en entier
!    (termes non-nuls).
!  - dans le mode "matrice complète" (MC) tous les processeurs connaissent
!    toute la matrice ASTER. Chaque processeur initialise sa partie de la
!    matrice PETSc (ie le bloc de lignes A(low2:high2-1))
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: nsmdi, nsmhc, nz, nvalm, nlong
    integer(kind=8) :: jdval1, jdval2, jvalm, jvalm2
    integer(kind=8) :: k, ilig, nzdeb, nzfin, bs, ieq1, ieq2
    integer(kind=8) :: iterm, jterm, neq2
    integer(kind=8) :: jrefn, jdeeq, numno1, numno2, nucmp1, nucmp2, rang
!
    character(len=19) :: nomat, nosolv
    character(len=16), parameter :: idxi1 = '&&APMAMC.IDXI1__', idxi2 = '&&APMAMC.IDXI2__'
    character(len=16), parameter :: trans1 = '&&APMAMC.TRANS1_', trans2 = '&&APMAMC.TRANS2_'
    character(len=14) :: nonu
    character(len=8)  :: noma
!
    aster_logical :: lmnsy, ldebug
!
    mpi_int :: mrank, msize
!
    real(kind=8) :: valm
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=4), pointer :: smhc(:) => null()
!
    PetscInt, pointer :: v_dxi1(:) => null()
    PetscInt, pointer :: v_dxi2(:) => null()
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscInt :: low2, high2, jcol1, jcol2
    PetscInt :: mm, nn
    PetscErrorCode ::  ierr
    PetscInt, parameter :: ione = 1
    Mat :: a
!----------------------------------------------------------------
    call jemarq()
!
!   -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
    bs = tblocs(kptsc)
    ASSERT(bs .ge. 1)

    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', vi4=smhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
!   neq2: nb total de ddls
    neq2 = nsmdi
!   nz: nombre de termes non-nuls (dans la partie triangulaire superieure
!                                             ou inferieure de la matrice)
    nz = smdi(neq2)

    ASSERT(mod(neq2, bs) .eq. 0)
!
!   Adresses needed to get the stiffness matrix wrt nodes and dof numbers (see below)
    ldebug = .false.
    if (ldebug) then
        call jeveuo(nonu//'.NUME.REFN', 'L', jrefn)
        noma = zk24(jrefn) (1:8)
        call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
    end if
!
!   la matrice est-elle symetrique ?
!   ---------------------------------
    call jelira(nomat//'.VALM', 'NMAXOC', nvalm)
    if (nvalm .eq. 1) then
        lmnsy = .false.
    else if (nvalm .eq. 2) then
        lmnsy = .true.
    else
        ASSERT(.false.)
    end if

!   les valeurs de la partie triangulaire superieure de la matrice sont stockees
!   dans valm
!   -----------------------------------------------------------------------------
    call jeveuo(jexnum(nomat//'.VALM', 1_8), 'L', jvalm)
    call jelira(jexnum(nomat//'.VALM', 1_8), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
!   si la matrice n'est pas symetrique, on a aussi besoin des valeurs de
!   la partie triangulaire inferieure
    if (lmnsy) then
        call jeveuo(jexnum(nomat//'.VALM', 2_8), 'L', jvalm2)
        call jelira(jexnum(nomat//'.VALM', 2_8), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if

!   low2  Donne la premiere ligne stockee localement
!   high2 Donne la premiere ligne stockee par le processus de (rang+1)
!   *ATTENTION* ces indices commencent a zero (convention C de PETSc)
!   -----------------------------------------------------------------------------
    call MatGetOwnershipRange(a, low2, high2, ierr)
    ASSERT(ierr .eq. 0)
!
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect(idxi1, 'V V S', neq2, vi4=v_dxi1)
    call wkvect(idxi2, 'V V S', neq2, vi4=v_dxi2)
#else
    call wkvect(idxi1, 'V V I', neq2, vi=v_dxi1)
    call wkvect(idxi2, 'V V I', neq2, vi=v_dxi2)
#endif
    call wkvect(trans1, 'V V R', neq2, jdval1)
    call wkvect(trans2, 'V V R', neq2, jdval2)

!
!   Le bloc de lignes A(low2:high2-1) est compose de trois sous-blocs:
!            (                               )
!   low2     ( x x x \ o o o | v v v v v v v )
!            ( x x x x \ o o | v v v v v v v )
!   high2-1  ( x x x x x \ o | v v v v v v v )
!            (                               )
!   - bloc C (x) = lower(A(low2:high2-1))
!   - bloc D (o) = upper(A(low2:high2-1,low2:high2-1))
!   - bloc E (v) = A(low2:high2-1,high2:neq2)
!
!   upper(A) est stockee au format CSC.
!   Si A n'est pas symetrique, on stocke egalement
!   lower(A),  au format CSR
!--------------------------------------------------------------------------------

!   -- On commence par s'occuper des blocs C et D
!      Indices C : jcol2
!      Indices F : jcol1, ilig
!------------------------------------------------

    do jcol2 = low2, high2-ione
!       -- Les termes non-nuls de A(1:jcol2,jcol2) sont stockes dans valm (nzdeb:nzfin)
!          Si A n'est pas symetrique, les termes non-nuls de A(jcol2,1:jcol2) sont stockes
!          dans valm2 (nzdeb:nzfin)
        iterm = 0
        jterm = 0
        jcol1 = jcol2+ione
        if (jcol1 .eq. 1) then
            nzdeb = 1
        else
            nzdeb = smdi(jcol1-1)+1
        end if
        nzfin = smdi(jcol1)
        do k = nzdeb, nzfin
!       -- ilig : indice ligne (fortran) du terme courant dans la matrice Aster
            ilig = smhc(k)
! ======
! Bloc C
! ======
!           -- Lecture de la ligne C(jcol2,:)
!           -- Compteur de termes dans la ligne jcol2 de C
            jterm = jterm+1
!           -- si A n'est pas symetrique, on lit valm2
            if (lmnsy) then
                valm = zr(jvalm2-1+k)
            else
!           -- si A est symetrique, on lit valm1
                valm = zr(jvalm-1+k)
            end if
            zr(jdval2+jterm-1) = valm
!           -- on stocke l'indice C de la ligne, c'est
!              l'indice de la colonne transposee
            v_dxi2(jterm) = to_petsc_int(ilig-1)
!           Writings to get the stiffness matrix wrt nodes and dof numbers
            if (ldebug) then
                numno1 = zi(jdeeq+2*(ilig-1))
                numno2 = zi(jdeeq+2*(jcol1-1))
                nucmp1 = zi(jdeeq+2*(ilig-1)+1)
                nucmp2 = zi(jdeeq+2*(jcol1-1)+1)
                ieq1 = 0
                ieq2 = 0
                if (numno1 .eq. 0) ieq1 = ilig
                if (numno2 .eq. 0) ieq2 = jcol1
                write (11+rang, *) numno2, nucmp2, numno1, nucmp1, valm, ieq2, ieq1
            end if
! ======
! bloc D
! ======
!           -- il est lu en colonne depuis valm
            if (ilig .ge. (low2+1)) then
!              -- Compteur de termes dans la colonne jcol2 de D
                iterm = iterm+1
                valm = zr(jvalm-1+k)
                zr(jdval1+iterm-1) = valm
!               -- on stocke l'indice C de la ligne
                v_dxi1(iterm) = to_petsc_int(ilig-1)
!               Writings to get the stiffness matrix wrt nodes and dof numbers
                if (ldebug) then
                    numno1 = zi(jdeeq+2*(ilig-1))
                    numno2 = zi(jdeeq+2*(jcol1-1))
                    nucmp1 = zi(jdeeq+2*(ilig-1)+1)
                    nucmp2 = zi(jdeeq+2*(jcol1-1)+1)
                    ieq1 = 0
                    ieq2 = 0
                    if (numno1 .eq. 0) ieq1 = ilig
                    if (numno2 .eq. 0) ieq2 = jcol1
                    write (11+rang, *) numno1, nucmp1, numno2, nucmp2, valm, ieq1, ieq2
                end if
            end if
        end do

!       -- On enleve un terme dans la ligne C(jcol2,:): c'est le terme diagonal
!         que l'on a stocke deux fois (pour C et pour D)
        jterm = jterm-1

!       -- Valeurs de D => on envoie les valeurs de la colonne jcol2
        mm = to_petsc_int(iterm)
        call MatSetValues(a, mm, v_dxi1(1:mm), ione, [to_petsc_int(jcol2)], &
                          zr(jdval1-1+1:jdval1-1+mm), INSERT_VALUES, ierr)
        ASSERT(ierr .eq. 0)

!       -- Valeurs de C => on envoie les valeurs de la ligne jcol2
        nn = to_petsc_int(jterm)
        call MatSetValues(a, ione, [to_petsc_int(jcol2)], nn, v_dxi2(1:nn), &
                          zr(jdval2-1+1:jdval2-1+nn), INSERT_VALUES, ierr)
        ASSERT(ierr .eq. 0)
!
    end do

!  -- Ensuite on finit par le bloc hors diagonal E
!      Indices C : jcol2
!      Indices F : jcol1, ilig
!  -------------------------------------------------
!
!   -- On lit colonne par colonne upper(A( :,high2:))
    do jcol2 = high2, to_petsc_int(neq2-1)
        iterm = 0
        jcol1 = jcol2+ione
        ASSERT(jcol1 .ge. 2)
        nzdeb = smdi(jcol1-1)+1
        nzfin = smdi(jcol1)
        do k = nzdeb, nzfin
            ilig = smhc(k)
!           -- On ignore les lignes avant low2
            if (ilig .lt. (low2+1)) then
                continue
!           -- On lit et on stocke A(low2+1:high2,jcol2)= E(:,jcol2)
            else if (ilig .le. high2) then
                iterm = iterm+1
                valm = zr(jvalm-1+k)
                zr(jdval1+iterm-1) = valm
                v_dxi1(iterm) = to_petsc_int(ilig-1)
!               Writings to get the stiffness matrix wrt nodes and dof numbers
                if (ldebug) then
                    numno1 = zi(jdeeq+2*(ilig-1))
                    numno2 = zi(jdeeq+2*(jcol1-1))
                    nucmp1 = zi(jdeeq+2*(ilig-1)+1)
                    nucmp2 = zi(jdeeq+2*(jcol1-1)+1)
                    ieq1 = 0
                    ieq2 = 0
                    if (numno1 .eq. 0) ieq1 = ilig
                    if (numno2 .eq. 0) ieq2 = jcol1
                    write (11+rang, *) numno1, nucmp1, numno2, nucmp2, valm, ieq1, ieq2
                end if
            else
!               -- On ignore les lignes après high2
                exit
            end if
        end do
!       -- Valeurs de E => on envoie les valeurs de la colonne jcol2
        mm = to_petsc_int(iterm)
        nn = to_petsc_int(jcol2)
        call MatSetValues(a, mm, v_dxi1(1:mm), ione, [nn], &
                          zr(jdval1-1+1:jdval1-1+nn), INSERT_VALUES, ierr)
    end do

    call jelibe(nonu//'.SMOS.SMDI')
    call jelibe(nonu//'.SMOS.SMHC')
    call jelibe(jexnum(nomat//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(nomat//'.VALM', 2_8))

!   -- menage :
    call jedetr(idxi1)
    call jedetr(idxi2)
    call jedetr(trans1)
    call jedetr(trans2)

! -- Close the logical unit dedicated to dump the matrix
    if (ldebug) flush (11+rang)

    call jedema()

#else
    integer(kind=8) :: idummy
    idummy = kptsc
#endif
!
end subroutine
