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

subroutine apmamh(kptsc)
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
#include "asterfort/crnustd.h"
#include "asterfort/jeexin.h"
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
!    matrice PETSc (ie le bloc de lignes A(low:high-1))
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: nz, nvalm, nlong
    integer(kind=8) :: jsmdi, jsmhc, jdval1, jdval2, jvalm, jvalm2
    integer(kind=8) :: k, iligl, jcoll, nzdeb, nzfin, nuno1, nucmp1, nuno2, nbproc, numno1, numno2
    integer(kind=8) :: jcolg, iligg, jnugll, nucmp2, procol, jprddl
    integer(kind=8) :: jnequ, nloc, jdeeq, prolig, rang
    integer(kind=8) :: jterm, jcolg4(1), iterm
    integer(kind=8) :: iret, step_dbg
    mpi_int :: mrank, msize
!
    character(len=19) :: nomat, nosolv
    character(len=16) :: idxi1, idxi2, trans1, trans2
    character(len=14) :: nonu
!
    aster_logical :: lmnsy, lgive
    aster_logical:: ldebug
    integer(kind=8), save :: nstep = 0
!
    real(kind=8) :: valm, valm2
!
    PetscInt, pointer :: v_dxi1(:) => null()
    PetscInt, pointer :: v_dxi2(:) => null()
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()
!
    parameter(idxi1='&&APMAMC.IDXI1__')
    parameter(idxi2='&&APMAMC.IDXI2__')
    parameter(trans1='&&APMAMC.TRANS1_')
    parameter(trans2='&&APMAMC.TRANS2_')
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscErrorCode ::  ierr
    PetscInt :: low, high, mm, nn
    PetscInt, parameter :: un = 1
    Mat :: a
!----------------------------------------------------------------
    call jemarq()
!
!   -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', jsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
    call jeveuo(nonu//'.NUME.NULG', 'L', jnugll)
    call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
    call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    nloc = zi(jnequ)
    nz = zi(jsmdi-1+nloc)
!
!   Adresses needed to get the stiffness matrix wrt nodes and dof numbers (see below)
    step_dbg = -1
    ldebug = ASTER_FALSE .and. nstep == step_dbg
    nstep = nstep+1
    if (ldebug) then
        print *, "DEBUG IN APMAMH"
        call jeexin(nonu//'.NUME.NULS', iret)
        if (iret == 0) then
            call crnustd(nonu)
        end if
        call jeveuo(nonu//'.NUME.NULS', 'L', vi=v_nuls)
        call jeveuo(nonu//'.NUME.DEEG', 'L', vi=v_deeg)
    end if
!
    call jelira(nomat//'.VALM', 'NMAXOC', nvalm)
    if (nvalm .eq. 1) then
        lmnsy = ASTER_FALSE
    else if (nvalm .eq. 2) then
        lmnsy = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(jexnum(nomat//'.VALM', 1_8), 'L', jvalm)
    call jelira(jexnum(nomat//'.VALM', 1_8), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
    if (lmnsy) then
        call jeveuo(jexnum(nomat//'.VALM', 2_8), 'L', jvalm2)
        call jelira(jexnum(nomat//'.VALM', 2_8), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
!
!     low donne la premiere ligne stockee localement
!     high donne la premiere ligne stockee par le processus (rang+1)
!     ATTENTION ces indices commencent a zero (convention C de PETSc)
    call MatGetOwnershipRange(a, low, high, ierr)
    ASSERT(ierr .eq. 0)
!
#if ASTER_PETSC_INT_SIZE == 4
    call wkvect(idxi1, 'V V S', nloc, vi4=v_dxi1)
    call wkvect(idxi2, 'V V S', nloc, vi4=v_dxi2)
#else
    call wkvect(idxi1, 'V V I', nloc, vi=v_dxi1)
    call wkvect(idxi2, 'V V I', nloc, vi=v_dxi2)
#endif
    call wkvect(trans1, 'V V R', nloc, jdval1)
    call wkvect(trans2, 'V V R', nloc, jdval2)
!
    iterm = 0
    jterm = 0
!
!   On commence par s'occuper du nombres de NZ par ligne
!   dans le bloc diagonal
!    call MatSetOption(a, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr)
    do jcoll = 1, nloc
        if (jcoll == 1) then
            nzdeb = 1
        else
            nzdeb = zi(jsmdi+jcoll-2)+1
        end if
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
            valm = zr(jvalm+k-1)
            valm2 = valm
            if (lmnsy) valm2 = zr(jvalm2+k-1)
            nuno1 = 0
            if (zi(jdeeq+(iligl-1)*2) .gt. 0) then
                nuno1 = 1
            end if
            if (nuno1 .ne. 0 .and. nuno2 .ne. 0) then
                if (prolig .eq. rang) then
                    jterm = jterm+1
                    zr(jdval2+jterm-1) = valm
                    v_dxi2(jterm) = to_petsc_int(iligg)
!                   Writings to get the stiffness matrix wrt nodes and dof numbers
                    if (ldebug .and. valm .ne. 0.d0) then
                        numno1 = v_deeg(2*(iligl-1)+1)
                        nucmp1 = v_deeg(2*(iligl-1)+2)
                        numno2 = v_deeg(2*(jcoll-1)+1)
                        nucmp2 = v_deeg(2*(jcoll-1)+2)
                        write (11+rang, *) numno1, nucmp1, numno2, nucmp2, valm, &
                            v_nuls(iligl), v_nuls(jcoll), prolig, procol, rang
                    end if
                    if (procol .eq. rang) then
                        if (iligg .ne. jcolg) then
                            iterm = iterm+1
                            zr(jdval1+iterm-1) = valm2
                            v_dxi1(iterm) = to_petsc_int(iligg)
!                           Writings to get the stiffness matrix wrt nodes and dof numbers
                            if (ldebug .and. valm2 .ne. 0.d0) then
                                numno1 = v_deeg(2*(iligl-1)+1)
                                nucmp1 = v_deeg(2*(iligl-1)+2)
                                numno2 = v_deeg(2*(jcoll-1)+1)
                                nucmp2 = v_deeg(2*(jcoll-1)+2)
                                write (11+rang, *) numno2, nucmp2, numno1, nucmp1, valm2, &
                                    v_nuls(jcoll), v_nuls(iligl), prolig, procol, rang
                            end if
                        end if
                    end if
                else if (procol .eq. rang) then
                    iterm = iterm+1
                    zr(jdval1+iterm-1) = valm2
                    v_dxi1(iterm) = to_petsc_int(iligg)
!                   Writings to get the stiffness matrix wrt nodes and dof numbers
                    if (ldebug .and. valm2 .ne. 0.d0) then
                        numno1 = v_deeg(2*(iligl-1)+1)
                        nucmp1 = v_deeg(2*(iligl-1)+2)
                        numno2 = v_deeg(2*(jcoll-1)+1)
                        nucmp2 = v_deeg(2*(jcoll-1)+2)
                        write (11+rang, *) numno2, nucmp2, numno1, nucmp1, valm2, &
                            v_nuls(jcoll), v_nuls(iligl), prolig, procol, rang
                    end if
                end if
            else
!               Si on est sur un ddl de Lagrange et qu'on possede le ddl d'en face
!               ou que les deux ddl sont de Lagrange, on doit donner le terme
                lgive = (nuno1 .eq. 0 .and. procol .eq. rang) .or. &
                        (nuno2 .eq. 0 .and. prolig .eq. rang) .or. &
                        (nuno1 .eq. 0 .and. nuno2 .eq. 0)
                if (lgive) then
                    jterm = jterm+1
                    zr(jdval2+jterm-1) = valm
                    v_dxi2(jterm) = to_petsc_int(iligg)
!                   Writings to get the stiffness matrix wrt nodes and dof numbers
                    if (ldebug .and. valm .ne. 0.d0) then
                        numno1 = v_deeg(2*(iligl-1)+1)
                        nucmp1 = v_deeg(2*(iligl-1)+2)
                        numno2 = v_deeg(2*(jcoll-1)+1)
                        nucmp2 = v_deeg(2*(jcoll-1)+2)
                        write (11+rang, *) numno1, nucmp1, numno2, nucmp2, valm, &
                            v_nuls(iligl), v_nuls(jcoll)
                        !, iligg, jcolg
                    end if
                    if (iligg .ne. jcolg) then
                        iterm = iterm+1
                        zr(jdval1+iterm-1) = valm2
                        v_dxi1(iterm) = to_petsc_int(iligg)
!                        Writings to get the stiffness matrix wrt nodes and dof numbers
                        if (ldebug .and. valm2 .ne. 0.d0) then
                            numno1 = v_deeg(2*(iligl-1)+1)
                            nucmp1 = v_deeg(2*(iligl-1)+2)
                            numno2 = v_deeg(2*(jcoll-1)+1)
                            nucmp2 = v_deeg(2*(jcoll-1)+2)
                            write (11+rang, *) numno2, nucmp2, numno1, nucmp1, valm2, &
                                v_nuls(jcoll), v_nuls(iligl)
                            !, jcolg, iligg
                        end if
                    end if
                end if
            end if
        end do
        jcolg4(1) = jcolg
        mm = to_petsc_int(jterm)
!       Ici zi4(jdxi2) donne le numero de ligne
!       Donc on donne ici le bloc triangulaire superieur
        call MatSetValues(a, mm, v_dxi2(1:mm), un, [to_petsc_int(jcolg4)], &
                          zr(jdval2:jdval2+mm), ADD_VALUES, ierr)
        nn = to_petsc_int(iterm)
!       on donne ici le bloc triangulaire inferieur
        call MatSetValues(a, un, [to_petsc_int(jcolg4)], nn, v_dxi1(1:nn), &
                          zr(jdval1:jdval1+nn), ADD_VALUES, ierr)
        iterm = 0
        jterm = 0
    end do
!
    call jelibe(nonu//'.SMOS.SMDI')
    call jelibe(nonu//'.SMOS.SMHC')
    call jelibe(jexnum(nomat//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(nomat//'.VALM', 2_8))
!
!     ON N'OUBLIE PAS DE DETRUIRE LES TABLEAUX
!     APRES AVOIR ALLOUE CORRECTEMENT
    call jedetr(idxi1)
    call jedetr(idxi2)
    call jedetr(trans1)
    call jedetr(trans2)
!   Close the logical unit dedicated to dump the matrix
    if (ldebug) flush (11+rang)
!
    call jedema()
!
#endif
!
end subroutine
