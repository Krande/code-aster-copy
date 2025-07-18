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

subroutine apmamd(kptsc)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use petsc_data_module

    implicit none
! person_in_charge: nicolas.sellenet at edf.fr
! aslint:disable=
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
    integer(kind=8) :: kptsc
!----------------------------------------------------------------
!
!  REMPLISSAGE DE LA MATRICE PETSC (INSTANCE NUMERO KPTSC)
!  DANS LE CAS MATR_DISTRIBUEE
!
!  En entrée : la matrice locale ASTER
!  En sortie : les valeurs de la matrice PETSc sont remplies à
!              partir des valeurs de la matrice ASTER
!
!  Rq :
!  - la matrice PETSc n'a pas de stockage symétrique: que la matrice
!    ASTER soit symétrique ou non, la matrice PETSc est stockée en entier
!    (termes non-nuls).
!  - dans le mode "matrice distribuée" chaque processeur connait une partie de
!    la matrice ASTER, la matrice ASTER "locale" Aloc.
!    Il initialise sa partie de matrice PETSc (i.e. le bloc de lignes A(low:high-1)).
!    La matrice locale ASTER et la matrice locale PETSc sont différentes.
!    Lors du MatSetValues, le processeur local envoie aux autres processeurs
!    les valeurs dont ils ont besoin et il récupère les valeurs
!    lui permettant d'initialiser son bloc. C'est PETSc qui gère les
!    communications entre processeurs (qui possède quoi etc ), et cette gestion
!    est cachée. Ici chaque proc envoie *toutes* les valeurs qu'il possède, sans
!    savoir à qui.
!
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: nsmdi, nsmhc, nz, nvalm, nlong
    integer(kind=8) :: jsmdi, jsmhc, jdval1, jdval2, jvalm, jvalm2
    integer(kind=8) :: k, iligl, jcoll, nzdeb, nzfin
    integer(kind=8) :: iterm, jterm, jcolg, iligg, jnugll
    integer(kind=8) :: jnequ, nloc, nglo, jnequl
    PetscInt :: mm, nn
!
    character(len=19) :: nomat, nosolv
    character(len=16), parameter :: idxi1 = '&&APMAMD.IDXI1__', idxi2 = '&&APMAMD.IDXI2__'
    character(len=16), parameter :: trans1 = '&&APMAMD.TRANS1_', trans2 = '&&APMAMD.TRANS2_'
    character(len=14) :: nonu
!
    aster_logical :: lmnsy
!
    real(kind=8) :: valm
!
    PetscInt, pointer :: v_dxi1(:) => null()
    PetscInt, pointer :: v_dxi2(:) => null()
!
!
!----------------------------------------------------------------
!     Variables PETSc
    PetscInt :: neql, neqg
    PetscErrorCode ::  ierr
    PetscInt, parameter :: one = 1
    Mat :: a
!----------------------------------------------------------------
    call jemarq()
!
!     -- LECTURE DU COMMUN
    nomat = nomat_courant
    nonu = nonu_courant
    nosolv = nosols(kptsc)
    a = ap(kptsc)
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', jsmdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    call jeveuo(nonu//'.NUME.NEQU', 'L', jnequ)
    call jeveuo(nonu//'.NUML.NEQU', 'L', jnequl)
    call jeveuo(nonu//'.NUML.NLGP', 'L', jnugll)
    nloc = zi(jnequl)
    nglo = zi(jnequ)
    neql = to_petsc_int(nloc)
    neqg = to_petsc_int(nglo)
    nz = zi(jsmdi-1+nloc)
!
! La matrice Aster est-elle symétrique ?
    call jelira(nomat//'.VALM', 'NMAXOC', nvalm)
    if (nvalm .eq. 1) then
        lmnsy = .false.
    else if (nvalm .eq. 2) then
        lmnsy = .true.
    else
        ASSERT(.false.)
    end if
! Vérification de la cohérence entre le(s) tableau(x) stockant les
! valeurs de la matrice nomat et sa structure creuse (telle que définie
! dans nonu)
    call jeveuo(jexnum(nomat//'.VALM', 1_8), 'L', jvalm)
    call jelira(jexnum(nomat//'.VALM', 1_8), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
    if (lmnsy) then
        call jeveuo(jexnum(nomat//'.VALM', 2_8), 'L', jvalm2)
        call jelira(jexnum(nomat//'.VALM', 2_8), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
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
!  Recopie de la matrice
!  C'est PETSc qui s'occupe de la recopie des termes vers
!  le bon processeur
!
! Envoi de Aloc(1,1)
    call MatSetValue(a, to_petsc_int(zi(jnugll)-1), to_petsc_int(zi(jnugll)-1), zr(jvalm), &
                     ADD_VALUES, ierr)
    ASSERT(ierr == 0)
!
    do jcoll = 2, nloc
        nzdeb = zi(jsmdi+jcoll-2)+1
        nzfin = zi(jsmdi+jcoll-1)
! Indice colonne global (F) de la colonne locale jcoll
        jcolg = zi(jnugll+jcoll-1)
        do k = nzdeb, nzfin
            iligl = zi4(jsmhc-1+k)
            iligg = zi(jnugll-1+iligl)
! Compteur de termes sur la colonne locale jcoll
            iterm = iterm+1
            valm = zr(jvalm-1+k)
! Stockage dans val1 de A(iligg,jcolg)
            zr(jdval1+iterm-1) = valm
! et de son indice ligne global (C)
            v_dxi1(iterm) = to_petsc_int(iligg-1)
! On passe à la *ligne* jcoll
            if (iligg .ne. jcolg) then
! Attention, il ne faut pas stocker le terme diagonal A(jcolg, jcolg)
! qui a déjà été rencontré dans la *colonne* jcoll
! Compteur de termes sur la ligne jcoll
                jterm = jterm+1
                if (.not. lmnsy) then
! si la matrice ASTER est symétrique
! la ligne jcoll est la transposée de la colonne jcoll
! on reprend la valeur lue depuis valm
                    valm = zr(jvalm-1+k)
                else
! si la matrice ASTER n'est pas symétrique
! on lit les termes de la ligne jcoll depuis valm2
                    valm = zr(jvalm2-1+k)
                end if
! on stocke dans val2
                zr(jdval2+jterm-1) = valm
! avec l'indice colonne global (C) correspondant
                v_dxi2(jterm) = to_petsc_int(iligg-1)
            end if
        end do
! Envoi de la colonne jcolg
        mm = to_petsc_int(iterm)
        call MatSetValues(a, mm, v_dxi1(1:mm), one, [to_petsc_int(jcolg-1)], &
                          zr(jdval1-1+1:jdval1-1+mm), ADD_VALUES, ierr)
        ASSERT(ierr == 0)
! Envoi de la ligne jcolg
        nn = to_petsc_int(jterm)
        call MatSetValues(a, one, [to_petsc_int(jcolg-1)], nn, v_dxi2(1:nn), &
                          zr(jdval2-1+1:jdval2-1+nn), ADD_VALUES, ierr)
        ASSERT(ierr == 0)
        iterm = 0
        jterm = 0
    end do
!
    call jelibe(nonu//'.SMOS.SMDI')
    call jelibe(nonu//'.SMOS.SMHC')
    call jelibe(jexnum(nomat//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(nomat//'.VALM', 2_8))
!
    call jedetr(idxi1)
    call jedetr(idxi2)
    call jedetr(trans1)
    call jedetr(trans2)
!
    call jedema()
!
#else
    integer(kind=8) :: idummy
    idummy = kptsc
#endif
!
end subroutine
