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

subroutine elg_calc_matm_red(matas1, matas2, bas1)
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use elg_data_module
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcncon.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/settco.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=19) :: matas1, matas2
    character(len=1) :: bas1
!--------------------------------------------------------------
! BUT :
!   Créer la matrice aster matas2 correspondant à la matrice
!   Kproj de matas1 ainsi que son nume_ddl
!
! in/jxin   : matas1 : sd_matr_asse
! in/jxout  : matas2 : sd_matr_asse "reduite" (sans lagranges)
! in        : bas1 : 'G'/'V' (pour la création de matas2)
!
! Remarque : on crée egalement un nume_ddl (sous-terrain) pour
!            matas2.
!---------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!================================================================
    character(len=1) :: kbid
    character(len=8) :: ma, nomgds
    character(len=14) :: nu1, nu2
    character(len=19) :: ligrmo
    integer(kind=8) :: iterm, ico, kdiag, neq2, nnz2, n1
    integer(kind=8) :: ilig, jcol, jnzcol, jsmdi2, ndiag, k, jsmhc2
    integer(kind=8) :: jvalm, jvalm2, j1, nbnl, nbno
    integer(kind=8) :: nbnom, jdeeq2, nec
    aster_logical :: lsym
    PetscInt :: nterm, mm, nn
    PetscErrorCode :: ierr
    PetscInt, allocatable :: irow(:)
    real(kind=8), allocatable :: vrow(:)
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=24), pointer :: refa(:) => null()
!----------------------------------------------------------------
    call jemarq()
!
!
!     1. On dupplique les 2 SD => NU2, MATAS2 :
!     -----------------------------------------
    call dismoi('NOM_NUME_DDL', matas1, 'MATR_ASSE', repk=nu1)
    call gcncon('_', nu2)
!
    call copisd('MATR_ASSE', bas1, matas1, matas2)
    call copisd('NUME_DDL', bas1, nu1, nu2)
!
!
!     2. On corrige ce qui doit l'etre :
!     -----------------------------------
    call jedetr(matas2//'.CONL')
!
!     2.1 On corrige MATAS2.REFA :
!     -----------------------------
    call jeveuo(matas2//'.REFA', 'E', vk24=refa)
    refa(2) = nu2
    ASSERT(refa(8) .eq. 'ASSE')
    lsym = (refa(9) .eq. 'MS')
    if (refa(10) .ne. 'NOEU') call utmess('F', 'ELIMLAGR_6')
    ASSERT(refa(11) .eq. 'MPI_COMPLET')
    call settco(matas2, 'MATR_ASSE_ELIM_R')
!
!     -- la matrice MATAS2 n'est pas concernée par ELIM_LAGR :
    refa(19) = ' '
    refa(20) = matas1
!
!
!     2.2 On calcule MATAS2.VALM, NU2.NUME.SMDI et NU2.SMOS.SMHC :
!     -------------------------------------------------------------
!
!     -- calcul de :
!       neq2 : nombre de ddls de MATAS2
!       nnz2 : nombre de termes dans .VALM(1)
    call MatGetSize(elg_context(ke)%kproj, mm, nn, ierr)
    ASSERT(ierr == 0)
    neq2 = to_aster_int(mm)
    if (neq2 .eq. 0) call utmess('F', 'ELIMLAGR_7')
!
!     -- on parcourt la matrice Kproj pour repérer ses termes non nuls
    allocate (irow(neq2))
    allocate (vrow(neq2))
    call wkvect('&&ELG_CALC_MATM_RED.NZCO', 'V V I', neq2, jnzcol)
    ndiag = 0
    nnz2 = 0

    do ilig = 0, neq2-1
        call MatGetRow(elg_context(ke)%kproj, to_petsc_int(ilig), nterm, irow(1), vrow(1), &
                       ierr)
        ASSERT(ierr == 0)
        do k = 1, nterm
            jcol = irow(k)
            if (jcol .eq. ilig) ndiag = ndiag+1
            if (jcol .ge. ilig) then
                nnz2 = nnz2+1
                zi(jnzcol-1+jcol+1) = zi(jnzcol-1+jcol+1)+1
            end if
        end do
        call MatRestoreRow(elg_context(ke)%kproj, to_petsc_int(ilig), nterm, irow(1), vrow(1), &
                           ierr)
        ASSERT(ierr == 0)
    end do
    ASSERT(ndiag .eq. neq2)
!
!
!     -- allocation de MATAS2.VALM :
!     ------------------------------
    call jedetr(matas2//'.VALM')
    if (lsym) then
!       cas symétrique
        call jecrec(matas2//'.VALM', bas1//' V R', 'NU', 'DISPERSE', 'CONSTANT', 1_8)
        call jecroc(jexnum(matas2//'.VALM', 1_8))
        call jeecra(matas2//'.VALM', 'LONMAX', nnz2, kbid)
        call jeveuo(jexnum(matas2//'.VALM', 1_8), 'E', jvalm)
    else
!       cas non-symétrique
        call jecrec(matas2//'.VALM', bas1//' V R', 'NU', 'DISPERSE', 'CONSTANT', 2_8)
        call jeecra(matas2//'.VALM', 'LONMAX', nnz2, kbid)
        call jecroc(jexnum(matas2//'.VALM', 1_8))
        call jeveuo(jexnum(matas2//'.VALM', 1_8), 'E', jvalm)
        call jecroc(jexnum(matas2//'.VALM', 2_8))
        call jeveuo(jexnum(matas2//'.VALM', 2_8), 'E', jvalm2)
    end if
!
!     -- calcul de NU2.SMOS.SMDI :
!     ----------------------------
    call jedetr(nu2//'.SMOS.SMDI')
    call wkvect(nu2//'.SMOS.SMDI', bas1//' V I', neq2, jsmdi2)
    ico = 0
    do ilig = 1, neq2
        ico = ico+zi(jnzcol-1+ilig)
        zi(jsmdi2-1+ilig) = ico
    end do
!
!     -- calcul de NU2.SMOS.SMHC et remplissage de matas2.VALM  :
!     -----------------------------------------------------------
    call jedetr(nu2//'.SMOS.SMHC')
    call wkvect(nu2//'.SMOS.SMHC', bas1//' V S', nnz2, jsmhc2)
    call jerazo('&&ELG_CALC_MATM_RED.NZCO', neq2, 1_8)
    iterm = 0
    do ilig = 0, neq2-1
        call MatGetRow(elg_context(ke)%kproj, to_petsc_int(ilig), nterm, irow(1), vrow(1), &
                       ierr)
        ASSERT(ierr == 0)
        do k = 1, nterm
            jcol = irow(k)
!           partie triangulaire supérieure dans tous les cas
            if (jcol .ge. ilig) then
                n1 = zi(jnzcol-1+jcol+1)+1
                if (jcol .eq. 0) then
                    kdiag = 0
                else
                    kdiag = zi(jsmdi2-1+jcol)
                end if
                zi4(jsmhc2-1+kdiag+n1) = to_petsc_int(ilig+1)
                zr(jvalm-1+kdiag+n1) = vrow(k)
                zi(jnzcol-1+jcol+1) = n1
            end if
!           partie triangulaire inférieure dans tous le cas non-symétrique
            if (jcol .le. ilig .and. .not. lsym) then
                zr(jvalm2+iterm) = vrow(k)
                iterm = iterm+1
            end if
        end do
        call MatRestoreRow(elg_context(ke)%kproj, to_petsc_int(ilig), nterm, irow(1), vrow(1), &
                           ierr)
        ASSERT(ierr == 0)
    end do
!
!
!     2.3  Correction des autres objets de NU2 :
!     --------------------------------------------
    call jeveuo(nu2//'.NUME.REFN', 'E', j1)
    ma = zk24(j1-1+1) (1:8)
    nomgds = zk24(j1-1+2) (1:8)
    call dismoi('NB_EC', nomgds, 'GRANDEUR', repi=nec)
!
!     nu2.REFN :
!     ----------
    call jeveuo(nu2//'.NUME.REFN', 'E', j1)
    zk24(j1-1+4) = 'ELIM_LAGR'
!
!     nu2.NEQU :
!     ----------
    call jeveuo(nu2//'.NUME.NEQU', 'E', j1)
    zi(j1-1+1) = neq2
!
!     nu2.SMDE :
!     ----------
    call jeveuo(nu2//'.SMOS.SMDE', 'E', j1)
    zi(j1-1+1) = neq2
    zi(j1-1+2) = nnz2
!
!     nu2.DELG :
!     ----------
    call jedetr(nu2//'.NUME.DELG')
    call wkvect(nu2//'.NUME.DELG', bas1//' V I', neq2, j1)
    do k = 1, neq2
        zi(j1-1+k) = 0
    end do
!
!     nu2.NUEQ :
!     -----------
    call jedetr(nu2//'.NUME.NUEQ')
    call wkvect(nu2//'.NUME.NUEQ', bas1//' V I', neq2, j1)
    do k = 1, neq2
        zi(j1-1+k) = k
    end do
!
!     nu2.LILI :
!     ----------
    call jedetr(nu2//'.NUME.LILI')
    call jecreo(nu2//'.NUME.LILI', bas1//' N K24')
    call jeecra(nu2//'.NUME.LILI', 'NOMMAX', 2_8, '  ')
    call jecroc(jexnom(nu2//'.NUME.LILI', '&MAILLA'))
!
!     -- pour éviter de se planter dans vpstor.f (ligne 160)
!        qui fait un dismoi('nom_modele',...) indésirable :
    call jenuno(jexnum(nu1//'.NUME.LILI', 2_8), ligrmo)
    call jecroc(jexnom(nu2//'.NUME.LILI', ligrmo))
!
!
!     nu2.DEEQ :
!     ----------
    call jeveuo(nu1//'.NUME.DEEQ', 'L', vi=deeq)
    call jedetr(nu2//'.NUME.DEEQ')
    call wkvect(nu2//'.NUME.DEEQ', bas1//' V I', 2*neq2, jdeeq2)
!
!     nu2.PRNO (calculé à partir de .DEEQ):
!     -------------------------------------
    call jedetr(nu2//'.NUME.PRNO')
    call jecrec(nu2//'.NUME.PRNO', bas1//' V I ', 'NU', 'CONTIG', 'VARIABLE', 1_8)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
    call dismoi('NB_NL_MAILLA', ma, 'MAILLAGE', repi=nbnl)
    ASSERT(nbnl .eq. 0)
    nbnom = nbno+nbnl
    call jeecra(jexnum(nu2//'.NUME.PRNO', 1_8), 'LONMAX', nbnom*(nec+2), kbid)
!
    deallocate (irow, vrow)
    call jedetr('&&ELG_CALC_MATM_RED.NZCO')
    call jedema()
#else
    character(len=1) :: kdummy
    kdummy = matas1(1:1)//matas2(1:1)//bas1(1:1)
    call utmess('F', 'ELIMLAGR_1')
#endif
!
end subroutine
