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
!
subroutine apetsc(action, solvez, matasz, rsolu, vcinez, &
                  nbsol, istop, iret)
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use petsc_data_module
    use elg_module
!
    implicit none
!
#include "jeveux.h"
#include "asterc/aster_petsc_initialize.h"
#include "asterc/aster_petsc_finalize.h"
#include "asterfort/apmain.h"
#include "asterfort/apldlt.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mtmchc.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dcopy.h"
!
    character(len=*) :: action, solvez, matasz, vcinez
    real(kind=8) :: rsolu(*)
    integer(kind=8) :: nbsol, istop, iret
!-----------------------------------------------------------------------
! BUT : ROUTINE D'INTERFACE ENTRE CODE_ASTER ET LA BIBLIOTHEQUE PETSC
!       DE RESOLUTION DE SYSTEMES LINEAIRES.
!
! IN  : ACTION
!     /'DETR_MAT': POUR DETRUIRE L'INSTANCE PETSC ASSOCIEE A UNE MATRICE
!     /'PRERES'  : POUR CONSTRUIRE LE PRECONDITIONNEUR
!                 (ATTENTION EN // LA CONSTRUCTION DE CERTAINS PC EST
!                  RETARDEE)
!     /'RESOUD'  : POUR RESOUDRE LE SYSTEME LINEAIRE
!     /'ELIM_LAGR'  : CALCULE LES MATRICES NECESSAIRES A
!                     LA FONCTIONNALITE ELIM_LAGR='OUI'
!     /'FIN'     : POUR FERMER DEFINITIVEMENT PETSC
!                  NECESSAIRE POUR DECLENCHER L'AFFICHAGE DU PROFILING
!
! IN  : SOLVEU   (K19) : NOM DE LA SD SOLVEUR
!                       (SI ACTION=PRERES/ELIM_LAGR[+/-]R)
! IN  : MATASS   (K19) : NOM DE LA MATR_ASSE
!                       (SI ACTION=PRERES/RESOUD)
! I/O : RSOLU      (R) : EN ENTREE : VECTEUR SECOND MEMBRE (REEL)
!                        EN SORTIE : VECTEUR SOLUTION      (REEL)
!                       (SI ACTION=RESOUD)
! IN  : VCINE    (K19) : NOM DU CHAM_NO DE CHARGEMENT CINEMATIQUE
!                       (SI ACTION=RESOUD)
! IN  : NBSOL      (I) : NOMBRE DE SYSTEMES A RESOUDRE
! IN  : ISTOP      (I) : COMPORTEMENT EN CAS D'ERREUR
! OUT : IRET       (I) : CODE_RETOUR
!                      /  0 : OK
!                      /  1 : NOMBRE MAX D'ITERATIONS ATTEINT
!-----------------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!     VARIABLES LOCALES
    integer(kind=8) :: iprem, k, nglo, kdeb, jnequ, ier2
    integer(kind=8) :: kptsc, lslvo
    integer(kind=8) :: np, i
    real(kind=8) :: r8
!
    logical :: mat_not_recorded
    character(len=19) :: solveu, matas, vcine, kbid
    character(len=14) :: nu
    character(len=4) :: etamat
    character(len=1) :: rouc
    real(kind=8), pointer :: travail(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=80), pointer :: slvo(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    character(len=2500) :: myopt
    !
!----------------------------------------------------------------
!
!     Variables PETSc
    PetscErrorCode :: ierr
    PetscScalar :: sbid
    PetscOffset :: offbid
    PetscReal :: rbid
    PetscBool :: initialized
    blas_int :: b_incx, b_incy, b_n
!
!----------------------------------------------------------------
!   INITIALISATION DE PETSC A FAIRE AU PREMIER APPEL
    save iprem
    data iprem/0/
!----------------------------------------------------------------
    call jemarq()
!
    solveu = solvez
    matas = matasz
    vcine = vcinez
    iret = 0
!
!
!   0. FERMETURE DE PETSC DANS FIN
!   ------------------------------
    if (action .eq. 'FIN') then
!       petsc a-t-il ete initialise ?
        if (iprem .eq. 1) then
            call aster_petsc_finalize()
!           on ne verifie pas le code retour car on peut
!           se retrouver dans fin suite a une erreur dans l'initialisation
            iprem = 0
        end if
        goto 999
    end if
!
!   Vérification du nom de la matrice
    ASSERT(matas .ne. ' ')
!
!   PETSc a-t-il ete initialise ?
    call PetscInitialized(initialized, ierr)
    ASSERT(ierr .eq. 0)
!
!   Récupération des options PETSc
    myopt = ''
    if (action .ne. 'DETR_MAT') then
        call jeveuo(solveu//'.SLVO', 'L', vk80=slvo)
        call jeveuo(solveu//'.SLVI', 'L', vi=slvi)
        lslvo = slvi(9)
        do i = 1, lslvo
            myopt(80*(i-1)+1:80*i) = slvo(i)
        end do
    end if
!
    if (iprem .eq. 0 .or. .not. initialized) then
!     --------------------
!        -- quelques verifications sur la coherence Aster / Petsc :
        ASSERT(kind(rbid) .eq. kind(r8))
        ASSERT(kind(sbid) .eq. kind(r8))
        ASSERT(kind(offbid) .eq. kind(np))
!
        ier2 = 0
        call aster_petsc_initialize(myopt, ier2)
        ierr = to_petsc_int(ier2)
        if (ierr .ne. 0) call utmess('F', 'PETSC_1')
        do k = 1, nmxins
            ap(k) = PETSC_NULL_MAT
            kp(k) = PETSC_NULL_KSP
            user_ksp(k) = ASTER_FALSE
            nomats(k) = ' '
            nosols(k) = ' '
            nonus(k) = ' '
            tblocs(k) = -1
            options(k) = ' '
        end do
        xlocal = PETSC_NULL_VEC
        xglobal = PETSC_NULL_VEC
        xscatt = PETSC_NULL_VECSCATTER
        spsomu = ' '
        spmat = ' '
        spsolv = ' '
        iprem = 1
    end if
!
!   1. On ne veut pas de matrice complexe :
!   ----------------------------------------
    call jelira(matas//'.VALM', 'TYPE', cval=rouc)
    if (rouc .ne. 'R') call utmess('F', 'PETSC_2')
!   nglo est le nombre total de degrés de liberté
    call dismoi('NOM_NUME_DDL', matas, 'MATR_ASSE', repk=nu)
    call jeveuo(nu//'.NUME.NEQU', 'L', jnequ)
    nglo = zi(jnequ)
!
!  2. On recherche l'identifiant de l'image PETSc
!  de la matrice matas
!
    kptsc = get_mat_id(matas)
    mat_not_recorded = (kptsc == 0)
!
!
    if (action .eq. 'DETR_MAT') then
        if (mat_not_recorded) then
! On n'a pas cree d'image PETSc de la matrice => rien à detruire !
        else
! L'image PETSc de la a matrice est stockée dans le tableau ap,
! a l'indice kptsc => on la détruit
            kbid = repeat(" ", 19)
            call apmain(action, kptsc, [0.d0], kbid, 0_8, &
                        iret)
        end if
        goto 999
    end if
!
!   3. Quelques verifications et petites actions :
!   ----------------------------------------------
!
    if (action .eq. 'PRERES') then
        call mat_record(matas, solveu, kptsc, myopt)
    else if (action .eq. 'RESOUD') then
        kptsc = get_mat_id(matas)
        ASSERT(nbsol .ge. 1)
        ASSERT((istop .eq. 0) .or. (istop .eq. 2))
!
    else if (action .eq. 'ELIM_LAGR') then
        call mat_record(matas, solveu, kptsc)
    end if
!
!   4. Si LDLT_INC, il faut renumeroter la matrice (RCMK) :
!   --------------------------------------------------------
!
    call apldlt(kptsc, action, 'PRE', rsolu, vcine, &
                nbsol)
!
!   5. Verifications + elimination des ddls (affe_char_cine)
!   --------------------------------------------------------
    call jeveuo(nomat_courant//'.REFA', 'E', vk24=refa)
    if (action .eq. 'PRERES' .or. nomat_courant .eq. '&&apldlt.matr') then
!
!       -- Verification que la matrice n'a pas deja ete factorisee
        etamat = refa(8) (1:4)
        if (etamat .eq. 'DECT') then
            call utmess('A', 'PETSC_4')
            goto 999
        else
            refa(8) = 'DECT'
        end if
!
!        -- elimination des ddls (affe_char_cine)
!        ASSERT(refa(3).ne.'ELIMF')
        if (refa(3) .eq. 'ELIML') call mtmchc(nomat_courant, 'ELIMF')
        ASSERT(refa(3) .ne. 'ELIML')
!
    else if (action .eq. 'ELIM_LAGR') then
        call build_elg_context(nomat_courant)
        iret = 0
        goto 999
    end if
!
!   5. APPEL DE PETSC :
!   -------------------
    if (action .eq. 'RESOUD') then
        AS_ALLOCATE(vr=travail, size=nglo)
        do k = 1, nbsol
            kdeb = (k-1)*nglo+1
            b_n = to_blas_int(nglo)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rsolu(kdeb), b_incx, travail, b_incy)
            call apmain(action, kptsc, travail, vcine, istop, &
                        iret)
            b_n = to_blas_int(nglo)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, travail, b_incx, rsolu(kdeb), b_incy)
        end do
        AS_DEALLOCATE(vr=travail)
    else
        call apmain(action, kptsc, rsolu, vcine, istop, &
                    iret)
    end if
!
!   6. Si LDLT_INC, il faut revenir a la numerotation initiale :
!   -------------------------------------------------------------
    call apldlt(kptsc, action, 'POST', rsolu, vcine, &
                nbsol)
!
999 continue
    call jedema()
#else
    character(len=1) :: kdummy
    real(kind=8) :: rdummy
    integer(kind=8) :: idummy
    idummy = nbsol+istop+iret
    rdummy = rsolu(1)
    kdummy = action//solvez//matasz//vcinez
    call utmess('F', 'FERMETUR_10')
#endif
end subroutine
