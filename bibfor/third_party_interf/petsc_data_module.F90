! --------------------------------------------------------------------
! Copyright (C) 2016 - 2025 - EDF R&D - www.code-aster.org
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

module petsc_data_module
!
#include "asterf_petsc.h"
    use aster_petsc_module
!
    implicit none
! aslint:disable=
    private
#ifdef ASTER_HAVE_PETSC
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!--------------------------------------------------------------------
!
!
! person_in_charge natacha.bereux at edf.fr
!
! Variables globales PETSc pour la définition d'un système linéaire
! au format PETSc
!
! On ne peut pas creer directement un tableau de pointer,
! il faut passer par un artifice (type derive) :

    type p_int4
        sequence
        integer(kind=4), pointer :: pi4(:)
    end type

    integer(kind=8), parameter, public :: nmxins = 5
    character(len=19), public  :: nomats(nmxins), nosols(nmxins), nomat_courant
    character(len=14), public  :: nonus(nmxins), nonu_courant
    character(len=2500), public  :: options(nmxins)
    Mat, public :: ap(nmxins)
    KSP, public :: kp(nmxins)
    Vec, public :: b, x
    aster_logical, public :: user_ksp(nmxins)
! Les variables suivantes sont utilisees par les preconditionneurs multigrille
    integer(kind=4), public :: tblocs(nmxins)
!
!----------------------------------------------------------------
! Variables globales pour la définition d'un preconditionneur
! simple precision ldlt_sp
    character(len=19), public :: spsomu, spmat, spsolv
    Vec, public :: xlocal, xglobal
    VecScatter, public :: xscatt
!----------------------------------------------------------------
!
    public :: get_mat_id, mat_record
!
contains
!
! Cette fonction renvoie l'identifiant (kptsc) permettant d'accéder
! à la matrice au format PETSc créée à partir de la matr_asse de nom
! matas.
! Si on n'a pas créé de matrice PETSc, la fonction renvoie 0.
!
    function get_mat_id(matas) result(kptsc)
        !
        ! Dummy arguments
        character(len=19), intent(in) :: matas
        integer(kind=8) :: kptsc
        ! Local variables
        character(len=14) :: nu
        integer(kind=8) :: jnequ, nglo, k
        PetscInt :: m, n
        PetscErrorCode :: ierr
        !
        call jemarq()
!
!   On cherche si la matrice est deja enregistree :
!   -------------------------------------------------
!   On teste le nom de la matrice, celui du nume_ddl,
!   et la taille des matrices aster et petsc
        call dismoi('NOM_NUME_DDL', matas, 'MATR_ASSE', repk=nu)
        call jeveuo(nu//'.NUME.NEQU', 'L', jnequ)
        nglo = zi(jnequ)
! Valeur par défaut de kptsc
        kptsc = 0
!
        do k = 1, nmxins
            if ((nomats(k) .eq. matas) .and. (nonus(k) .eq. nu)) then
! si de plus le clone PETSc a ete cree, on verifie que les dimensions
! des matrices aster et petsc sont coherentes
                if (ap(k) .ne. PETSC_NULL_MAT) then
                    call MatGetSize(ap(k), m, n, ierr)
                    ASSERT(ierr .eq. 0)
                    ASSERT(m .eq. n)
                    ASSERT(nglo .le. n)
                end if
! la verification a ete effectuee avec succes, on renvoie k
                kptsc = k
            end if
        end do
        call jedema()
        !
    end function get_mat_id
!
!
! Retourne un identifiant libre pour stocker une nouvelle
! matrice. Si on ne trouve pas d'identifiant, on renvoie 0
    function get_new_mat_id() result(kptsc)
        !
        ! Dummy arguments
        integer(kind=8) :: kptsc
        ! Local variables
        integer(kind=8) :: k
!
!   Y-a-t-il encore une place libre ? Calcul de kptsc :
!   ---------------------------------------------------
        kptsc = 0
        do k = 1, nmxins
            if (nomats(k) .eq. ' ') then
                kptsc = k
                exit
            end if
        end do
    end function get_new_mat_id
!
! La routine mat_record enregistre la matrice matas
! i.e. determine son identifiant kptsc
! et note son nom dans les tableaux
    subroutine mat_record(matas, solveu, kptsc, user_opt)
        ! Dummy arguments
        character(len=19), intent(in) :: matas, solveu
        integer(kind=8), intent(out)          :: kptsc
        character(len=2500), intent(in), optional :: user_opt
        ! Local variables
        character(len=19) :: nu
        !
        call dismoi('NOM_NUME_DDL', matas, 'MATR_ASSE', repk=nu)
        !
        ! Verification : est-ce que la matrice est deja enregistree ?
        kptsc = get_mat_id(matas)
        if (kptsc == 0) then
            kptsc = get_new_mat_id()
            if (kptsc == 0) then
                call utmess('F', 'PETSC_3')
            end if
            !
            ASSERT(nomats(kptsc) .eq. ' ')
            ASSERT(nosols(kptsc) .eq. ' ')
            ASSERT(nonus(kptsc) .eq. ' ')
            ASSERT(options(kptsc) .eq. ' ')
            !
            nomats(kptsc) = matas
            nonus(kptsc) = nu(1:14)
            nosols(kptsc) = solveu
            if (present(user_opt)) then
                options(kptsc) = user_opt
            end if
        end if
!
    end subroutine mat_record
#endif
end module petsc_data_module
