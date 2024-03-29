! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine matass2petsc(matasz, petscMatz, iret)
!
!
! person_in_charge: natacha.bereux at edf.fr
!
#include "asterf_types.h"
#ifdef ASTER_HAVE_PETSC
#include "asterf_petsc.h"
#endif
    use aster_petsc_module
    use petsc_data_module
!
    implicit none
!
#include "asterf.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/apetsc.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/crsvfm.h"
#include "asterfort/jemarq.h"
#include "asterfort/detrsd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
!
! arguments
    character(len=*), intent(in) :: matasz
#ifdef ASTER_HAVE_PETSC
    PetscErrorCode, intent(out) :: iret
    Mat, intent(out) :: petscMatz
#else
    integer, intent(out) :: petscMatz, iret
#endif
!
!-----------------------------------------------------------------------
! BUT : ROUTINE D'INTERFACE ENTRE CODE_ASTER ET LA BIBLIOTHEQUE PETSC
!       DE RESOLUTION DE SYSTEMES LINEAIRES.
!
! IN  : MATASS   (K19) : NOM DE LA MATR_ASSE
!                       (SI ACTION=PRERES/RESOUD)
! OUT : petscMatz   (I) : MATRICE PETSC
! OUT : IRET       (I) : CODE_RETOUR
!-----------------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!----------------------------------------------------------------
!
!     VARIABLES LOCALES
    integer :: k, ierror
    integer ::  ibid
    real(kind=8) :: rbid
!
    character(len=24), dimension(:), pointer :: slvk => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=19) :: solvbd, matas
!
!----------------------------------------------------------------
!
!     Variables PETSc
    PetscErrorCode :: ierr

    call jemarq()
!
    matas = matasz
    rbid = 0.d0

!   -- Creation d'un solveur bidon
    solvbd = '&&MAT2PET'
   call crsvfm(solvbd, matas, 'D', rank='L', pcpiv=50, usersmz='IN_CORE', blreps=rbid, renumz=' ', &
                redmpi=-9999)
    call jeveuo(solvbd//'.SLVK', 'L', vk24=slvk)
    slvk(2) = 'SANS'

!   -- Effacement si déjà factorisée
    call jeveuo(matas//'.REFA', 'E', vk24=refa)
    refa(8) = ' '

!   -- Conversion de matass vers petsc
    call apetsc('PRERES', solvbd, matas, [0.d0], ' ', &
                0, ibid, ierror)
    k = get_mat_id(matas)
    call MatDuplicate(ap(k), MAT_COPY_VALUES, petscMatz, ierr)
    ASSERT(ierr .eq. 0)

!   -- Nettoyage
!   Destruction des objets petsc
    call apetsc('DETR_MAT', solvbd, matas, [0.d0], ' ', &
                0, ibid, ierror)
!   Destruction du solveur bidon
    call detrsd('SOLVEUR', solvbd)
    iret = 0

    call jedema()
#else
    petscMatz = 0
    iret = 0
    call utmess('F', 'FERMETUR_10')
#endif
end subroutine
