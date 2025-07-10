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

subroutine elg_calcxl_modal(x1, omega2, ke_mass, vlag)
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use elg_data_module

    implicit none
! person_in_charge: nicolas.tardieu at edf.fr
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
!----------------------------------------------------------------
! Calcul des coefficients de Lagrange pour ELIM_LAGR='OUI'
!     IN  : Vec X1    (solution sur les dds physiques)
!
!     OUT : Vec VLAG  (solution sur les dds Lagrange "1")
!
!     On résout  vlag = A'\ Y de la façon suivante :
!       * calcul de CCT = A*A'
!       * calcul de AY = A*Y
!       * résolution (CG) vlag = A*A' \ AY
!
!     Rq: La méthode originale utilise une factorisation QR
!     préalable de A. Le calcul des multiplicateurs s'effectue alors
!     par des descentes remontées:
!      Les coefficients de Lagrange (1 et 2) sont
!     calculés par :
!       L = (R'*R) \ A*(b - B*x)
!----------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
    Vec :: x1, vlag
    real(kind=8) :: omega2
    PetscInt :: ke_mass
!
!================================================================
    mpi_int :: mpicomm, nbproc, rang
    integer(kind=8) :: ifm, niv
    aster_logical :: info
    PetscInt :: n1, n2, n3
    KSPConvergedReason :: reason
    PetscErrorCode :: ierr
    Vec :: mx, kx, bkx
    PetscReal :: aster_petsc_real

!----------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
    info = niv .eq. 2
    !
!   -- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET_WORLD', mpicomm)
    call asmpi_info(rank=rang, size=nbproc)
!
    aster_petsc_real = PETSC_DEFAULT_REAL
!
!     -- dimensions :
!       n1 : # ddls physiques
!       n2 : # lagranges "1"
!     ----------------------------------------------------------
    call MatGetSize(elg_context(ke)%matc, n2, n1, ierr)
!   le système est-il bien  sur-déterminé ?
    ASSERT(n1 > n2)
!
!     -- vérifs :
    call VecGetSize(x1, n3, ierr)
    ASSERT(n3 .eq. n1)
    call VecGetSize(vlag, n3, ierr)
    ASSERT(n3 .eq. n2)
!
!
!     -- calcul de MX = M*x :
    call VecDuplicate(x1, mx, ierr)
    ASSERT(ierr == 0)
    call MatMult(elg_context(ke_mass)%matb, x1, mx, ierr)
    ASSERT(ierr == 0)

!     -- calcul de KX = K*x :
    call VecDuplicate(x1, kx, ierr)
    ASSERT(ierr == 0)
    call MatMult(elg_context(ke)%matb, x1, kx, ierr)
    ASSERT(ierr == 0)

!     -- calcul de KX = KX - omega2*MX :
    call VecAXPY(kx, -omega2, mx, ierr)
    ASSERT(ierr == 0)

!     -- calcul de BKX = B*KX :
    call VecDuplicate(elg_context(ke)%vecc, bkx, ierr)
    ASSERT(ierr == 0)
    call MatMult(elg_context(ke)%matc, kx, bkx, ierr)
    ASSERT(ierr == 0)
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! On résout :
!  * vlag = -A A'\ A(K-omega**2*M) X
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   -- Solve the linear system
    call KSPSolve(elg_context(ke)%ksp, bkx, vlag, ierr)
    ASSERT(ierr == 0)
    call VecScale(vlag, -1.D0, ierr)
    ASSERT(ierr == 0)
!
!  Check the reason why KSP solver ended
!
    call KSPGetConvergedReason(elg_context(ke)%ksp, reason, ierr)
    ASSERT(ierr == 0)
    if (reason < 0) then
        call utmess('F', 'ELIMLAGR_8')
    end if
!
!  Free work space.  All PETSc objects should be destroyed when they
!  are no longer needed.
    call VecDestroy(bkx, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(mx, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(kx, ierr)
    ASSERT(ierr == 0)
!
    call jedema()
!
#else
    integer(kind=8) :: x1, vlag, ke_mass, omega2
    vlag = x1
    ASSERT(.false.)
#endif
!
end subroutine
