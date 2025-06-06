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

subroutine elg_calc_solu(matas1, nsecm, rsolu2, rsolu1, omega2, ke_mass)
#include "asterf_types.h"
#include "asterf_petsc.h"
!
    use aster_petsc_module
    use elg_data_module
    implicit none
! person_in_charge: natacha.bereux at edf.fr
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elg_allocvr.h"
#include "asterfort/elg_calcxl.h"
#include "asterfort/elg_calcxl_modal.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nudlg2.h"
#include "asterfort/utmess.h"
!
    character(len=19) :: matas1
    integer(kind=8) :: nsecm
    real(kind=8) :: rsolu2(*), rsolu1(*)
    real(kind=8), optional, intent(in)::omega2
#ifdef ASTER_HAVE_PETSC
    PetscInt, optional, intent(in):: ke_mass
!--------------------------------------------------------------
! BUT :
!   calculer les solutions complètes (RSOLU1) correspondant aux
!   solutions réduites (RSOLU2)
!
! IN  : MATAS1 : sd_matr_asse avec ses conditions dualisées
!                à eliminer
! IN  : NSECM  :  nombre de solutions
! IN  : RSOLU2(*)  : vecteur de réels de dimension nsecm*neq2
!                  (valeurs des solutions réduites)
! IN  : MODAL_  :  vrai si on reconstruit un calcul modal
! OUT : RSOLU1(*)  : vecteur de réels de dimension nsecm*neq1
!                  (valeurs des solutions complètes)
!---------------------------------------------------------------
!
!
!================================================================
    character(len=14) :: nu1, nu2
    character(len=19) :: matas2
    real(kind=8) :: val
    integer(kind=8) :: neq1, neq2, ico, ieq2
    integer(kind=8) :: k1, k2
    real(kind=8), pointer :: conl(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
    integer(kind=8), pointer :: dlg2(:) => null()
    character(len=24), pointer :: refa(:) => null()
    aster_logical::modal
    PetscErrorCode :: ierr
    PetscInt :: n1, n2, n3
    PetscScalar :: xx(1), p1
    PetscOffset :: xidx
    Vec :: x1, y, vlag, tmp1
!
!
!----------------------------------------------------------------
    call jemarq()
!
    modal = .false.
    if (present(omega2)) then
        modal = .true.
        ASSERT(present(ke_mass))
    end if
!
    call jeveuo(matas1//'.REFA', 'L', vk24=refa)
    matas2 = refa(19) (1:19)
    ASSERT(matas2 .ne. ' ')
!
    call dismoi('NOM_NUME_DDL', matas1, 'MATR_ASSE', repk=nu1)
    call dismoi('NOM_NUME_DDL', matas2, 'MATR_ASSE', repk=nu2)
    call dismoi('NB_EQUA', nu1, 'NUME_DDL', repi=neq1)
    call dismoi('NB_EQUA', nu2, 'NUME_DDL', repi=neq2)
    call jeveuo(nu1//'.NUME.DELG', 'L', vi=delg)
    call nudlg2(nu1)
    call jeveuo(nu1//'.NUME.DLG2', 'L', vi=dlg2)
!
! à faire ....
    ASSERT(nsecm .eq. 1)
!
!     -- dimensions n1, n2, n3 :
    call MatGetSize(elg_context(ke)%tfinal, n1, n3, ierr)
    ASSERT(ierr == 0)
    call MatGetSize(elg_context(ke)%matc, n2, n1, ierr)
    ASSERT(ierr == 0)
    ASSERT(neq2 .eq. n3)
    ASSERT(neq1 .eq. n1+2*n2)
!
!
!     allocation et remplissage de Y = RSOLU2
    call elg_allocvr(y, to_aster_int(n3))
    call VecGetArray(y, xx, xidx, ierr)
    ASSERT(ierr == 0)
    do ieq2 = 1, neq2
        xx(xidx+ieq2) = rsolu2(ieq2)
    end do
    call VecRestoreArray(y, xx, xidx, ierr)
    ASSERT(ierr == 0)
!
!
!     Calcul de TMP1 =  T*Y :
    call elg_allocvr(tmp1, to_aster_int(n1))
    call MatMult(elg_context(ke)%tfinal, y, tmp1, ierr)
    ASSERT(ierr == 0)
!
!     Calcul de X1= x0 + T*Y :
    call elg_allocvr(x1, int(n1, 8))
    call VecCopy(elg_context(ke)%vx0, x1, ierr)
    ASSERT(ierr == 0)
    p1 = 1.d0
    call VecAXPY(x1, p1, tmp1, ierr)
    ASSERT(ierr == 0)
!
!     -- on recopie X1 dans RSOLU1 :
    call VecGetArray(x1, xx, xidx, ierr)
    ASSERT(ierr == 0)
    ico = 0
    do k1 = 1, neq1
        if (delg(k1) .eq. 0) then
            ico = ico+1
            rsolu1(k1) = xx(xidx+ico)
        end if
    end do
    ASSERT(ico .eq. n1)
    call VecRestoreArray(x1, xx, xidx, ierr)
    ASSERT(ierr == 0)
!
!
!     calcul des coefficients de Lagrange :
    call elg_allocvr(vlag, to_aster_int(n2))
    if (modal) then
        call elg_calcxl_modal(x1, omega2, ke_mass, vlag)
    else
        call elg_calcxl(x1, vlag)
    end if
    !
!     -- on recopie VLAG dans RSOLU1 :
!        remarque : il faut diviser VLAG par 2 (2 lagranges)
!                   Lagrange "1"  et "2" :
    call jeveuo(matas1//'.CONL', 'L', vr=conl)
    call VecGetArray(vlag, xx, xidx, ierr)
    ASSERT(ierr == 0)
    ico = 0
    do k1 = 1, neq1
        if (delg(k1) .eq. -1) then
            ico = ico+1
            if (modal) then
                val = xx(xidx+ico)/2.d0
            else
                val = xx(xidx+ico)*conl(k1)/2.d0
            end if
            rsolu1(k1) = val
! k2 lagrange "2" associé au lagrange "1" k1
            k2 = dlg2(k1)
            ASSERT(k2 .gt. 0)
            rsolu1(k2) = val
        end if
    end do
    ASSERT(ico .eq. n2)
    call VecRestoreArray(vlag, xx, xidx, ierr)
    ASSERT(ierr == 0)
!
!     -- ménage :
    call VecDestroy(y, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(tmp1, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(x1, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(tmp1, ierr)
    ASSERT(ierr == 0)
!
    call jedema()
#else
    character(len=1) :: kdummy
    integer(kind=8) :: idummy
    integer(kind=8), optional, intent(in):: ke_mass
    real(kind=8) :: rdummy
    kdummy = matas1(1:1)
    idummy = nsecm+ke_mass
    rdummy = rsolu1(1)+rsolu2(1)+omega2
    call utmess('F', 'ELIMLAGR_1')
#endif
!
end subroutine
