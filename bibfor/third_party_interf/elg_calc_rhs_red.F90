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

subroutine elg_calc_rhs_red(matas1, nsecm, secm, solu2)
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
#include "asterfort/elg_calcx0.h"
#include "asterfort/elg_allocvr.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: matas1
    integer(kind=8) :: nsecm
    real(kind=8) :: secm(*)
    character(len=24) :: solu2
!--------------------------------------------------------------
! BUT :
!   calculer le(s) second-membre(s) réduits correspondant à
!   un (ou plusieurs) second-membre(s) complets.
!
! IN  : MATAS1 : sd_matr_asse avec ses conditions dualisées
!                à eliminer
! IN  : NSECM  :  nombre de seconds membres
! IN  : SECM(*)  : vecteur de réels de dimension nsecm*neq1
!                  (valeurs de(s) second-membre(s) complets)
! IN/JXOUT : SOLU2  : objet jeveux qui contiendra le(s)
!                     seconds membre(s) réduit(s)
!---------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!================================================================
    character(len=1) :: kbid
    character(len=14) :: nu1, nu2
    character(len=19) :: matas2
    integer(kind=8) ::  ibid
    integer(kind=8) :: neq1, neq2, icob, icoc, ieq1, ieq2
    integer(kind=8) :: jsolu2
    character(len=24), pointer :: refa(:) => null()
    real(kind=8), pointer :: conl(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
    PetscErrorCode :: ierr
    PetscInt :: n1, n2, n3
    PetscScalar :: xx(1), m1
    PetscOffset :: xidxb, xidxc, xidxb2
    Vec :: bx0, vecb2, vectmp
!
!
!----------------------------------------------------------------
    call jemarq()
!
    call jeveuo(matas1//'.REFA', 'L', vk24=refa)
    matas2 = refa(19) (1:19)
    ASSERT(matas2 .ne. ' ')
!
    call dismoi('NOM_NUME_DDL', matas1, 'MATR_ASSE', repk=nu1)
    call dismoi('NOM_NUME_DDL', matas2, 'MATR_ASSE', repk=nu2)
    call dismoi('NB_EQUA', nu1, 'NUME_DDL', repi=neq1)
    call dismoi('NB_EQUA', nu2, 'NUME_DDL', repi=neq2)
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
!     -- allocation de VecB, VecC, VecB2 :
!     ---------------------------------------------
    call elg_allocvr(elg_context(ke)%vecb, int(n1, 8))
    call elg_allocvr(elg_context(ke)%vecc, int(n2, 8))
    call elg_allocvr(vecb2, int(n3, 8))
!
!
!     -- calcul de VecB et VecC (extraits de SECM) :
!     ------------------------------------------------
    call VecGetArray(elg_context(ke)%vecb, xx, xidxb, ierr)
    ASSERT(ierr == 0)
    call VecGetArray(elg_context(ke)%vecc, xx, xidxc, ierr)
    ASSERT(ierr == 0)
    call jeveuo(nu1//'.NUME.DELG', 'L', vi=delg)
    call jeveuo(matas1//'.CONL', 'L', vr=conl)
!
    icob = 0
    icoc = 0
    do ieq1 = 1, neq1
        if (delg(ieq1) .eq. 0) then
            icob = icob+1
            xx(xidxb+icob) = secm(ieq1)
        else if (delg(ieq1) .eq. -1) then
            icoc = icoc+1
            xx(xidxc+icoc) = secm(ieq1)*conl(ieq1)
        end if
    end do
    call VecRestoreArray(elg_context(ke)%vecb, xx, xidxb, ierr)
    ASSERT(ierr == 0)
    call VecRestoreArray(elg_context(ke)%vecc, xx, xidxc, ierr)
    ASSERT(ierr == 0)
!
!
!     -- calcul de Vx0 = A \ VecC
    call VecDuplicate(elg_context(ke)%vecb, elg_context(ke)%vx0, ierr)
    ASSERT(ierr == 0)
!
    call elg_calcx0()
!
!     -- calcul de BX0 = B*Vx0 :
    call VecDuplicate(elg_context(ke)%vecb, bx0, ierr)
    ASSERT(ierr == 0)
    call MatMult(elg_context(ke)%matb, elg_context(ke)%vx0, bx0, ierr)
    ASSERT(ierr == 0)
!
!
!     -- calcul de VecTmp = b - B*Vx0 :
    m1 = -1.d0
    call VecDuplicate(elg_context(ke)%vecb, vectmp, ierr)
    ASSERT(ierr == 0)
    call VecCopy(elg_context(ke)%vecb, vectmp, ierr)
    ASSERT(ierr == 0)
    call VecAXPY(vectmp, m1, bx0, ierr)
    ASSERT(ierr == 0)
!
!     -- calcul de VecB2 = T'*(b - B*Vx0) :
    call MatMultTranspose(elg_context(ke)%tfinal, vectmp, vecb2, ierr)
    ASSERT(ierr == 0)
!
!     -- recopie de VecB2 dans SOLU2 :
    call wkvect(solu2, 'V V R', neq2, jsolu2)
    call VecGetArray(vecb2, xx, xidxb2, ierr)
    ASSERT(ierr == 0)
    do ieq2 = 1, neq2
        zr(jsolu2-1+ieq2) = xx(xidxb2+ieq2)
    end do
    call VecRestoreArray(vecb2, xx, xidxb2, ierr)
    ASSERT(ierr == 0)
!
!
!
!     -- ménage :
    call VecDestroy(vectmp, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(bx0, ierr)
    ASSERT(ierr == 0)
    call VecDestroy(vecb2, ierr)
    ASSERT(ierr == 0)
!
    call jedema()
#else
    character(len=1) :: kdummy
    integer(kind=8) :: idummy
    real(kind=8) :: rdummy
    kdummy = matas1(1:1)//solu2(1:1)
    idummy = nsecm
    rdummy = secm(1)
    call utmess('F', 'ELIMLAGR_1')
#endif
!
end subroutine
