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

subroutine ldlt_matr(matas1, matas2, kperm, basp)
    implicit none
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jelira.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/ldlt_renum.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/jeexin.h"
!-----------------------------------------------------------------------
!  FONCTION  :
!  Creation de la matr_asse matas2 par copie de matas1 en renumerotant
!  les lignes et colonnes pour que LDLT soit plus efficace.
!
!  Retourne :
!  * matas2 (sd_matr_asse)
!  * l'objet jeveux kperm qui donne la permutation entre matas1
!    et matas2
!
!  Remarques :
!    * matas2 (et son nume_ddl) sont crees systematiquement sur la base volatile
!    * l'objet kperm est cree sur la base 'basp'
!-----------------------------------------------------------------------
    character(len=19), intent(in) :: matas1, matas2
    character(len=24), intent(in) :: kperm
    character(len=1), intent(in) :: basp
!-----------------------------------------------------------------------
    integer(kind=8) :: nvalm, neq, nnz, kterm1, kterm2, jcol1, jcol2, ilig1, ilig2
    integer(kind=8) :: k2, k2_1, k2_2, k, iexi
    character(len=24), pointer :: refa1(:) => null()
    character(len=24), pointer :: refa2(:) => null()
    character(len=1) :: typ1
    integer(kind=8), pointer :: smde1(:) => null()
    integer(kind=8), pointer :: smde2(:) => null()
    integer(kind=8), pointer :: perm(:) => null()
    integer(kind=8), pointer :: smdi1(:) => null()
    integer(kind=8), pointer :: smdi2(:) => null()
    integer(kind=8), pointer :: ccid(:) => null()
    integer(kind=8), pointer :: cci2(:) => null()
    real(kind=8), pointer :: conl(:) => null()
    real(kind=8), pointer :: con2(:) => null()
    integer(kind=4), pointer :: smhc1(:) => null()
    integer(kind=4), pointer :: smhc2(:) => null()
    real(kind=8), pointer :: valm1_1(:) => null()
    real(kind=8), pointer :: valm1_2(:) => null()
    real(kind=8), pointer :: valm2_1(:) => null()
    real(kind=8), pointer :: valm2_2(:) => null()
    aster_logical :: lsym, lmd, perm_diag
    character(len=14) nu1, nu2
!-----------------------------------------------------------------------
    call jemarq()

!   1. On fabrique la nouvelle numerotation (nu2) :
!   ----------------------------------------------------
    call dismoi('NOM_NUME_DDL', matas1, 'MATR_ASSE', repk=nu1)
    nu2 = '&&ldlt_matr.nu'
    call ldlt_renum(nu1, nu2, kperm, basp)
    call jeveuo(kperm, 'L', vi=perm)

!   2. On fabrique la nouvelle matrice (matas2) :
!   ----------------------------------------------------
    call jeveuo(matas1//'.REFA', 'L', vk24=refa1)
    lmd = (refa1(11) .eq. 'MATR_DISTR')
    ASSERT(.not. lmd)
    call copisd('MATR_ASSE', 'V', matas1, matas2)
    call jeveuo(matas2//'.REFA', 'E', vk24=refa2)
    refa2(2) = nu2

    call jeveuo(nu1//'.SMOS.SMDE', 'L', vi=smde1)
    call jeveuo(nu2//'.SMOS.SMDE', 'L', vi=smde2)
    neq = smde1(1)
    nnz = smde1(2)
    ASSERT(smde1(3) .eq. 1)
    ASSERT(smde2(1) .eq. neq)
    ASSERT(smde2(2) .eq. nnz)
    call jeveuo(nu1//'.SMOS.SMDI', 'L', vi=smdi1)
    call jeveuo(nu2//'.SMOS.SMDI', 'L', vi=smdi2)
    call jeveuo(nu1//'.SMOS.SMHC', 'L', vi4=smhc1)
    call jeveuo(nu2//'.SMOS.SMHC', 'L', vi4=smhc2)

!   2.1 Recopie de matas1.VALM dans matas2.VALM :
!   -------------------------------------------------
    call jelira(matas1//'.VALM', 'TYPE', cval=typ1)
    ASSERT(typ1 .eq. 'R')

    call jelira(matas1//'.VALM', 'NMAXOC', nvalm)
    lsym = nvalm .eq. 1
    call jeveuo(jexnum(matas1//'.VALM', 1), 'L', vr=valm1_1)
    call jeveuo(jexnum(matas2//'.VALM', 1), 'E', vr=valm2_1)
    valm2_1(:) = 0.d0
    if (.not. lsym) then
        call jeveuo(jexnum(matas1//'.VALM', 2), 'L', vr=valm1_2)
        call jeveuo(jexnum(matas2//'.VALM', 2), 'E', vr=valm2_2)
        valm2_2(:) = 0.d0
    end if

    jcol1 = 1
    do kterm1 = 1, nnz
        if (smdi1(jcol1) .lt. kterm1) jcol1 = jcol1+1
        ilig1 = smhc1(kterm1)

        jcol2 = perm(jcol1)
        ilig2 = perm(ilig1)

        if (ilig2 .gt. jcol2) then
            perm_diag = .true._1
            jcol2 = ilig2
            ilig2 = perm(jcol1)
        else
            perm_diag = .false._1
        end if

        if (jcol2 .eq. 1) then
            k2_1 = 0
        else
            k2_1 = smdi2(jcol2-1)
        end if
        k2_2 = smdi2(jcol2)
        kterm2 = 0
        do k2 = k2_1+1, k2_2
            if (smhc2(k2) .eq. ilig2) then
                kterm2 = k2
                exit
            end if
        end do
        ASSERT(kterm2 .ge. 1 .and. kterm2 .le. nnz)

        if (perm_diag) then
            if (.not. lsym) then
                ASSERT(valm2_1(kterm2) .eq. 0.d0)
                valm2_1(kterm2) = valm1_2(kterm1)
                ASSERT(valm2_2(kterm2) .eq. 0.d0)
                valm2_2(kterm2) = valm1_1(kterm1)
            else
                ASSERT(valm2_1(kterm2) .eq. 0.d0)
                valm2_1(kterm2) = valm1_1(kterm1)
            end if
        else
            ASSERT(valm2_1(kterm2) .eq. 0.d0)
            valm2_1(kterm2) = valm1_1(kterm1)
            if (.not. lsym) then
                ASSERT(valm2_2(kterm2) .eq. 0.d0)
                valm2_2(kterm2) = valm1_2(kterm1)
            end if
        end if
    end do

!   2.2 Il faut modifier les objets .CCID et .CONL pour tenir compte de perm :
!   --------------------------------------------------------------------------
    if (refa2(3) .ne. ' ') then
        call jeveuo(matas2//'.CCID', 'E', vi=ccid)
        ASSERT(size(ccid) .eq. neq+1)
        AS_ALLOCATE(vi=cci2, size=neq+1)
        cci2(:) = ccid(:)
        do k = 1, neq
            ccid(perm(k)) = cci2(k)
        end do
        ccid(neq+1) = cci2(neq+1)
        AS_DEALLOCATE(vi=cci2)
    end if

    call jeexin(matas2//'.CONL', iexi)
    if (iexi .gt. 0) then
        call jeveuo(matas2//'.CONL', 'E', vr=conl)
        ASSERT(size(conl) .eq. neq)
        AS_ALLOCATE(vr=con2, size=neq)
        con2(:) = conl(:)
        do k = 1, neq
            conl(perm(k)) = con2(k)
        end do
        AS_DEALLOCATE(vr=con2)
    end if

    call jedema()

end subroutine
