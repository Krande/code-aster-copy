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

subroutine rldlg3(metres, lmat, xsol, cxsol, nbsol)
    implicit none
#include "jeveux.h"

#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/rldlc8.h"
#include "asterfort/rldlr8.h"
#include "asterfort/rlduc8.h"
#include "asterfort/rldur8.h"
#include "asterfort/rlfc16.h"
#include "asterfort/rltfr8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"

    character(len=*) :: metres
    integer(kind=8) :: lmat, nbsol
    real(kind=8) :: xsol(*)
    complex(kind=8) :: cxsol(*)

    integer(kind=8) :: typvar, typsym, iexi
    character(len=14) :: nu
    character(len=19) :: mat19, stolci
    integer(kind=8), pointer :: m2lc(:) => null()
    real(kind=8), pointer :: tempor(:) => null()
    complex(kind=8), pointer :: tempoc(:) => null()
    integer(kind=8) :: jscbl, jscdi, jschc, nbbloc, neq, isol, iequa
!-----------------------------------------------------------------------
    call jemarq()
    neq = zi(lmat+2)
    typvar = zi(lmat+3)
    typsym = zi(lmat+4)
    mat19 = zk24(zi(lmat+1)) (1:19)

    if (metres .eq. 'LDLT') then
!     ---------------------------------
        call jelira(mat19//'.UALF', 'NMAXOC', nbbloc)
        call mtdsc2(mat19, 'SCHC', 'L', jschc)
        call mtdsc2(mat19, 'SCDI', 'L', jscdi)
        call mtdsc2(mat19, 'SCBL', 'L', jscbl)

!       -- Permutation de xsol pour tenir compte de .M2LC :
!       ----------------------------------------------------
        call dismoi('NOM_NUME_DDL', mat19, 'MATR_ASSE', repk=nu)
        stolci = nu//'.SLCS'
        call jeexin(stolci//'.M2LC', iexi)
        if (iexi .gt. 0) then
            call jeveuo(stolci//'.M2LC', 'L', vi=m2lc)
            if (typvar .eq. 1) then
                AS_ALLOCATE(vr=tempor, size=neq)
                do isol = 1, nbsol
                    do iequa = 1, neq
                        tempor(m2lc(iequa)) = xsol(iequa+(isol-1)*neq)
                    end do
                    do iequa = 1, neq
                        xsol(iequa+(isol-1)*neq) = tempor(iequa)
                    end do
                end do
            else
                AS_ALLOCATE(vc=tempoc, size=neq)
                do isol = 1, nbsol
                    do iequa = 1, neq
                        tempoc(m2lc(iequa)) = cxsol(iequa+(isol-1)*neq)
                    end do
                    do iequa = 1, neq
                        cxsol(iequa+(isol-1)*neq) = tempoc(iequa)
                    end do
                end do
            end if
        end if

        if (typvar .eq. 1) then
!           --- SYSTEME REELLE ---

!         -- CAS D'UNE MATRICE SYMETRIQUE
            if (typsym .eq. 1) then
                call rldlr8(mat19, zi(jschc), zi(jscdi), zi(jscbl), neq, &
                            nbbloc, xsol, nbsol)

!        -- CAS D'UNE MATRICE NON_SYMETRIQUE
            else if (typsym .eq. 0) then
                call rldur8(mat19, zi(jschc), zi(jscdi), zi(jscbl), neq, &
                            nbbloc/2, xsol, nbsol)
            end if

        else if (typvar .eq. 2) then
!         -- SYSTEME COMPLEXE (SYMETRIQUE) ---
            if (typsym .eq. 1) then
                call rldlc8(mat19, zi(jschc), zi(jscdi), zi(jscbl), neq, &
                            nbbloc, cxsol, nbsol)
!        -- CAS D'UNE MATRICE NON_SYMETRIQUE
            else if (typsym .eq. 0) then
                call rlduc8(mat19, zi(jschc), zi(jscdi), zi(jscbl), neq, &
                            nbbloc/2, cxsol, nbsol)
            end if
        end if

!       -- permutation de xsol pour tenir compte de .M2LC :
!       ----------------------------------------------------
        if (iexi .gt. 0) then
            if (typvar .eq. 1) then
                do isol = 1, nbsol
                    do iequa = 1, neq
                        tempor(iequa) = xsol(m2lc(iequa)+(isol-1)*neq)
                    end do
                    do iequa = 1, neq
                        xsol(iequa+(isol-1)*neq) = tempor(iequa)
                    end do
                end do
                AS_DEALLOCATE(vr=tempor)
            else
                do isol = 1, nbsol
                    do iequa = 1, neq
                        tempoc(iequa) = cxsol(m2lc(iequa)+(isol-1)*neq)
                    end do
                    do iequa = 1, neq
                        cxsol(iequa+(isol-1)*neq) = tempoc(iequa)
                    end do
                end do
                AS_DEALLOCATE(vc=tempoc)
            end if
        end if

    else if (metres .eq. 'MULT_FRONT') then
!     ------------------------------------
        if (typvar .eq. 1) then
            call rltfr8(mat19, neq, xsol, nbsol, typsym)
        else if (typvar .eq. 2) then
            call rlfc16(mat19, neq, cxsol, nbsol, typsym)
        end if
    else
        ASSERT(.false.)
    end if

    call jedema()

end subroutine
