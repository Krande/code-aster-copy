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

subroutine extdia(matr, numddl, icode, diag)
    implicit none
! 15/03/91    G.JACQUART AMV/P61 47 65 49 41
!***********************************************************************
!
!     FONCTION : EXTRACTION DE LA DIAGONALE D'UNE MATRICE
!
!-----------------------------------------------------------------------
!    MATR   /I/ : NOM DE LA MATRICE
!    NUMDDL /I/ : NUMEROTATION ASSOCIEE A MATR
!    ICODE  /I/ : 2 SI CALCUL TRANSITOIRE DIRECT
!                 1 SI SOUS-STRUCTURATION DYNAMIQUE TRANSITOIRE
!                   SANS DOUBLE PROJECTION
!                 0 SINON
!    DIAG   /O/ : VECTEUR CONTENANT LA DIAGONALE DE MATR
!-----------------------------------------------------------------------
!   Note : if the matrix is complex, return the moduli of the diagonal terms
!
!
!
!
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtdscr.h"
#include "asterfort/typddl.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=24) :: numddl
    character(len=8) :: matr
    real(kind=8) :: diag(*)
    integer(kind=8) :: icode
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    aster_logical :: iscmplx
    integer(kind=8) :: idia, j, jbloc, k
    integer(kind=8) :: l, lmat, nbacti, nbbloq, nblagr, nbliai, neq
    character(len=24) :: typmatr
    integer(kind=8), pointer :: vtypddl(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: smde(:) => null()

!
!-----------------------------------------------------------------------
    call jemarq()
    call mtdscr(matr)
    call jeveuo(matr//'           .&INT', 'L', lmat)
!
    call gettco(matr, typmatr)
    iscmplx = (typmatr(1:9) .eq. 'MATR_ASSE') .and. (typmatr(16:16) .eq. 'C')
!
    call jeveuo(numddl(1:8)//'      .SMOS.SMDE', 'L', vi=smde)
    neq = smde(1)
!
    AS_ALLOCATE(vi=vtypddl, size=neq)
    call typddl('ACTI', numddl(1:8), neq, vtypddl, nbacti, &
                nbbloq, nblagr, nbliai)
    if (icode .eq. 2) then
        if (nbliai .gt. 0) then
            call utmess('F', 'UTILITAI_76')
        end if
    end if
!
    call jeveuo(numddl(1:8)//'      .SMOS.SMDI', 'L', vi=smdi)
    k = 0
    l = 0
    call jeveuo(jexnum(matr//'           .VALM', 1), 'L', jbloc)
    if (.not. (iscmplx)) then
        do j = 1, neq
            k = k+1
            if (vtypddl(k) .ne. 0) then
                idia = smdi(k)
                l = l+1
                diag(l) = zr(jbloc-1+idia)
            else if (icode .eq. 0 .or. icode .eq. 2) then
                l = l+1
                diag(l) = 0.d0
            end if
        end do
    else
        do j = 1, neq
            k = k+1
            if (vtypddl(k) .ne. 0) then
                idia = smdi(k)
                l = l+1
                diag(l) = sqrt(real(zc(jbloc-1+idia))**2+imag(zc(jbloc-1+idia))**2)
            else if (icode .eq. 0 .or. icode .eq. 2) then
                l = l+1
                diag(l) = 0.d0
            end if
        end do
    end if
    AS_DEALLOCATE(vi=vtypddl)
!
    call jedema()
end subroutine
