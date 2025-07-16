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

subroutine matr_asse_scale(matasz, lvect, rvect)
! person_in_charge: nicolas.tardieu at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/isParallelMatrix.h"

    character(len=*), intent(in) :: matasz
    real(kind=8), intent(in) :: lvect(*), rvect(*)
!-----------------------------------------------------------------------
! Goal : scale the assembly matrix by right and left multiplication by
!        diagonal matrices stored as vectors.
!
!  matas : name of the assembly matrix
!  lvect : left array
!  rvect : right array
!
!-----------------------------------------------------------------------

    aster_logical :: l_parallel_matrix, lmd, lsym
    character(len=1) :: ktyp
    character(len=14) :: nonu
    character(len=19) :: matass, kstoc
    integer(kind=8) :: n, n1, nz, kterm, jnlogl, nsmdi, jsmhc, nsmhc, jcolg
    integer(kind=8) :: jprddl, nvale, jvale, nlong, coltmp, ier, iligg, jcoll, iligl, jval2
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    !-------------------------------------------------------------------
    call jemarq()
!
    matass = matasz
!
    call jeveuo(matass//'.REFA', 'E', vk24=refa)
    nonu = refa(2) (1:14)
!
    ! force to recompute the blocks to account for BCs
    if (refa(3) .eq. 'ELIMF') then
        refa(3) = 'ELIML'
        call jedetr(matass//'.CCLL')
        call jedetr(matass//'.CCVA')
        call jedetr(matass//'.CCII')
    end if
!
    ! storage of the matrix
    kstoc = nonu//'.SMOS'
    call jeexin(kstoc//'.SMDI', ier)
    if (ier .eq. 0) then
        call utmess('F', 'ALGELINE3_8', sk=matass)
    end if
!
    lmd = (refa(11) .eq. 'MATR_DISTR')
!
    l_parallel_matrix = isParallelMatrix(matass)
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    if (lmd) then
        call jelira(nonu//'.NUML.DELG', 'LONMAX', n1)
        call jeveuo(nonu//'.NUML.NULG', 'L', jnlogl)
    else
        call jelira(nonu//'.NUME.DELG', 'LONMAX', n1)
        jnlogl = 0
        if (l_parallel_matrix) then
            call jeveuo(nonu//'.NUME.PDDL', 'L', jprddl)
        end if
    end if
    ASSERT(n1 .eq. nsmdi)
!     --- CALCUL DE N
    n = nsmdi
!     --- CALCUL DE NZ
    nz = smdi(n)
!
    ASSERT(nz .le. nsmhc)
    call jelira(matass//'.VALM', 'NMAXOC', nvale)
    if (nvale .eq. 1) then
        lsym = .true.
    else if (nvale .eq. 2) then
        lsym = .false.
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(jexnum(matass//'.VALM', 1), 'L', jvale)
    call jelira(jexnum(matass//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
    if (.not. lsym) then
        call jeveuo(jexnum(matass//'.VALM', 2), 'L', jval2)
        call jelira(jexnum(matass//'.VALM', 2), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
!
    call jelira(jexnum(matass//'.VALM', 1), 'TYPE', cval=ktyp)
    if (ktyp .ne. 'R') call utmess('F', 'ALGELINE_9')

    ! scale the matrix
    jcoll = 1
    do kterm = 1, nz
!
!       --- PARTIE TRIANGULAIRE SUPERIEURE
        if (smdi(jcoll) .lt. kterm) jcoll = jcoll+1
        iligl = zi4(jsmhc-1+kterm)
        if (lmd) then
            iligg = zi(jnlogl+iligl-1)
            jcolg = zi(jnlogl+jcoll-1)
        else
            iligg = iligl
            jcolg = jcoll
        end if
        if ((.not. lsym) .and. (iligg .ge. jcolg)) then
            coltmp = jcolg
            jcolg = iligg
            iligg = coltmp
        end if
        zr(jvale-1+kterm) = zr(jvale-1+kterm)*rvect(jcoll)*lvect(iligl)
!
!        --- PARTIE TRIANGULAIRE INFERIEURE
        if ((.not. lsym) .and. (iligg .ne. jcolg)) then
            zr(jval2-1+kterm) = zr(jval2-1+kterm)*rvect(iligl)*lvect(jcoll)
        end if
!
    end do
!
!
    call jedema()
end subroutine
