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

subroutine editgd(ncmp, nedit, dg, ncmpmx, ctype, &
                  jnocmp, jncmp, jvalv, jvale)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/jacopo.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: ncmp
    integer(kind=8), intent(in) :: nedit
    integer(kind=8), intent(inout) :: dg(*)
    integer(kind=8), intent(in) :: ncmpmx
    character(len=8), intent(in) :: ctype
    integer(kind=8), intent(in) :: jnocmp
    integer(kind=8), intent(in) :: jncmp
    integer(kind=8), intent(in) :: jvalv
    integer(kind=8), intent(in) :: jvale
! ----------------------------------------------------------------------
!     entrees:
!     ncmp  : nombre de cmp a stocker
!
!     sorties:
!     dg  : descripteur_grandeur a mettre a jour
!
! ----------------------------------------------------------------------
    integer(kind=8) :: ior
!
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: i, j, iec, reste, code, deb2, debgd
    integer(kind=8) :: lshift
    character(len=8) :: nomcmp
    character(len=24) :: valk
!-----------------------------------------------------------------------
    integer(kind=8) :: ico, indgd, ncmp2
!-----------------------------------------------------------------------
    debgd = (nedit-1)*ncmpmx
!
!   -- on compte le nombre de cmps a noter reellement :
    ncmp2 = 0
    do i = 1, ncmp
        if (zk8(jncmp-1+i) (1:1) .ne. ' ') ncmp2 = ncmp2+1
    end do
!
    indgd = 0
    ico = 0
    do i = 1, ncmpmx
        nomcmp = zk8(jnocmp-1+i)
        j = indik8(zk8(jncmp), nomcmp, 1, ncmp)
        if (j .ne. 0) then
            ico = ico+1
            indgd = indgd+1
            iec = (i-1)/30+1
            reste = i-30*(iec-1)
            code = lshift(1, reste)
            dg(iec) = ior(dg(iec), code)
            deb2 = debgd+indgd
            call jacopo(1, ctype, jvalv+j-1, jvale+deb2-1)
        end if
    end do
    if (ico .ne. ncmp2) then
        call utmess('F+', 'CALCULEL6_68')
        do i = 1, ncmpmx
            nomcmp = zk8(jnocmp-1+i)
            valk = nomcmp
            call utmess('F+', 'CALCULEL6_69', sk=valk)
        end do
        call utmess('F+', 'CALCULEL6_70')
        do i = 1, ncmp
            nomcmp = zk8(jncmp-1+i)
            valk = nomcmp
            call utmess('F+', 'CALCULEL6_71', sk=valk)
        end do
        call utmess('F', 'VIDE_1')
    end if
!
end subroutine
