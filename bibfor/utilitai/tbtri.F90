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

subroutine tbtri(ndim, tabint, tabchi, tabchr, tabchk)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!     FONCTION:
!     RANGEMENT DES ENTIERS DU TABLEAU TABCHAI
!     RANGEMENT DES REELS DU TABLEAU TABCHAR
!     RANGEMENT DES CHAINES DE CARACTERES DU TABLEAU TABCHAK
!     DANS L'ORDRE CROISSANT.
!-----------------------------------------------------------------------
! IN  NDIM   : I  : DIMENSION DU TABLEAU TABCHA.
! IN  TABCHI : I  : TABLEAU CONTENANT DES ENTIERS A RANGER
!                   DANS L'ORDRE CROISSANT.
! IN  TABCHR : R  : TABLEAU CONTENANT DES REELS A RANGER
!                   DANS L'ORDRE CROISSANT.
! IN  TABCHK : K  : TABLEAU CONTENANT DES CHAINES DE CARACTERES A RANGER
!                   DANS L'ORDRE CROISSANT.
! OUT TABINT : I  : TABLEAU D'ENTIERS CONTENANT LES POSITIONS
!                   DANS LE TABLEAU  TABCHA DANS L'ORDRE CROISSANT.
!-----------------------------------------------------------------------
!
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in), optional, target :: tabchi(*)
    real(kind=8), intent(in), optional, target :: tabchr(*)
    character(len=*), intent(in), optional, target :: tabchk(*)
    integer(kind=8), intent(out), optional, target :: tabint(*)
!
! ----------------------------------------------------------------------
    integer(kind=8) :: imin, j0, j1, i, j
    integer(kind=8), pointer :: masq(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    ASSERT(AU_MOINS_UN3(tabchi, tabchr, tabchk))
!
!     --- ON DEMASQUE TOUS LES ELEMENTS DU TABLEAU A TRIER ---
!
    AS_ALLOCATE(vi=masq, size=ndim)
!
    j0 = 1
    do i = 1, ndim
!        --- RECHERCHE DU PREMIER ELEMENT NON MASQUE ---
        do j = j0, ndim
            if (masq(j) .eq. 0) then
                j1 = j
                goto 22
            end if
        end do
!
22      continue
!
!        -- RECHERCHE DU PLUS PETIT ELEMENT NON MASQUE --
        j0 = j1
        imin = j1
        if (present(tabchi)) then
            do j = j0+1, ndim
                if (masq(j) .eq. 0 .and. tabchi(j) .lt. tabchi(imin)) imin = j
            end do
        else if (present(tabchr)) then
            do j = j0+1, ndim
                if (masq(j) .eq. 0 .and. tabchr(j) .lt. tabchr(imin)) imin = j
            end do
        else if (present(tabchk)) then
            do j = j0+1, ndim
                if (masq(j) .eq. 0 .and. tabchk(j) .lt. tabchk(imin)) imin = j
            end do
        end if
!
!        -- RANGEMENT DU IEME ELEMENT ET MISE A JOUR DU MASQUE --
        tabint(i) = imin
        masq(imin) = 1
!
    end do
!
    AS_DEALLOCATE(vi=masq)
!
    call jedema()
!
end subroutine
