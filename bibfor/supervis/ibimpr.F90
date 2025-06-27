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

subroutine ibimpr()
    implicit none
!     DEFINITION DES UNITES LOGIQUES DES IMPRESSIONS
!     ------------------------------------------------------------------
!     Conserver la cohérence avec iunifi et ulopen.
!     L'unité 15 est ouverte dans debut.
!
#include "asterfort/uldefi.h"
!
    integer(kind=8) :: mximpr
    parameter(mximpr=3)
    character(len=16) :: nompr(mximpr)
    integer(kind=8) :: unitpr(mximpr)
    character(len=1) :: autpr(mximpr)
    integer(kind=8) :: passe
    save passe
    data passe/0/
    data nompr/'MESSAGE', 'RESULTAT', 'ERREUR'/
    data unitpr/6, 8, 9/
    data autpr/'N', 'O', 'N'/
    integer(kind=8) :: i
!   ------------------------------------------------------------------
    passe = passe+1
!
! --- DEFINITION DES UNITES STANDARDS
    if (passe .eq. 1) then
        do i = 1, mximpr
            call uldefi(unitpr(i), ' ', nompr(i), 'A', 'N', &
                        autpr(i))
        end do
    end if
!
end subroutine
