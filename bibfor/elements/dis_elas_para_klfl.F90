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
!
subroutine dis_elas_para_klfl(for_discret, kp, kt1, kt2, utotxyz, klvp, flp)
!
    use te0047_type
    implicit none
#include "asterf_types.h"
#include "asterfort/diklvraid.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    real(kind=8), intent(in)      :: kp, kt1, kt2
    real(kind=8), intent(in)      :: utotxyz(3)
    real(kind=8), allocatable, intent(out) :: klvp(:)
    real(kind=8), intent(out)     :: flp(3)
!
! --------------------------------------------------------------------------------------------------
!
! Calcul de la matrice tangente et des forces d'un élément élastique parallèle à un élément discret
!
! --------------------------------------------------------------------------------------------------
! in  :
!       for_discret : type dérivé de l'élément discret
!       kp          : raideur selon l'axe x du repère local du discret élastique en parallèle
!       kt1         : raideur selon l'axe y du repère local du discret élastique en parallèle
!       kt2         : raideur selon l'axe z du repère local du discret élastique en parallèle
!       utotxyz     : déplacement total dans le repère local du discret élastique en parallèle
!
! out :
!       klvp        : matrice tangente dans le repère local du discret élastique en parallèle
!       flp         : vecteur force dans le repère local du discret élastique en parallèle
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: donatien.rossat at edf.fr
!
    integer(kind=8) :: nc, nn, neq, nsym
    real(kind=8) :: raide(6)
! --------------------------------------------------------------------------------------------------

!   Matrice de raideur de l'élément élastique dans le repère local
    nn = for_discret%nno
    nc = for_discret%nc
    neq = nc*nn
    nsym = neq*(neq+1)/2
    raide(1:6) = 0.d0
    raide(1) = kp
    raide(2) = kt1
    if (for_discret%ndim .eq. 3) then
        raide(3) = kt2
    end if
    allocate (klvp(nsym))
    klvp(1:nsym) = 0.d0
    call diklvraid(for_discret%nomte, klvp, raide)

!   Vecteur force de l'élément élastique dans le repère local
    flp(1:3) = 0.d0
    flp(1) = kp*utotxyz(1)
    flp(2) = kt1*utotxyz(2)
    if (for_discret%ndim .eq. 3) then
        flp(3) = kt2*utotxyz(3)
    end if

end subroutine
