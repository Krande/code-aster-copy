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

subroutine dikpkt(imater, nomphe, kp, kt)
    implicit none
    integer(kind=8), intent(in)     :: imater
    character(len=16), intent(in)   :: nomphe
    real(kind=8), intent(out)       :: kp, kt
!
#include "asterfort/rcvala.h"
!
! --------------------------------------------------------------------------------------------------
! Récupération des raideurs associées à un élément élastique parallèle à un élément discret
! --------------------------------------------------------------------------------------------------
!
!   in:
!        imater : adresse du matériau codé
!        nomphe : nom du phénomène dans DEFI_MATERIAU (e.g. "DIS_CONTACT") pour lire (KP,KT)
!
!   out:
!        kp     : raideur selon l'axe x du repère local de l'élément
!        kt     : raideur selon les axes y et z du repère local de l'élément
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: donatien.rossat at edf.fr
!
    integer(kind=8), parameter :: nbres = 2
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=8) :: nomres(nbres)
!
    data nomres/'KP', 'KT'/
!
    valres(:) = 0.d0
    call rcvala(imater, ' ', nomphe, 0, ' ', &
                [0.0d0], nbres, nomres, valres, codres, &
                0, nan='NON')
    if (codres(1) .ne. 0) then
        kp = 0.d0
    else
        kp = valres(1)
    end if
    if (codres(2) .ne. 0) then
        kt = 0.d0
    else
        kt = valres(2)
    end if
!
!
end subroutine
