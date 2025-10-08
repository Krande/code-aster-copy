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

subroutine dikpkt(imater, nomphe, kp, kt1, kt2)
    implicit none
    integer(kind=8), intent(in)     :: imater
    character(len=*), intent(in)   :: nomphe
    real(kind=8), intent(out)       :: kp, kt1, kt2
!
#include "asterc/r8prem.h"
#include "asterfort/rcvala.h"
#include "asterfort/utmess.h"
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
!        kt1    : raideur selon l'axe y du repère local de l'élément
!        kt2    : raideur selon l'axe z du repère local de l'élément
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: donatien.rossat at edf.fr
!
    integer(kind=8), parameter :: nbres = 4
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    real(kind=8) :: kt
    character(len=8) :: nomres(nbres)
    character(len=32) :: messak(1)
!
    data nomres/'KP', 'KT', 'KT1', 'KT2'/
!
    ! Lecture des paramètres de raideur de l'élément élastique en parallèle
    valres(:) = 0.d0
    call rcvala(imater, ' ', nomphe, 0, ' ', &
                [0.0d0], nbres, nomres, valres, codres, &
                0, nan='NON')
    ! Raideur normale
    if (codres(1) .ne. 0) then
        kp = 0.d0
    else
        kp = valres(1)
    end if
    ! Raideurs tangentielles
    if (codres(2) .ne. 0) then
        kt = 0.d0
    else
        kt = valres(2)
    end if
    if (codres(3) .ne. 0) then
        kt1 = 0.d0
    else
        kt1 = valres(3)
    end if
    if (codres(4) .ne. 0) then
        kt2 = 0.d0
    else
        kt2 = valres(4)
    end if
    if ((abs(kt) .gt. r8prem()) .and. ((abs(kt1) .gt. r8prem()) &
                                       .or. (abs(kt2) .gt. r8prem()))) then
        messak(1) = nomphe
        call utmess('F', 'DISCRETS_38', nk=1, valk=messak)
    else if (abs(kt) .gt. r8prem()) then
        kt1 = kt
        kt2 = kt
    end if
!
!
end subroutine
