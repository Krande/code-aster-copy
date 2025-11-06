! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
!
interface
    subroutine cafelsqp_gap(cequi, effm, ht, enrobi, enrobs, kt, eys, &
                            phiinf, phisup, dnsinf, dnssup, sigmsi, sigmss, &
                            sigmci, sigmcs, alpha, etat, unite_m, fctm, &
                            wfins, wfini)
        real(kind=8), intent(in) :: cequi
        real(kind=8), intent(in) :: effm
        real(kind=8), intent(in) :: ht
        real(kind=8), intent(in) :: enrobi
        real(kind=8), intent(in) :: enrobs
        real(kind=8), intent(in) :: kt
        real(kind=8), intent(in) :: eys
        real(kind=8), intent(in) :: phiinf
        real(kind=8), intent(in) :: phisup
        real(kind=8), intent(in) :: dnsinf
        real(kind=8), intent(in) :: dnssup
        real(kind=8), intent(in) :: sigmsi
        real(kind=8), intent(in) :: sigmss
        real(kind=8), intent(in) :: sigmci
        real(kind=8), intent(in) :: sigmcs
        real(kind=8), intent(in) :: alpha
        integer(kind=8), intent(in) :: etat
        real(kind=8), intent(in) :: unite_m, fctm
        real(kind=8), intent(out) :: wfins
        real(kind=8), intent(out) :: wfini
    end subroutine cafelsqp_gap
end interface
