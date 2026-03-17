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

subroutine cafelsqp_gap(cequi, effm, ht, enrobi, enrobs, kt, eys, &
                        phiinf, phisup, dnsinf, dnssup, sigmsi, sigmss, &
                        sigmci, sigmcs, alpha, etat, unite_m, fctm, &
                        wfins, wfini)

    implicit none

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
    real(kind=8), intent(out)  :: wfins
    real(kind=8), intent(out)  :: wfini
!
    real(kind=8) :: ScTracMAX, ScTracMIN, xAN, f1, f2, f3, f4
    real(kind=8) :: hceff, rhoeff, SrmaxSUP, SrmaxINF, DESUP, DEINF

    f1 = 0.8
    f4 = 0.425

    if ((alpha .ge. 0) .AND. (alpha .le. 1)) then
        if (effm .gt. 0) then
            xAN = alpha*(ht-enrobi)
        elseif (effm .lt. 0) then
            xAN = alpha*(ht-enrobs)
        end if
    else
        xAN = -1
    end if
    ScTracMAX = min(sigmci, sigmcs)
    ScTracMIN = max(sigmci, sigmcs)
    if (ScTracMIN .gt. 0) then
        ScTracMIN = 0
    end if
    if (ScTracMAX .ge. 0) then
        f2 = 0.5
    else
        f2 = (ScTracMIN+ScTracMAX)/(2*ScTracMAX)
    end if
    if ((sigmss .ge. 0) .OR. (dnssup .eq. 0) .OR. (etat .eq. 1)) then
        wfins = 0
    else
        if ((enrobs*unite_m) .lt. 25) then
            f3 = 3.4
        else
            f3 = 3.4*((25.0/(enrobs*unite_m))**(2.0/3.0))
        end if
        if (xAN .ne. -1) then
            hceff = min(2.5*enrobs, 0.5*ht, (ht-xAN)/3.0)
        else
            hceff = min(2.5*enrobs, 0.5*ht)
        end if
        rhoeff = dnssup/hceff
        SrmaxSUP = f3*enrobs+f1*f2*f4*(phisup/rhoeff)
        if (xAN .ne. -1) then
            SrmaxSUP = min(SrmaxSUP, 1.3*(ht-xAN))
        else
            SrmaxSUP = min(SrmaxSUP, 1.3*ht)
        end if
        DESUP = (-sigmss-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
        DESUP = max(DESUP, -0.6*sigmss/eys)
        wfins = SrmaxSUP*DESUP
    end if
    if ((sigmsi .ge. 0) .OR. (dnsinf .eq. 0) .OR. (etat .eq. 1)) then
        wfini = 0
    else
        if ((enrobi*unite_m) .lt. 25) then
            f3 = 3.4
        else
            f3 = 3.4*((25.0/(enrobi*unite_m))**(2.0/3.0))
        end if
        if (xAN .ne. -1) then
            hceff = min(2.5*enrobi, 0.5*ht, (ht-xAN)/3.0)
        else
            hceff = min(2.5*enrobi, 0.5*ht)
        end if
        rhoeff = dnsinf/hceff
        SrmaxINF = f3*enrobi+f1*f2*f4*(phiinf/rhoeff)
        if (xAN .ne. -1) then
            SrmaxINF = min(SrmaxINF, 1.3*(ht-xAN))
        else
            SrmaxINF = min(SrmaxINF, 1.3*ht)
        end if
        DEINF = (-sigmsi-kt*(fctm/rhoeff)*(1.0+cequi*rhoeff))/eys
        DEINF = max(DEINF, -0.6*sigmsi/eys)
        wfini = SrmaxINF*DEINF
    end if

end subroutine
