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
! aslint: disable=W1504
!
subroutine lcjohm(imate, lSigm, lMatr, lVari, kpi, npg, &
                  nomail, addeme, advico, ndim, dimdef, &
                  dimcon, nbvari, defgem, defgep, varim, &
                  varip, sigm, sigp, drde, ouvh, &
                  retcom)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: imate, kpi, npg, addeme, advico, ndim, dimdef, dimcon, nbvari
    real(kind=8) :: defgem(dimdef), varim(nbvari), sigm(dimcon)
    character(len=8) :: nomail
    aster_logical, intent(in) :: lSigm, lMatr, lVari
    integer(kind=8) :: retcom
    real(kind=8) :: defgep(dimdef), varip(nbvari), sigp(dimcon)
    real(kind=8) :: drde(dimdef, dimdef), ouvh
!
! - VARIABLES LOCALES
!
    integer(kind=8) :: i, kpg, spt
    real(kind=8) :: kni, umc, gamma, kt, clo, para(4), valr(2), tmecn, tmecs
    character(len=8) :: fami, poum
    integer(kind=8) :: icodre(18)
    character(len=8), parameter :: ncra1(4) = (/'K    ', 'DMAX ', 'GAMMA', 'KT   '/)
!
! - RECUPERATION DES PARAMETRES MATERIAU
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'JOINT_BANDIS', 0, ' ', [0.d0], &
                4, ncra1, para, icodre, 1)
    kni = para(1)
    umc = para(2)
    gamma = para(3)
    kt = para(4)
!
! - MISE A JOUR FERMETURE
    clo = 0.d0
    ouvh = varim(advico)
    clo = umc-ouvh
    clo = clo-defgep(addeme)+defgem(addeme)
!
! - Internal state variable
!
    if (lVari) then
        ouvh = umc-clo
        varip(advico) = ouvh
    end if
!
! - Stresses
!
    if (lSigm) then
        if ((clo .gt. umc) .or. (clo .lt. -1.d-3)) then
            valr(1) = clo
            valr(2) = umc
            call utmess('A', 'ALGORITH17_11', sk=nomail, nr=2, valr=valr)
            retcom = 1
            goto 999
        end if
        do i = 1, dimcon
            sigp(i) = 0.d0
        end do
        sigp(1) = sigm(1)-kni/(1-clo/umc)**gamma*(varim(advico)-varip(advico))
        do i = 2, ndim
            sigp(i) = sigm(i)+kt*(defgep(addeme+1)-defgem(addeme+1))
        end do
    end if

!
! - CALCUP OPERATEUR TANGENT
!
    if (lMatr .and. (kpi .le. npg)) then
        tmecn = kni/(1-clo/umc)**gamma
        tmecs = kt
        drde(addeme, addeme) = tmecn
        do i = 2, ndim
            drde(addeme+i-1, addeme+i-1) = tmecs
        end do
    end if
!
999 continue
!
end subroutine
