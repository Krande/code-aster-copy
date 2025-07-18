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
subroutine acedat(motfac, in, npara, sec, exp, &
                  tab, car)
    implicit none
    integer(kind=8) :: in, npara(*)
    character(len=*) :: motfac, sec(*), exp(*), tab(*), car(*)
!     AFFE_CARA_ELEM
!     INITIALISATION DES PARAMETRES ET DES DATAS
! ----------------------------------------------------------------------
!
!     --- MOT CLE FACTEUR "POUTRE"-------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, k, nbobar, nbopou, ncerba, ncerpo, ngenba
    integer(kind=8) :: ngenpo, nrecba, nrecpo, nsecba, nsecpo, ntseba, ntsepo
!
!-----------------------------------------------------------------------
    parameter(ngenpo=15, nrecpo=6, ncerpo=2)
    parameter(nsecpo=2, ntsepo=3, nbopou=44)
    character(len=8) :: exppou(nbopou), tabpou(nbopou)
    character(len=8) :: carpou(ngenpo*(nsecpo+1), ntsepo)
    character(len=16) :: secpou(ntsepo)
!
!     --- MOT CLE FACTEUR "BARRE" -------------------------------------
    parameter(ngenba=1, nrecba=6, ncerba=2)
    parameter(nsecba=1, ntseba=3, nbobar=8)
    character(len=8) :: expbar(nbobar), tabbar(nbobar), carbar(nrecba, ntseba)
    character(len=16) :: secbar(ntseba)
!
!     --- POUTRE -------------------------------------------------------
    data tabpou/'A1      ', 'IY1     ', 'IZ1     ', 'AY1     ',&
     &             'AZ1     ', 'EY1     ', 'EZ1     ', 'JX1     ',&
     &             'RY1     ', 'RZ1     ', 'RT1     ',&
     &             'A2      ', 'IY2     ', 'IZ2     ', 'AY2     ',&
     &             'AZ2     ', 'EY2     ', 'EZ2     ', 'JX2     ',&
     &             'RY2     ', 'RZ2     ', 'RT2     ', 'TVAR    ',&
     &             'HY1     ', 'HZ1     ', 'EPY1    ', 'EPZ1    ',&
     &             'HY2     ', 'HZ2     ', 'EPY2    ', 'EPZ2    ',&
     &             'R1      ', 'EP1     ', 'R2      ', 'EP2     ',&
     &             'TSEC    ', 'AI1     ', 'AI2     ', 'JG1     ',&
     &             'JG2     ', 'IYR21   ', 'IYR22   ', 'IZR21   ',&
     &             'IZR22   '/
    data exppou/'A       ', 'IY      ', 'IZ      ', 'AY      ',&
     &             'AZ      ', 'EY      ', 'EZ      ', 'JX      ',&
     &             'RY      ', 'RZ      ', 'RT      ',&
     &             'A       ', 'IY      ', 'IZ      ', 'AY      ',&
     &             'AZ      ', 'EY      ', 'EZ      ', 'JX      ',&
     &             'RY      ', 'RZ      ', 'RT      ', 'TVAR    ',&
     &             'HY      ', 'HZ      ', 'EPY     ', 'EPZ     ',&
     &             'HY      ', 'HZ      ', 'EPY     ', 'EPZ     ',&
     &             'R       ', 'EP      ', 'R       ', 'EP      ',&
     &             'TSEC    ', 'AI      ', 'AI      ', 'JG      ',&
     &             'JG      ', 'IYR2    ', 'IYR2    ', 'IZR2    ',&
     &             'IZR2    '/
    data secpou(1)/'GENERALE  '/
    data(carpou(i, 1), i=1, ngenpo*(nsecpo+1))/&
     &             'A       ', 'IY      ', 'IZ      ', 'AY      ',&
     &             'AZ      ', 'EY      ', 'EZ      ', 'JX      ',&
     &             'RY      ', 'RZ      ', 'RT      ', 'AI      ',&
     &             'JG      ', 'IYR2    ', 'IZR2    ',&
     &             'A1      ', 'IY1     ', 'IZ1     ', 'AY1     ',&
     &             'AZ1     ', 'EY1     ', 'EZ1     ', 'JX1     ',&
     &             'RY1     ', 'RZ1     ', 'RT1     ', 'AI1     ',&
     &             'JG1     ', 'IYR21   ', 'IZR21   ',&
     &             'A2      ', 'IY2     ', 'IZ2     ', 'AY2     ',&
     &             'AZ2     ', 'EY2     ', 'EZ2     ', 'JX2     ',&
     &             'RY2     ', 'RZ2     ', 'RT2     ', 'AI2     ',&
     &             'JG2     ', 'IYR22   ', 'IZR22   '/
    data secpou(2)/'RECTANGLE'/
    data(carpou(i, 2), i=1, nrecpo*(nsecpo+1))/&
     &             'H       ', 'HY      ', 'HZ      ',&
     &             'EP      ', 'EPY     ', 'EPZ     ',&
     &             'H1      ', 'HY1     ', 'HZ1     ',&
     &             'EP1     ', 'EPY1    ', 'EPZ1    ',&
     &             'H2      ', 'HY2     ', 'HZ2     ',&
     &             'EP2     ', 'EPY2    ', 'EPZ2    '/
    data secpou(3)/'CERCLE    '/
    data(carpou(i, 3), i=1, ncerpo*(nsecpo+1))/&
     &             'R       ', 'EP      ',&
     &             'R1      ', 'EP1     ',&
     &             'R2      ', 'EP2     '/
!
!     --- BARRE --------------------------------------------------------
    data expbar/'A', 'HY', 'HZ', 'EPY', 'EPZ', 'R', 'EP', 'TSEC'/
    data tabbar/'A1', 'HY1', 'HZ1', 'EPY1', 'EPZ1', 'R1', 'EP1', 'TSEC'/
    data secbar(1)/'GENERALE'/
    data(carbar(i, 1), i=1, ngenba)/'A'/
    data secbar(2)/'RECTANGLE'/
    data(carbar(i, 2), i=1, nrecba)/'H', 'HY', 'HZ', 'EP', 'EPY', 'EPZ'/
    data secbar(3)/'CERCLE '/
    data(carbar(i, 3), i=1, ncerba)/'R', 'EP'/
!     ------------------------------------------------------------------
    if (motfac .eq. 'POUTRE') then
        npara(1) = nsecpo
        npara(2) = ntsepo
        npara(3) = nbopou
        npara(4) = nsecpo*ngenpo
        npara(5) = nsecpo*ngenpo
        npara(6) = ngenpo
        npara(7) = nrecpo
        npara(8) = ncerpo
        if (in .eq. 0) goto 999
        do i = 1, ntsepo
            sec(i) = secpou(i)
        end do
        do i = 1, nbopou
            exp(i) = exppou(i)
            tab(i) = tabpou(i)
        end do
        k = 0
        do i = 1, ngenpo*(nsecpo+1)
            k = k+1
            car(k) = carpou(i, 1)
        end do
        do i = 1, nrecpo*(nsecpo+1)
            k = k+1
            car(k) = carpou(i, 2)
        end do
        k = 2*ngenpo*(nsecpo+1)
        do i = 1, ncerpo*(nsecpo+1)
            k = k+1
            car(k) = carpou(i, 3)
        end do
    else if (motfac .eq. 'BARRE') then
        npara(1) = nsecba
        npara(2) = ntseba
        npara(3) = nbobar
        npara(4) = nsecba*nrecba
        npara(5) = nsecba*nrecba
        npara(6) = ngenba
        npara(7) = nrecba
        npara(8) = ncerba
        if (in .eq. 0) goto 999
        do i = 1, ntseba
            sec(i) = secbar(i)
        end do
        do i = 1, nbobar
            exp(i) = expbar(i)
            tab(i) = tabbar(i)
        end do
        k = 0
        do i = 1, ngenba
            k = k+1
            car(k) = carbar(i, 1)
        end do
        k = nrecba
        do i = 1, nrecba
            k = k+1
            car(k) = carbar(i, 2)
        end do
        k = 2*nrecba
        do i = 1, ncerba
            k = k+1
            car(k) = carbar(i, 3)
        end do
    end if
!
999 continue
end subroutine
