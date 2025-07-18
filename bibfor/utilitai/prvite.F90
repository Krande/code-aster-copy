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
subroutine prvite(vec1, long, ip1, ip2, itp)
    implicit none
!      CALCUL DU PROFIL DE VITESSE
! ----------------------------------------------------------------------
!
#include "jeveux.h"
    integer(kind=8) :: ip(3), long, ip1, ip2, itp
    real(kind=8) :: angle(71), vite(71), angl, vec1(long)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ij, j, k
    integer(kind=8) :: kk, long2
    real(kind=8) :: alfa, beta
!-----------------------------------------------------------------------
    data(angle(i), i=1, 17)/&
     &   0.d0, 15.3d0, 30.d0, 40.d0, 48.5d0, 61.3d0,&
     &  75.d0, 83.6d0, 90.d0, 96.4d0, 105.d0, 118.7d0,&
     & 131.5d0, 140.d0, 150.d0, 164.7d0, 180.d0/
!
    data(angle(i), i=18, 42)/&
     &   0.d0, 14.5d0, 29.d0, 33.49d0, 33.51d0, 39.d0,&
     &  47.4d0, 60.d0, 70.99d0, 71.01d0, 75.4d0, 84.d0,&
     &  90.d0, 96.d0, 104.6d0, 108.99d0, 109.01d0, 120.d0,&
     & 132.6d0, 141.d0, 146.49d0, 146.51d0, 151.d0, 165.d0,&
     & 180.d0/
!
    data(angle(i), i=43, 71)/&
     &   0.d0, 14.8d0, 25.7d0, 25.72d0, 30.7d0, 39.6d0,&
     &  48.d0, 51.42d0, 51.44d0, 60.3d0, 75.2d0, 77.13d0,&
     &  77.15d0, 84.d0, 90.d0, 96.d0, 102.85d0, 102.87d0,&
     & 104.8d0, 119.7d0, 128.56d0, 128.58d0, 132.d0, 140.4d0,&
     & 149.3d0, 154.28d0, 154.3d0, 165.2d0, 180.d0/
!
    data(vite(j), j=1, 17)/&
     &   0.d0, 1.08d0, 1.77d0, 1.72d0, 1.33d0, 0.95d0,&
     &   0.67d0, 0.59d0, 0.58d0, 0.59d0, 0.67d0, 0.95d0,&
     &   1.33d0, 1.72d0, 1.77d0, 1.08d0, 0.d0/
!
    data(vite(j), j=18, 42)/&
     &   0.d0, 0.2d0, 0.7d0, 1.d0, 0.18d0, 0.32d0,&
     &   0.48d0, 0.61d0, 0.7d0, 0.16d0, 0.29d0, 0.5d0,&
     &   0.54d0, 0.5d0, 0.29d0, 0.16d0, 0.7d0, 0.61d0,&
     &   0.48d0, 0.32d0, 0.18d0, 1.d0, 0.7d0, 0.2d0,&
     &   0.d0/
!
    data(vite(j), j=43, 71)/&
     &   0.d0, 0.48d0, 0.69d0, 0.27d0, 0.54d0, 0.89d0,&
     &   1.16d0, 1.24d0, 0.02d0, 0.08d0, 0.79d0, 0.86d0,&
     &   0.17d0, 0.49d0, 0.55d0, 0.49d0, 0.17d0, 0.86d0,&
     &   0.79d0, 0.08d0, 0.02d0, 1.24d0, 1.16d0, 0.89d0,&
     &   0.54d0, 0.27d0, 0.69d0, 0.48d0, 0.d0/
!
    data(ip(ij), ij=1, 3)/1, 18, 43/
!
!
    long2 = long/2
!
    do kk = 1, long2
        if (kk .gt. ip1 .and. kk .lt. ip2) then
            angl = 180.d0*(vec1(kk)-vec1(ip1))/(vec1(ip2)-vec1(ip1))
            k = ip(itp)
20          continue
            if (angl .gt. angle(k+1)) then
                k = k+1
                goto 20
            end if
            alfa = (vite(k+1)-vite(k))/(angle(k+1)-angle(k))
            beta = (vite(k)*angle(k+1)-(vite(k+1)*angle(k)))/(angle(k+1)-angle(k))
            vec1(long2+kk) = alfa*angl+beta
        else
            vec1(long2+kk) = 0.d0
        end if
    end do
end subroutine
