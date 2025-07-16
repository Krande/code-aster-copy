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
subroutine fgrmax(ncyc, sigmin, sigmax, smin, smax)
!      COMPTAGE DES CYCLES PAR LA METHODE RAINFLOW AVEC LE CYCLE MAX
!      AU DEBUT DE CHARGEMENT
!       ----------------------------------------------------------------
!      IN  NCYC    NOMBRE  DE  CYCLE
!      IN  SIGMAX  CONTRAINTES MAXIMALES DES CYCLES APRES RAINFLOW
!          SIGMIN  CONTRAINTES MINIMALES DES CYCLES APRES RAINFLOW
!      OUT SMAX  CONTRAINTES MAXIMALES DES CYCLES APRES RAINFLOW_MAX
!          SMIN  CONTRAINTES MINIMALES DES CYCLES APRES RAINFLOW
!       ----------------------------------------------------------------
    implicit none
#include "asterfort/infniv.h"
    real(kind=8) :: sigmax(*), sigmin(*)
    real(kind=8) :: ampmax, smax(*), smin(*)
    integer(kind=8) :: ncyc, cycmax
!       ----------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifm, niv
!-----------------------------------------------------------------------
    call infniv(ifm, niv)
!
    cycmax = 1
    ampmax = sigmax(1)-sigmin(1)
!
    do i = 2, ncyc
        if ((sigmax(i)-sigmin(i)) .gt. ampmax) then
            ampmax = sigmax(i)-sigmin(i)
            cycmax = i
        end if
    end do
!
    smin(1) = sigmin(cycmax)
    smax(1) = sigmax(cycmax)
!
    do i = 2, ncyc
        if (i .lt. cycmax) then
            smin(i) = sigmin(i-1)
            smax(i) = sigmax(i-1)
        else if (i .gt. cycmax) then
            smin(i) = sigmin(i)
            smax(i) = sigmax(i)
        end if
    end do
    if (cycmax .eq. ncyc) then
        smin(ncyc) = sigmin(ncyc-1)
        smax(ncyc) = sigmax(ncyc-1)
    end if
!
!     --- IMPRESSION DES PICS EXTRAITS DE LA FONCTION ----
    if (niv .eq. 2) then
        write (ifm, *)
        write (ifm, '(1X,A)') 'PICS APRES LE COMPTAGE RAINFLOW_MAX'
        write (ifm, *)
        write (6, *) 'NOMBRE DE CYCLES = ', ncyc
        write (ifm, *)
        write (ifm, '(1X,A)') '     CHARGEMENT_MAX     CHARGEMENT_MIN'
        write (ifm, *)
        write (ifm, '(2(1X,E18.6))') (smax(i), smin(i), i=1, ncyc)
!         DO 106 I = 1,NCYC
!             WRITE (IFM,'(2(1X,E18.6))'), SMAX(I),SMIN(I)
! 106     CONTINUE
!
    end if
!
end subroutine
