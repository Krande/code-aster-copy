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
subroutine fmampc(nbfonc, nbptot, sigm, rampmx)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/fmdevi.h"
    integer(kind=8) :: nbfonc, nbptot
    real(kind=8) :: sigm(*), rampmx
!     NBFONC  : IN  : NOMBRE DE FONCTIONS (6 EN 3D 4 EN 2D)
!     NBPTOT  : IN  : NOMBRE DE PAS DE TEMPS DE CALCUL
!     SIGM    : IN  : VECTEUR DES CONTRAINTES EN TOUS LES PAS DE TEMPS
!     RAMPMX  : OUT : VALEUR AMPLITUDE DE CISSION
!     -----------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: i1, i2, j
    real(kind=8) :: sig(6), rampc
    real(kind=8), pointer :: dev(:) => null()
!     ------------------------------------------------------------------
!
!------- CALCUL DU DEVIATEUR -------
!
    AS_ALLOCATE(vr=dev, size=nbfonc*nbptot)
    call fmdevi(nbfonc, nbptot, sigm, dev)
!
! -------- CALCUL AMPLITUDE DE CISSION ------
!
    rampmx = 0.d0
    do i1 = 1, nbptot-1
        do i2 = i1+1, nbptot
            do j = 1, nbfonc
                sig(j) = dev(1+(i2-1)*nbfonc+j-1)-dev(1+(i1-1)*nbfonc+j-1)
            end do
            if (nbfonc .eq. 6) then
                rampc = ( &
                        sig(1)*sig(1)+sig(2)*sig(2)+sig(3)*sig(3))/2.d0+sig(4)*sig(4)+sig(5)&
                       &*sig(5)+sig(6)*sig(6 &
                                           )
            else if (nbfonc .eq. 4) then
                rampc = (sig(1)*sig(1)+sig(2)*sig(2)+sig(3)*sig(3))/2.d0+sig(4)*sig(4)
            end if
            if (rampc .gt. rampmx) rampmx = rampc
        end do
    end do
    rampmx = 1.d0/2.d0*sqrt(rampmx)
!
    AS_DEALLOCATE(vr=dev)
!
end subroutine
