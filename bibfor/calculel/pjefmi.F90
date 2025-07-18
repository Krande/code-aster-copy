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
subroutine pjefmi(elrefp, nnop, coor, xg, ndim, &
                  x1, x2, lext, xmi, distv)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
    character(len=8) :: elrefp
    integer(kind=8) :: nnop, ndim
    real(kind=8) :: coor(ndim*nnop)
    real(kind=8) :: xg(ndim), x1(ndim), x2(ndim), xmi(ndim)
    real(kind=8), intent(out) :: distv
!
! ----------------------------------------------------------------------
! but :
!   Determiner les meilleures coordonnees barycentriques entre x1 et x2
!   permet de verifier que la routine reereg a ameliore la precision
!   des coordonnees barycentriques grossieres de la 1ere etape de
!   proj_champ.
!   Retourner la distance entre le point et son projete.
! ----------------------------------------------------------------------
!
!
! in  elrefp : type de l'element
! in  nnop   : nombre de noeuds de l'element
! in  coor   : coordonnees ds espace reel des noeuds de l'element
! in  xg     : coordonnees du point dans l'espace reel
! in  ndim   : dimension de l'espace
! in  x1     : coordonnees du point 1 dans l'espace para de l'element
!              (resultat de la projection "grossiere")
! in  x2     : coordonnees du point 2 dans l'espace para de l'element
!              (resultat du raffinement reereg.f)
! in  lext   : le point x2 est "exterieur" a la maille
! out xmi    : "x mieux" : recopie de x1 ou x2 (selon le cas)
! out distv  : distance entre le point 1 et son projete
!
! ----------------------------------------------------------------------
    integer(kind=8) :: nbnomx
    parameter(nbnomx=27)
    real(kind=8) :: xr1(3), xr2(3), d1, d2
    real(kind=8) :: ff(nbnomx)
    integer(kind=8) :: k, idim, ino, nno
    aster_logical :: lext
! ----------------------------------------------------------------------
!
!     -- Si le point est exterieur, on ne tient pas compte de x2
!        => on choisit xmi=x1
!     -----------------------------------------------------------
!
!
!   -- calcul de xr1 : geometrie reelle de x1 :
!   --------------------------------------------
    call elrfvf(elrefp, x1, ff, nno)
    ASSERT(nno .eq. nnop)
    xr1(1:ndim) = 0.d0
    do idim = 1, ndim
        do ino = 1, nno
            xr1(idim) = xr1(idim)+ff(ino)*coor(ndim*(ino-1)+idim)
        end do
    end do
!
!
!   -- calcul de xr2 : geometrie reelle de x2 :
!   --------------------------------------------
    if (.not. lext) then
        call elrfvf(elrefp, x2, ff, nno)
        xr2(1:ndim) = 0.d0
        do idim = 1, ndim
            do ino = 1, nno
                xr2(idim) = xr2(idim)+ff(ino)*coor(ndim*(ino-1)+idim)
            end do
        end do
    end if
!
!
!   -- calcul de distv
!   -- quelle est la meilleure approximation de xg ?
!   -------------------------------------------------
    d1 = 0.d0
    d2 = 0.d0
    do k = 1, ndim
        d1 = d1+(xr1(k)-xg(k))**2
        if (.not. lext) d2 = d2+(xr2(k)-xg(k))**2
    end do
!
    if (lext) then
        xmi(1:ndim) = x1(1:ndim)
        distv = sqrt(d1)
    else
        if (d1 .le. d2) then
            xmi(1:ndim) = x1(1:ndim)
            distv = sqrt(d1)
        else
            xmi(1:ndim) = x2(1:ndim)
            distv = sqrt(d2)
        end if
    end if
!
end subroutine
