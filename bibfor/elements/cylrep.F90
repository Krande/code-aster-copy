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

subroutine cylrep(ndim, x, axe_z, orig, pgcyl, &
                  ipaxe)
!
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
!
    integer(kind=8), intent(in) :: ndim
    real(kind=8), dimension(:), intent(in) :: x, axe_z, orig
    real(kind=8), dimension(:, :), intent(inout) :: pgcyl
    integer(kind=8), intent(inout), optional :: ipaxe
!   ---------------------------------------------------------------------
!         CALCUL DE LA MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE CYLINDRIQUE
!         AU POINT X
!   ------------------------------------------------------------------
!     IN    NDIM                I  dimension du problème (2 ou 3)
!     IN    X(3)                R  coordonnées du point X
!     IN    AXE_Z(3), ORIG(3)   R  axe z et origine du repère cylindrique
!     INOUT PGCYL(3,3)          R  matrice de passage du repère global
!                                  au repère cylindrique
!     INOUT IPAXE              I  compteur (pt sur l'axe z => ipaxe = ipaxe + 1)
!
!     PGCYL est telle que V_{global} = PGCYL * V_{cyl}
!   -------------------------------------------------------------------
    real(kind=8) :: xnorm
    real(kind=8), dimension(3) :: axe_r, axe_t
!
    pgcyl(:, :) = 0.d0
!
!   Calcul du premier vecteur axe_r du repère cylindrique
!
    axe_r(:) = x(:)-orig(:)
!   Pour être conforme à chrpel
    if (ndim == 2) then
        axe_r(3) = 0.0d0
    end if
    axe_r(:) = axe_r(:)-dot_product(axe_r, axe_z)*axe_z(:)
!   Pour être conforme à chrpel
    if (ndim == 2) then
        axe_r(3) = 0.0d0
    end if
    call normev(axe_r, xnorm)
!
    if (xnorm .lt. r8prem()) then
!   si le point x appartient à l'axe z, alors l'axe r n'est pas défini par
!   le calcul précédent
!   on prend un axe arbitraire orthogonal
!   à l'axe z
        if (axe_z(1) .ne. 0.d0 .or. axe_z(2) .ne. 0.d0) then
            axe_r(1) = axe_z(2)
            axe_r(2) = -axe_z(1)
            axe_r(3) = 0.0d0
        else
!   si axe_z = ez, axe_r = ex
            axe_r(1) = 1.0d0
            axe_r(2) = 0.0d0
            axe_r(3) = 0.0d0
        end if
        if (present(ipaxe)) then
            ipaxe = ipaxe+1
        end if
    end if
!  Calcul du second vecteur axe_t du repère cylindrique
!  axe_t = axe_z x axe_r
    call provec(axe_z, axe_r, axe_t)
    call normev(axe_t, xnorm)
!   pgcyl:  (e_r, e_z, e_{\theta})
    pgcyl(:, 1) = axe_r
    pgcyl(:, 2) = axe_z
    pgcyl(:, 3) = axe_t
!  pour tests
! pgcyl(:,1) = axe_z
! pgcyl(:,2) = axe_t
! pgcyl(:,3) = -axe_r
!
end subroutine cylrep
