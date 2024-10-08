! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine lcpopl(loi, angmas, nmat, materd, materf, &
                  mod, deps, sigd, sigf, vind, &
                  vinf)
    implicit none
!     ROUTINE DE POST-TRAITEMENT POUR CERTAINES LOIS
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!     IN  NMAT    :  NOMBRE DE PARAMETRES MATERIAU INELASTIQUE
!         MATERD :  COEFFICIENTS MATERIAU A T
!         MATERF :  COEFFICIENTS MATERIAU A T+DT
!         MOD    :  TYPE DE MODELISATION
!         ANGMAS :  ANGLES NAUTIQUES (AFFE_CARA_ELEM)
!     OUT SIGF   :  CONTRAINTE A T+DT
!         VINF   :  VARIABLES INTERNES A T+DT
!     ----------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/hujori.h"
#include "asterfort/lgldcm.h"
#include "asterfort/utmess.h"
    integer :: nmat
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), sigf(*), vind(*), vinf(*)
    real(kind=8) :: angmas(3), sigd(6), deps(6)
    character(len=8) :: mod
    character(len=16) :: loi
!
    real(kind=8) :: bid66(6, 6), hill, dsig(6), nsig, neps
    real(kind=8) :: zero, un, deux, dix
    aster_logical :: reorie
    integer :: i, ndt
!
    parameter(ndt=6)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(dix=1.d1)
!
    if (loi(1:6) .eq. 'LAIGLE') then
        call lgldcm(nmat, materf, sigf, vinf)
    end if
!
! --  CONTRAINTES PLANES
    if (mod(1:6) .eq. 'C_PLAN') sigf(3) = 0.d0
!
    if (loi .eq. 'HAYHURST') then
        materd(1, 1) = materd(1, 1)*(1.0d0-vind(11))
        materf(1, 1) = materf(1, 1)*(1.0d0-vinf(11))
    end if
    if (loi .eq. 'VENDOCHAB') then
        materd(1, 1) = materd(1, 1)*(1.0d0-vind(9))
        materf(1, 1) = materf(1, 1)*(1.0d0-vinf(9))
    end if
!
end subroutine
