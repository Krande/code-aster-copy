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

subroutine lclbr2(fami, kpg, ksp, imate, compor, &
                  ndim, epsm, t, e, sigmt, &
                  sigmc, epsic, compn, gamma)
    implicit none
#include "asterfort/rcvalb.h"
    character(len=16) :: compor(*)
    character(len=*) :: fami
    integer(kind=8) :: imate, ndim, t(3, 3), kpg, ksp
    real(kind=8) :: epsm(6), e, sigmt, sigmc, gamma, compn, epsic
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT BETON REGLEMENTAIRE : INITIALISATION
!
! IN  COMPOR     : NOM DE LA LOI DE COMPORTEMENT
! IN  IMATE      : CODE MATERIAU
! IN  EPSM       : DEFORMATION AU TEMPS MOINS
! OUT T          : TENSEUR DE PLACEMENT (PASSAGE VECT -> MATRICE)
! OUT LAMBDA
! OUT DEUXMU
! OUT ALPHA
! OUT GAMMA
! OUT SEUIL
! ----------------------------------------------------------------------
    integer(kind=8) :: icodre(6)
    character(len=16) :: nomres(6)
    real(kind=8) :: valres(6)
!
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
!    LECTURE DES CARACTERISTIQUES DU MATERIAU
    nomres(1) = 'E'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                1, nomres, valres, icodre, 1)
    e = valres(1)
!    LECTURE DES CARACTERISTIQUES D'ENDOMMAGEMENT
    nomres(1) = 'D_SIGM_EPSI'
    nomres(2) = 'SYT'
    nomres(3) = 'SYC'
    nomres(4) = 'EPSC'
    nomres(5) = 'N'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'BETON_REGLE_PR', 0, ' ', [0.d0], &
                5, nomres, valres, icodre, 1)
    gamma = -e/valres(1)
    sigmt = valres(2)
    sigmc = -valres(3)
    epsic = -valres(4)
    compn = valres(5)
!
end subroutine
