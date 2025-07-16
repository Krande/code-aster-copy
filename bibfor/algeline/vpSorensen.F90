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
subroutine vpSorensen(mod45, matrAsse, matrGeom, &
                      optionModal, calcLevel, &
                      coefDimSpace, nbFreq, bande, &
                      eigsol)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterfort/freqom.h"
#include "asterfort/omega2.h"
#include "asterfort/vpcres.h"
!
    character(len=4), intent(in) :: mod45
    character(len=19), intent(in) :: matrAsse, matrGeom
    character(len=16), intent(in) :: optionModal, calcLevel
    integer(kind=8), intent(in) :: nbFreq
    real(kind=8), intent(in) :: bande(2)
    integer(kind=8), intent(in) :: coefDimSpace
    character(len=19), intent(in) :: eigsol
!
! --------------------------------------------------------------------------------------------------
!
! CREATION ET REMPLISSAGE DE LA SD EIGSOL pour GEP SYM REEL RESOLU VIA SORENSEN
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ibid, nbborn, nbvec2
    real(kind=8) :: r8bid, omecor, bandeLocal(2)
    character(len=1), parameter :: k1bid = 'R'
    character(len=16) :: typres
    character(len=16), parameter :: k16bid = " "
    character(len=19), parameter :: k19bid = " "
    character(len=19) :: matrA, matrB
! - Valeurs en dur
    character(len=8), parameter :: method = 'SORENSEN'
    integer(kind=8), parameter :: nbvect = 0
    integer(kind=8), parameter :: nbrss = 5
    real(kind=8), parameter :: alpha = 0.717d0
    integer(kind=8), parameter :: maxitr = 200
    real(kind=8), parameter :: tolsor = 0.d0, precsh = 5.d-2, fcorig = 1.d-2, precdc = 5.d-2
    character(len=8), parameter :: arret = 'NON'
    character(len=16), parameter :: sturm = 'NON'
    character(len=16), parameter :: modri2 = 'SANS'
    character(len=16), parameter :: stoper = 'NON'
!
! --------------------------------------------------------------------------------------------------
!

! - Set matrixes
    matrA = matrAsse
    matrB = matrGeom

! - Set parameters
    bandeLocal = bande
    ibid = isnnem()
    r8bid = r8vide()
    if (mod45 .eq. 'VIBR') then
        typres = 'DYNAMIQUE'
    else
        typres = 'MODE_FLAMB'
    end if
    if (optionModal .eq. 'BANDE') then
        nbborn = 2
    else
        nbborn = 1
    end if
    nbvec2 = coefDimSpace
    if (typres(1:9) .eq. 'DYNAMIQUE') then
        omecor = omega2(fcorig)
        bandeLocal(1) = freqom(bande(1))
        bandeLocal(2) = freqom(bande(2))
    else
        omecor = fcorig
    end if

!
    call vpcres(eigsol, typres, matrA, matrB, k19bid, &
                optionModal, method, modri2, arret, k19bid, &
                stoper, sturm, calcLevel, k1bid, k16bid, &
                nbFreq, nbvect, nbvec2, nbrss, nbborn, &
                ibid, ibid, ibid, ibid, maxitr, &
                bandeLocal, precsh, omecor, precdc, r8bid, &
                r8bid, r8bid, r8bid, r8bid, tolsor, alpha)
!
end subroutine
