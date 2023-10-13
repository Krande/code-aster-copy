! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine nlcomp(phenom, imate, icamas, ndim, coorpg, time, Kglo)
!.
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/matrRotLGTher.h"
!
    character(len=32), intent(in) :: phenom
    integer, intent(in) :: imate, icamas, ndim
    real(kind=8), intent(in) :: coorpg(3), time
    real(kind=8), intent(out) :: Kglo(3, 3)
!
    integer :: j, nbres, kpg, spt
    parameter(nbres=3)
    integer :: icodre(nbres)
    character(len=8) :: fami, poum
    character(len=16) :: nomres(nbres)
    real(kind=8) :: lambor(3), lambda
    real(kind=8) ::  p(3, 3), Kloc(3, 3), valres(nbres), valpar(nbres)
    aster_logical :: aniso
!
    Kglo = 0.d0
!
! ------- EVALUATION DE LA CONDUCTIVITE LAMBDA
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    valpar(1) = time
    if (phenom .eq. 'THER') then
        nomres(1) = 'LAMBDA'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'INST', [valpar], &
                    1, nomres, valres, icodre, 1)
        lambda = valres(1)
        aniso = ASTER_FALSE
    else if (phenom .eq. 'THER_ORTH') then
        nomres(1) = 'LAMBDA_L'
        nomres(2) = 'LAMBDA_T'
        nomres(3) = 'LAMBDA_N'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'INST', [valpar], &
                    3, nomres, valres, icodre, 1)
        lambor(1) = valres(1)
        lambor(2) = valres(2)
        lambor(3) = valres(3)
        aniso = ASTER_TRUE
    else
        call utmess('F', 'ELEMENTS2_63')
    end if
!
! ------- TRAITEMENT DE L ANISOTROPIE
!
    if (aniso) then
        call matrRotLGTher(aniso, icamas, ndim, coorpg, p)
        Kloc = transpose(p)
        do j = 1, ndim
            Kloc(j, 1:3) = lambor(j)*Kloc(j, 1:3)
        end do
        Kglo = matmul(p, Kloc)
    else
        do j = 1, ndim
            Kglo(j, j) = lambda
        end do
    end if
!
end subroutine
