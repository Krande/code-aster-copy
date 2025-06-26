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

subroutine nlcomp(phenom, fami, kpg, imate, ndim, coorpg, time, tp, Kglo, dtp_, fluglo_)
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
    character(len=16), intent(in) :: phenom
    character(len=8), intent(in) :: fami
    integer(kind=8), intent(in) :: imate, ndim, kpg
    real(kind=8), intent(in) :: coorpg(3), time, tp
    real(kind=8), intent(out) :: Kglo(3, 3)
    real(kind=8), optional, intent(in) :: dtp_(3)
    real(kind=8), optional, intent(out) :: fluglo_(3)
!
    integer(kind=8), parameter :: spt = 1
    character(len=8), parameter :: poum = "+"
    integer(kind=8) :: j, nbres
    parameter(nbres=3)
    integer(kind=8) :: icodre(nbres)
    character(len=16) :: nomres(nbres)
    real(kind=8) :: lambor(3), lambda
    real(kind=8) ::  p(3, 3), Kloc(3, 3), valres(1)
    aster_logical :: aniso
!
    Kglo = 0.d0
!
! ------- EVALUATION DE LA CONDUCTIVITE LAMBDA
!
    if (phenom .eq. 'THER') then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'INST', [time], &
                    1, 'LAMBDA', valres, icodre, 1)
        lambda = valres(1)
        aniso = ASTER_FALSE
    else if (phenom .eq. 'THER_ORTH') then
        nomres(1) = 'LAMBDA_L'
        nomres(2) = 'LAMBDA_T'
        nomres(3) = 'LAMBDA_N'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'INST', [time], &
                    3, nomres, lambor, icodre, 1)
        aniso = ASTER_TRUE
    else if (phenom .eq. 'THER_NL') then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'TEMP', [tp], &
                    1, 'LAMBDA', valres, icodre, 1)
        lambda = valres(1)
        aniso = ASTER_FALSE
    else if (phenom .eq. 'THER_NL_ORTH') then
        nomres(1) = 'LAMBDA_L'
        nomres(2) = 'LAMBDA_T'
        nomres(3) = 'LAMBDA_N'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'TEMP', [tp], &
                    3, nomres, lambor, icodre, 1)
        aniso = ASTER_TRUE
    else if (phenom .eq. 'THER_HYDR') then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, 'TEMP', [tp], &
                    1, 'LAMBDA', valres, icodre, 1)
        lambda = valres(1)
        aniso = ASTER_FALSE
    else
        call utmess('F', 'THERMIQUE1_1')
    end if
!
! ------- TRAITEMENT DE L ANISOTROPIE
!
    if (aniso) then
        call matrRotLGTher(aniso, ndim, coorpg, p)
        Kloc = transpose(p)
        do j = 1, ndim
            Kloc(j, 1:3) = lambor(j)*Kloc(j, 1:3)
        end do
        Kglo = matmul(p, Kloc)
        if (present(fluglo_)) then
            ASSERT(present(dtp_))
            fluglo_ = matmul(Kglo, dtp_)
        end if
    else
        do j = 1, ndim
            Kglo(j, j) = lambda
        end do
        if (present(fluglo_)) then
            ASSERT(present(dtp_))
            fluglo_ = lambda*dtp_
        end if
    end if
!
end subroutine
