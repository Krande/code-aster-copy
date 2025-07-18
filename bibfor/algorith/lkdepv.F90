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
subroutine lkdepv(nbmat, mater, depsv, ddepsv, dgamv, &
                  ddgamv)
!
    implicit none
#include "asterfort/lcdevi.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nbmat
    real(kind=8) :: mater(nbmat, 2), depsv(6), ddepsv(6)
    real(kind=8) :: dgamv, ddgamv(6)
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE--------------------------
! =================================================================
! --- BUT : DERIVEE DE LA DEFORMATION VISQUEUSE ET DU PARAMETRE
! ---- D ECROUISSAGE VISQUEUX
! =================================================================
! IN  : NBMAT :  NOMBRE DE PARAMETRES MATERIAU --------------------
! --- : MATER :  COEFFICIENTS MATERIAU A T+DT ---------------------
! ---         :  MATER(*,1) = CARACTERISTIQUES ELASTIQUES ---------
! ---         :  MATER(*,2) = CARACTERISTIQUES PLASTIQUES ---------
! --- : DT    : PAS DE TEMPS --------------------------------------
!     : DEPSV : DEFORMATIONS VISQUEUSES ---------------------------
! OUT : DDEPSV: DEFORMATIONS DEVIATORIQUES VISQUEUSES -------------
!     : DGAMV : PARAMETRE D ECROUISSAGE VISQUEUX ------------------
!     :DDGAMV : DERIVEE DU PARAMETRE D ECROUISSAGE VISQUEUX PAR  --
!     :         RAPPORT A DEPS ------------------------------------
! =================================================================
    common/tdim/ndt, ndi
    integer(kind=8) :: i, k, ndi, ndt
    real(kind=8) :: zero, un, mun, deux, trois
    real(kind=8) :: devia(6, 6), deviat(6, 6)
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(mun=-1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! =================================================================
! --- CALCUL DU DEVIATEUR DU TENSEUR DES DEFORMATIONS VISQUEUSES -
! =================================================================
!
    call lcdevi(depsv, ddepsv)
!
! =================================================================
! --- CALCUL DE DGAMV ------------------------------------
! =================================================================
!
    dgamv = 0.d0
!
    do i = 1, ndt
!
        dgamv = dgamv+ddepsv(i)**2
!
    end do
    dgamv = sqrt(deux/trois*dgamv)
! =================================================================
! --- MATRICE DE PROJECTION DEVIATORIQUE --------------------------
! =================================================================
!
    call r8inir(6*6, 0.d0, devia, 1)
!
    do i = 1, 3
        do k = 1, 3
            devia(i, k) = mun/trois
        end do
    end do
!
    do i = 1, ndt
        devia(i, i) = devia(i, i)+un
    end do
!
    deviat(1:ndt, 1:ndt) = transpose(devia(1:ndt, 1:ndt))
! =================================================================
! --- CALCUL DE DERIVEE DE DGAMV/DEPS ---------------------------
! =================================================================
!
    call r8inir(6, 0.d0, ddgamv, 1)
!
!      SI PAS DE VISCOSITE DGAMV=0
    if (dgamv .eq. zero) then
        do i = 1, ndt
            ddgamv(i) = zero
        end do
    else
        ddgamv(1:ndt) = matmul(deviat(1:ndt, 1:ndt), ddepsv(1:ndt))
        do i = 1, ndt
            ddgamv(i) = deux/trois*ddgamv(i)/dgamv
        end do
    end if
!
! =================================================================
end subroutine
