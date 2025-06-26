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
subroutine arljac(nno, ndim, dff, coor, invjac)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/matinv.h"
!
    integer(kind=8) :: nno, ndim
    real(kind=8) :: coor(ndim*nno)
    real(kind=8) :: dff(3, nno), invjac(3, 3)
!
! ----------------------------------------------------------------------
! CALCUL DE L'INVERSE DE LA JACOBIENNE EN XE
! ----------------------------------------------------------------------
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELT
! IN  COOR   : COORDONNEES DES NOEUDS DE L'ELEMENT
! IN  DFF    : DERIVEES DES FONCTION DES FORMES AU POINT XE
! IN  NDIM   : DIMENSION DE L'ESPACE
! OUT INVJAC : INVERSE DE LA JACONIENNE AU POINT XE
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k
    real(kind=8) :: jacobi(3, 3), temp(3, 3), det
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- JACOBIENNE EN XE
!
    jacobi(:, :) = 0.d0
    do i = 1, ndim
        do j = 1, ndim
            do k = 1, nno
                jacobi(i, j) = jacobi(i, j)+dff(j, k)*coor(ndim*(k-1)+i)
            end do
        end do
    end do
!
    if (ndim == 2) then
        jacobi(3, 3) = 1.d0
    else if (ndim == 1) then
        jacobi(3, 3) = 1.d0
        jacobi(2, 2) = 1.d0
    end if
!
! --- INVERSE DE LA JACOBIENNE
!
    call matinv('S', 3, jacobi, temp, det)
    do i = 1, 3
        do j = 1, 3
            invjac(i, j) = 0.d0
        end do
    end do
    do i = 1, ndim
        do j = 1, ndim
            invjac(i, j) = temp(i, j)
        end do
    end do
!
    call jedema()
!
end subroutine
