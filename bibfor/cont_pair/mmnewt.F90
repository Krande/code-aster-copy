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
subroutine mmnewt(type_elem, nb_node, nb_dim, elem_coor, pt_coor, &
                  iter_maxi, tole_maxi, ksi1, ksi2, tang_1, &
                  tang_2, error, l_reli_)
!
    implicit none
!
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mmfonf.h"
#include "asterfort/mmreli.h"
#include "asterfort/mmtang.h"
!
    character(len=8), intent(in) :: type_elem
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_dim
    real(kind=8), intent(in) :: elem_coor(27)
    real(kind=8), intent(in) :: pt_coor(3)
    integer(kind=8), intent(in) :: iter_maxi
    real(kind=8), intent(in) :: tole_maxi
    real(kind=8), intent(out) :: ksi1, ksi2
    real(kind=8), intent(out) :: tang_1(3), tang_2(3)
    integer(kind=8), intent(out) :: error
    aster_logical, intent(in), optional :: l_reli_
!
! --------------------------------------------------------------------------------------------------
!
! Contact (all methods)
!
! Projection of point on element (Newton algorithm) - Minimum distance
!
! --------------------------------------------------------------------------------------------------
!
! In  type_elem : element type
! In  nb_node   : number of nodes of element
! In  nb_dim    : dimension of element (2 or 3)
! In  elem_coor : coordinates of nodes of the element
! In  pt_coor   : coordinates of poitn to project
! In  iter_maxi : Newton algorithm - Maximum number of iterations
! In  tole_maxi : Newton algorithm - Tolerance
! Out ksi1      : first parametric coordinate of projection of point on element
! Out ksi2      : second parametric coordinate of projection of point on element
! Out tang_1    : first tangent of local basis for the projection of point on element
! Out tang_2    : second tangent of local basis for the projection of point on element
! Out error     : error code
!                  0  OK
!                  1  NON-OK
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: ino, idim, iter
    real(kind=8) :: ff(9), dff(2, 9), ddff(3, 9)
    real(kind=8) :: vect_posi(3)
    real(kind=8) :: matrix(2, 2), par11(3), par12(3), par22(3)
    real(kind=8) :: residu(2)
    real(kind=8) :: dksi1, dksi2
    real(kind=8) :: det, test, refe
    real(kind=8) :: alpha
    real(kind=8) :: tole_rela, tole_abso, tole_newt
    real(kind=8) :: dist, dist_mini, ksi1_mini, ksi2_mini
    aster_logical :: l_reli
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .le. 9)
    ASSERT(nb_dim .le. 3)
    ASSERT(nb_dim .ge. 2)

! - Initializations
    error = 0
    ksi1 = zero
    ksi2 = zero
    iter = 0
    tole_abso = tole_maxi/100.d0
    tole_rela = tole_maxi
    alpha = 1.d0
    dist_mini = r8gaem()
    if (present(l_reli_)) then
        l_reli = l_reli_
    else
        l_reli = ASTER_TRUE
    end if

! - Newton loop
20  continue
!
    vect_posi = zero
    tang_1 = zero
    tang_2 = zero
    par11 = zero
    par12 = zero
    par22 = zero
    matrix = zero
    residu = zero
    dksi1 = zero
    dksi2 = zero

! - Shape functions (and derivates) at current point
    call mmfonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                ff, dff, ddff)

! - Position vector of current point
    do idim = 1, 3
        do ino = 1, nb_node
            vect_posi(idim) = elem_coor(3*(ino-1)+idim)*ff(ino)+vect_posi(idim)
        end do
    end do

! - Local base
    call mmtang(nb_dim, nb_node, elem_coor, dff, tang_1, &
                tang_2)

! - Quantity to minimize
    do idim = 1, 3
        vect_posi(idim) = pt_coor(idim)-vect_posi(idim)
    end do
    dist = sqrt(vect_posi(1)*vect_posi(1)+vect_posi(2)*vect_posi(2)+vect_posi(3)*vect_posi(3))

! - Newton residual
    residu(1) = vect_posi(1)*tang_1(1)+vect_posi(2)*tang_1(2)+vect_posi(3)*tang_1(3)
    if (nb_dim .eq. 3) then
        residu(2) = vect_posi(1)*tang_2(1)+vect_posi(2)*tang_2(2)+vect_posi(3)*tang_2(3)
    end if

! - Local curvatures
    if (type_elem .ne. 'QU4' .and. type_elem .ne. 'QU8' .and. type_elem .ne. 'QU9') then
        do idim = 1, nb_dim
            do ino = 1, nb_node
                par11(idim) = elem_coor(3*(ino-1)+idim)*ddff(1, ino)+par11(idim)
                if (nb_dim .eq. 3) then
                    par22(idim) = elem_coor(3*(ino-1)+idim)*ddff(2, ino)+par22(idim)
                    par12(idim) = elem_coor(3*(ino-1)+idim)*ddff(3, ino)+par12(idim)
                end if
            end do
        end do
    end if

! - Tangent matrix (Newton)
    do idim = 1, 3
        matrix(1, 1) = -tang_1(idim)*tang_1(idim)+par11(idim)*vect_posi(idim)+matrix(1, 1)
        if (nb_dim .eq. 3) then
            matrix(1, 2) = -tang_2(idim)*tang_1(idim)+par12(idim)*vect_posi(idim)+matrix(1, 2)
            matrix(2, 1) = -tang_1(idim)*tang_2(idim)+par12(idim)*vect_posi(idim)+matrix(2, 1)
            matrix(2, 2) = -tang_2(idim)*tang_2(idim)+par22(idim)*vect_posi(idim)+matrix(2, 2)
        end if
    end do

! - System determinant
    if (nb_dim .eq. 2) then
        det = matrix(1, 1)
    else if (nb_dim .eq. 3) then
        det = matrix(1, 1)*matrix(2, 2)-matrix(1, 2)*matrix(2, 1)
    end if
!
    if (abs(det) .le. r8prem()) then
        error = 1
        goto 999
    end if

! - Solve system
    if (nb_dim .eq. 2) then
        dksi1 = -residu(1)/matrix(1, 1)
        dksi2 = 0.d0
    else if (nb_dim .eq. 3) then
        dksi1 = (matrix(2, 2)*(-residu(1))-matrix(1, 2)*(-residu(2)))/det
        dksi2 = (matrix(1, 1)*(-residu(2))-matrix(2, 1)*(-residu(1)))/det
    else
        ASSERT(ASTER_FALSE)
    end if

! - Line search
    if (l_reli) then
        call mmreli(type_elem, nb_node, nb_dim, elem_coor, pt_coor, &
                    ksi1, ksi2, dksi1, dksi2, alpha)
    else
        alpha = 1.d0
    end if

! - Update
    ksi1 = ksi1+alpha*dksi1
    ksi2 = ksi2+alpha*dksi2

! - Save values if Newton avoids
    if (dist .le. dist_mini) then
        dist_mini = dist
        ksi1_mini = ksi1
        ksi2_mini = ksi2
    end if

! - Convergence
    refe = (ksi1*ksi1+ksi2*ksi2)
    if (refe .le. tole_rela) then
        tole_newt = tole_abso
        test = sqrt(dksi1*dksi1+dksi2*dksi2)
    else
        tole_newt = tole_rela
        test = sqrt(dksi1*dksi1+dksi2*dksi2)/sqrt(refe)
    end if

! - Continue or not ?
    if ((test .gt. tole_newt) .and. (iter .lt. iter_maxi)) then
        iter = iter+1
        goto 20
    else if ((iter .ge. iter_maxi) .and. (test .gt. tole_newt)) then
        ksi1 = ksi1_mini
        ksi2 = ksi2_mini
        call mmfonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                    ff, dff, ddff)
        call mmtang(nb_dim, nb_node, elem_coor, dff, tang_1, &
                    tang_2)
        error = 1
    end if

! - End of loop
999 continue
!
    if (error .eq. 1) then
        write (6, *) 'POINT A PROJETER : ', pt_coor(1), pt_coor(2), pt_coor(3)
        write (6, *) 'MAILLE             ', type_elem, nb_node
!
        do ino = 1, nb_node
            write (6, *) '  NOEUD ', ino
            write (6, *) '   (X,Y,Z)', elem_coor(3*(ino-1)+1), &
                elem_coor(3*(ino-1)+2), &
                elem_coor(3*(ino-1)+3)
        end do
        write (6, *) 'KSI   : ', ksi1, ksi2
        write (6, *) 'ALPHA : ', alpha
    end if
!
end subroutine
