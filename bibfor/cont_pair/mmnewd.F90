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

subroutine mmnewd(type_elem, nb_node, nb_dim, elem_coor, pt_coor, &
                  iter_maxi, tole_maxi, proj_dire, ksi1, ksi2, &
                  tang_1, tang_2, error, dist_, ksi1_init, ksi2_init, beta_)
!
    implicit none
!
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/mmnonf.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmtang.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: type_elem
    integer(kind=8), intent(in) :: nb_node
    integer(kind=8), intent(in) :: nb_dim
    real(kind=8), intent(in) :: elem_coor(27)
    real(kind=8), intent(in) :: pt_coor(3)
    integer(kind=8), intent(in) :: iter_maxi
    real(kind=8), intent(in) :: tole_maxi
    real(kind=8), intent(in) :: proj_dire(3)
    real(kind=8), intent(out) :: ksi1
    real(kind=8), intent(out) :: ksi2

    real(kind=8), intent(out) :: tang_1(3)
    real(kind=8), intent(out) :: tang_2(3)
    integer(kind=8), intent(out) :: error
    real(kind=8), optional, intent(out) :: dist_
    real(kind=8), optional, intent(in) :: ksi1_init, ksi2_init
    real(kind=8), optional, intent(out) :: beta_
!
! --------------------------------------------------------------------------------------------------
!
! Contact (all methods)
!
! Projection of point on element (Newton algorithm) - By given direction
!
! --------------------------------------------------------------------------------------------------
!
! In  type_elem : element type
! In  nb_node   : number of nodes of element
! In  nb_dim    : dimension of element (2 or 3)
! In  elem_coor : coordinates of nodes of the element
! In  pt_coor   : coordinates of point to project
! In  iter_maxi : Newton algorithm - Maximum number of iterations
! In  tole_maxi : Newton algorithm - Tolerance
! In  proj_dire : direction of projection
! Out ksi1      : first parametric coordinate of projection of point on element
! Out ksi2      : second parametric coordinate of projection of point on element
! Out tang_1    : first tangent of local basis for the projection of point on element
! Out tang_2    : second tangent of local basis for the projection of point on element
! Out beta      : scalar value that multiply proj_dire to obtain the projected position
! Out error     : error code
!                  0  OK
!                  1  NON-OK
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: zero
    parameter(zero=0.d0)
    aster_logical, parameter :: debug = ASTER_FALSE
!
    integer(kind=8) :: ino, idim, iter
    real(kind=8) :: ff(9), dff(2, 9)
    real(kind=8) :: vect_posi(3)
    real(kind=8) :: matrix(3, 3)
    real(kind=8) :: residu(3)
    real(kind=8) :: dksi1, dksi2, dbeta
    real(kind=8) :: det, test, refe, beta
    real(kind=8) :: tole_rela, tole_abso, tole_newt
    real(kind=8) :: dist, dist_mini, ksi1_mini, ksi2_mini, beta_mini
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_node .le. 9)
    ASSERT(nb_dim .le. 3)
    ASSERT(nb_dim .ge. 2)
!
! - Initializations
!
    error = 0
    ksi1 = zero
    ksi2 = zero
    beta = 1.d0
    iter = 0
    tole_abso = tole_maxi/100.d0
    tole_rela = tole_maxi
    dist_mini = r8gaem()
!
    if (present(ksi1_init)) then
        ksi1 = ksi1_init
        ksi2 = ksi2_init
    end if

    if (present(beta_)) then
        beta_ = beta
    end if
!
! - Newton loop
!
20  continue
!
    vect_posi(1:3) = zero
    tang_1(1:3) = zero
    tang_2(1:3) = zero
    matrix(1:3, 1:3) = zero
    residu(1:3) = zero
    dksi1 = zero
    dksi2 = zero
    dbeta = zero
!
! - Shape functions (and derivates) at current point
!
    call mmnonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                ff)
    call mmdonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                dff)
!
! - Position vector of current point
!
    do idim = 1, 3
        do ino = 1, nb_node
            vect_posi(idim) = elem_coor(3*(ino-1)+idim)*ff(ino)+vect_posi(idim)
        end do
    end do
!
! - Local base
!
    call mmtang(nb_dim, nb_node, elem_coor, dff, tang_1, &
                tang_2)
!
! - Quantity to minimize
!
    do idim = 1, 3
        vect_posi(idim) = pt_coor(idim)-vect_posi(idim)
    end do
    dist = sqrt(vect_posi(1)*vect_posi(1)+vect_posi(2)*vect_posi(2)+vect_posi(3)*vect_posi(3))

    if (present(dist_)) then
        dist_ = dist
    end if
!
! - Newton residual
!
    do idim = 1, 3
        residu(idim) = vect_posi(idim)-beta*proj_dire(idim)
    end do
!
! - Tangent matrix (Newton)
!
    do idim = 1, 3
        matrix(idim, 1) = tang_1(idim)
        if (nb_dim .eq. 2) then
            matrix(idim, 2) = proj_dire(idim)
        elseif (nb_dim .eq. 3) then
            matrix(idim, 2) = tang_2(idim)
            matrix(idim, 3) = proj_dire(idim)
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
! - System determinant
!
    if (nb_dim .eq. 2) then
        det = matrix(1, 1)*matrix(2, 2)-matrix(1, 2)*matrix(2, 1)
    else if (nb_dim .eq. 3) then
        det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))- &
              matrix(2, 1)*(matrix(1, 2)*matrix(3, 3)-matrix(3, 2)*matrix(1, 3))+ &
              matrix(3, 1)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3))
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (det .eq. 0.d0) then
        error = 1
        goto 999
    end if
!
! - Solve system
!
    if (nb_dim .eq. 2) then
        dksi1 = (residu(1)*matrix(2, 2)-residu(2)*matrix(1, 2))/det
        dksi2 = 0.d0
        dbeta = (residu(2)*matrix(1, 1)-residu(1)*matrix(2, 1))/det
    else if (nb_dim .eq. 3) then
        dksi1 = (residu(1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))+ &
                 residu(2)*(matrix(3, 2)*matrix(1, 3)-matrix(1, 2)*matrix(3, 3))+ &
                 residu(3)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3)))/det
        dksi2 = (residu(1)*(matrix(3, 1)*matrix(2, 3)-matrix(2, 1)*matrix(3, 3))+ &
                 residu(2)*(matrix(1, 1)*matrix(3, 3)-matrix(3, 1)*matrix(1, 3))+ &
                 residu(3)*(matrix(2, 1)*matrix(1, 3)-matrix(2, 3)*matrix(1, 1)))/det
        dbeta = (residu(1)*(matrix(2, 1)*matrix(3, 2)-matrix(3, 1)*matrix(2, 2))+ &
                 residu(2)*(matrix(3, 1)*matrix(1, 2)-matrix(1, 1)*matrix(3, 2))+ &
                 residu(3)*(matrix(1, 1)*matrix(2, 2)-matrix(2, 1)*matrix(1, 2)))/det
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Update
!
    ksi1 = ksi1+dksi1
    ksi2 = ksi2+dksi2
    beta = beta+dbeta

    if (present(beta_)) then
        beta_ = beta
    end if
!
!
! - Save values if Newton avoids
!
    if (dist .le. dist_mini) then
        dist_mini = dist
        ksi1_mini = ksi1
        ksi2_mini = ksi2
        beta_mini = beta
    end if

!
! - Convergence
!
    refe = (ksi1*ksi1+ksi2*ksi2+beta*beta)
    if (refe .le. tole_rela) then
        tole_newt = tole_abso
        test = sqrt(dksi1*dksi1+dksi2*dksi2+dbeta*dbeta)
    else
        tole_newt = tole_rela
        test = sqrt(dksi1*dksi1+dksi2*dksi2+dbeta*dbeta)/sqrt(refe)
    end if
!
! - Continue or not ?
!
    if ((test .gt. tole_newt) .and. (iter .lt. iter_maxi)) then
        iter = iter+1
        goto 20
    else if ((iter .ge. iter_maxi) .and. (test .gt. tole_newt)) then
        ksi1 = ksi1_mini
        ksi2 = ksi2_mini
        call mmnonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                    ff)
        call mmdonf(nb_dim, nb_node, type_elem, ksi1, ksi2, &
                    dff)
        call mmtang(nb_dim, nb_node, elem_coor, dff, tang_1, &
                    tang_2)
        error = 1
    end if
!
! - End of loop
!
999 continue
!
    if (debug .and. error .eq. 1) then
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
        write (6, *) 'BETA  : ', beta
    end if
!
end subroutine
