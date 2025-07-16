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
subroutine rco3d_calcmat(nb_gauss, gauss_weight, gauss_coor, jac_det, &
                         ff_co, ff_3d, s, t, n, epai, crig, &
                         nno_co, nno_3d, skip, mat)
!
    use raco3d_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"

    real(kind=8), intent(in) :: epai, crig
    integer(kind=8), intent(in) :: nb_gauss, nno_co, nno_3d
    real(kind=8), intent(in) :: jac_det(NB_GAUSS_MAX)
    real(kind=8), intent(in) :: gauss_weight(NB_GAUSS_MAX)
    real(kind=8), intent(in) :: gauss_coor(2, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: ff_co(NB_NO_CO_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: ff_3d(NB_NO_3D_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(in) :: t(3, NB_GAUSS_MAX), n(3, NB_GAUSS_MAX), s(3)
    aster_logical, intent(in) :: skip(NB_GAUSS_MAX)
    real(kind=8), intent(out) :: mat(:, :)

! -------------------------------------------------------------------------------
!  SUBROUTINE: rco3d_calcmat
!
!  DESCRIPTION:
!  This subroutine calculates the elementary link matrix for a 3D-shell
!  connection, considering contributions from both shell and 3D elements.
!  It uses integration over Gauss points and incorporates rotation and thickness
!  effects to compute the final matrix.
!
!  INPUT PARAMETERS:
!  -----------------
!  nb_gauss         - IN    - I    - Number of Gauss points for integration.
!  gauss_weight     - IN    - R*8  - Weights associated with each Gauss point.
!                                     Dimensions: (10)
!  gauss_coor       - IN    - R*8  - Coordinates of Gauss points.
!                                     Dimensions: (2, 10)
!  jac_det          - IN    - R*8  - Determinant of the Jacobian matrix at each Gauss point.
!                                     Dimensions: (10)
!  ff_co            - IN    - R*8  - Shape functions for shell nodes.
!                                     Dimensions: (3, 10)
!  ff_3d            - IN    - R*8  - Shape functions for 3D nodes.
!                                     Dimensions: (8, 10)
!  s                - IN    - R*8  - Tangential vectors 1 at each Gauss point.
!                                     Dimensions: (3)
!  t                - IN    - R*8  - Tangential vectors 2 at each Gauss point.
!                                     Dimensions: (3, 10)
!  n                - IN    - R*8  - Normal vectors at each Gauss point.
!                                     Dimensions: (3, 10)
!  epai             - IN    - R*8  - Thickness of the shell.
!  nno_co           - IN    - I    - Number of shell nodes at the interface.
!  nno_3d           - IN    - I    - Number of 3D nodes at the interface.
!  skip             - IN   - LOG  - Logical flag indicating if the computation should be skipped.
!
!  OUTPUT PARAMETERS:
!  ------------------
!  mat              - OUT   - R*8  - Link elementary matrix.
!                                     Dimensions: (dynamic, determined during computation)
! -------------------------------------------------------------------------------

!
    real(kind=8) :: mat1(3, 6*nno_co), mat2(3, 6*nno_co+3*nno_3d)
    real(kind=8) :: rot(3, 3), mat3(3, 3), mat4(3, 3), mat5(3, 3)
    integer(kind=8) :: i, j, index, ddl_co, ddl_3d
    real(kind=8) :: epsilon

    mat1 = 0.0d0
    mat2 = 0.0d0
    mat3 = 0.0d0
    mat4 = 0.0d0
    mat5 = 0.0d0
    rot = 0.0d0

    epsilon = crig
    !
    ddl_co = 6
    ddl_3d = 3
    !
    do i = 1, nb_gauss

        if (skip(i)) then
            cycle
        end if

        ! Precompute the repeated index scaling factor

        rot(:, 1) = s(:)
        rot(:, 2) = t(:, i)
        rot(:, 3) = n(:, i)
        ! Iterate over the shell nodes
        do j = 1, nno_co
            mat3 = 0.0d0
            mat4 = 0.0d0
            mat5 = 0.0d0
            ! Fill the diagonal of mat3 with the current `ff_co(j, i)`
            mat3(1, 1) = ff_co(j, i)
            mat3(2, 2) = ff_co(j, i)
            mat3(3, 3) = ff_co(j, i)
            ! Assign `mat3` to the appropriate slice in `mat2`
            mat2(:, ddl_co*(j-1)+1:ddl_co*(j-1)+3) = mat3
            mat1(:, ddl_co*(j-1)+1:ddl_co*(j-1)+3) = mat3
            ! rotation ddls
            mat5(1, 2) = ff_co(j, i)
            mat5(2, 1) = -ff_co(j, i)
            mat5(3, 3) = epsilon*ff_co(j, i)
            mat4 = matmul(matmul(rot, mat5), transpose(rot))
            mat2(:, ddl_co*(j-1)+4:ddl_co*(j-1)+6) = &
                0.5d0*epai*gauss_coor(2, i)*mat4

            mat1(:, ddl_co*(j-1)+4:ddl_co*(j-1)+6) = &
                0.5d0*epai*gauss_coor(2, i)*mat4
        end do
        !
        index = 6*nno_co
        ! Iterate over the 3d nodes
        do j = 1, nno_3d
            mat3 = 0.0d0
            mat3(1, 1) = ff_3d(j, i)
            mat3(2, 2) = ff_3d(j, i)
            mat3(3, 3) = ff_3d(j, i)
            !
            mat2(:, index+ddl_3d*(j-1)+1:index+ddl_3d*(j-1)+3) = -mat3
        end do

        mat = mat+jac_det(i)*gauss_weight(i)*matmul(transpose(mat1), mat2)

    end do

    !do i = 1, 6*nno_co
    !   do j = 1, 6*nno_co+3*nno_3d
    !        write(*, '(F15.12)', advance="no") mat(i, j) ! Format with 2 decimal points
    !        if (j < 6*nno_co+3*nno_3d) then
    !            write(*, '(A)', advance="no") ", "
    !        else
    !            write(*, *) ","  ! Newline at the end of the row
    !        end if
    !    end do
    !end do
    !write(*,*) " END"
end subroutine
