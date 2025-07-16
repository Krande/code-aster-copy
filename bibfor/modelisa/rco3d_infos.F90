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
subroutine rco3d_infos(typmaco, typma3d, epai, j_geom, nb_gauss, gauss_coor, &
                       gauss_weight, jac_det, ff_co, ff_3d, s, t, n, skip)
!
    use raco3d_module, only: NB_GAUSS_MAX, det_jacob, find_parametric_coordinates, &
                             NB_NO_CO_MAX, NB_NO_3D_MAX
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfd2.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/elrfno.h"
#include "asterfort/provec.h"
#include "jeveux.h"
!

!
    character(len=8), intent(in) :: typmaco, typma3d
    real(kind=8), intent(in) :: epai
    integer(kind=8), intent(in) :: j_geom
    integer(kind=8), intent(out) :: nb_gauss
    real(kind=8), intent(out) :: gauss_coor(2, NB_GAUSS_MAX)
    real(kind=8), intent(out) :: jac_det(NB_GAUSS_MAX)
    real(kind=8), intent(out) :: gauss_weight(NB_GAUSS_MAX)
    real(kind=8), intent(out) :: ff_co(NB_NO_CO_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(out) :: ff_3d(NB_NO_3D_MAX, NB_GAUSS_MAX)
    real(kind=8), intent(out) :: t(3, NB_GAUSS_MAX), n(3, NB_GAUSS_MAX), s(3)
    aster_logical, intent(out) :: skip(NB_GAUSS_MAX)
!
! -------------------------------------------------------------------------------
!  SUBROUTINE: rco3d_infos
!
!  DESCRIPTION:
!  This subroutine computes and retrieves the geometric and integration
!  information required for a 3D-shell connection. It provides the Gauss points,
!  shape functions, Jacobian determinant, and coordinate transformations
!  necessary for the elementary link matrix.
!
!  INPUT PARAMETERS:
!  -----------------
!  typmaco          - IN    - K8   - Type of the shell mesh element.
!  typma3d          - IN    - K8   - Type of the 3D mesh element.
!  epai             - IN    - R*8  - Thickness of the shell.
!  j_geom           - IN    - I    - Geometry identifier used to query geometric data.
!
!  OUTPUT PARAMETERS:
!  ------------------
!  nb_gauss         - OUT   - I    - Number of Gauss points used for integration.
!  gauss_coor       - OUT   - R*8  - Coordinates of Gauss points in the reference element.
!                                     Dimensions: (2, 10)
!  gauss_weight     - OUT   - R*8  - Weights for integration at Gauss points.
!                                     Dimensions: (10)
!  jac_det          - OUT   - R*8  - Determinant of the Jacobian matrix at Gauss points.
!                                     Dimensions: (10)
!  ff_co            - OUT   - R*8  - Shape function values for shell nodes at Gauss points.
!                                     Dimensions: (3, 10)
!  ff_3d            - OUT   - R*8  - Shape function values for 3D nodes at Gauss points.
!                                     Dimensions: (8, 10)
!  t                - OUT   - R*8  - Tangential vectors at each Gauss point.
!                                     Dimensions: (3, 10)
!  n                - OUT   - R*8  - Normal vectors at each Gauss point.
!                                     Dimensions: (3, 10)
!  s                - OUT   - R*8  - Tangential vectors at each Gauss point.
!                                     Dimensions: (3)
!  skip             - OUT   - LOG  - Logical flag indicating if the computation should be skipped.
! -------------------------------------------------------------------------------

!
    real(kind=8) :: gscoo2(NB_GAUSS_MAX), gpt_wei(NB_GAUSS_MAX), x(1)
    real(kind=8) :: ff(max(NB_NO_3D_MAX, NB_NO_CO_MAX)), df(3, NB_NO_CO_MAX)
    real(kind=8) :: df_co(NB_NO_CO_MAX, NB_GAUSS_MAX)
    real(kind=8) :: ddf(3, 3, NB_NO_CO_MAX), ddf_co(NB_NO_CO_MAX, NB_GAUSS_MAX)
    real(kind=8) :: coor_pt_co(3, NB_NO_CO_MAX), coor_pt_3d(3, NB_NO_3D_MAX)
    real(kind=8) :: coor_gp_cartesian(3, NB_GAUSS_MAX)
    real(kind=8) :: v1(3), v2(3), v3(3)
    real(kind=8) :: magnitude, eps
    integer(kind=8) :: i, j, k, l, idx, ndim, nbg, nno_co, nno_3d, nno, dim, dimd
    character(len=8) :: elrefa_co, elrefa_3d

    real(kind=8) :: res(3)

    ! initializations
    t = 0.d0
    s = 0.d0
    n = 0.d0
    df = 0.0d0
    ddf = 0.0d0
    df_co = 0.d0
    ff_co = 0.d0
    ff_3d = 0.d0
    jac_det = 0.d0
    gauss_weight = 0.d0
    gauss_coor = 0.0d0
    coor_gp_cartesian = 0.d0
    coor_pt_3d = 0.d0
    coor_pt_co = 0.d0
    eps = 1.0d-10

    skip = .true.

    ! some information about the shell/3d interface element
    ndim = 1

    elrefa_co = typmaco(1:3)
    call elrfno(elrefa_co, nno_co)

    elrefa_3d = typma3d(1:3)
    call elrfno(elrefa_3d, nno_3d)

    if (elrefa_co(1:3) .eq. 'SE2') then
        call elraga(elrefa_co, 'FPG2', ndim, nbg, gscoo2, gpt_wei)
    end if

    if (elrefa_co(1:3) .eq. 'SE3') then
        call elraga(elrefa_co, 'FPG4', ndim, nbg, gscoo2, gpt_wei)
    end if

    ! coordinates of the points of the shell/3d interface element
    ndim = 3
    idx = 0
    do i = 1, nno_co
        do j = 1, ndim
            coor_pt_co(j, i) = zr(j_geom-1+idx+j)
        end do
        idx = idx+ndim
    end do

    do i = 1, nno_3d
        do j = 1, ndim
            coor_pt_3d(j, i) = zr(j_geom-1+idx+j)
        end do
        idx = idx+ndim
    end do

    ! normal vector to the surface element in the 3d section
    do i = 1, ndim
        v1(i) = coor_pt_3d(i, 2)-coor_pt_3d(i, 1)
        v2(i) = coor_pt_3d(i, 3)-coor_pt_3d(i, 1)
    end do
    call provec(v1, v2, s)
    magnitude = sqrt(sum(s(:)**2))
    s(:) = s(:)/magnitude

    ! retrieve gauss points in the reference configuration
    ! along with the shape functions and their derivatives
    nb_gauss = nbg*nbg
    dimd = NB_NO_CO_MAX*3*3
    !
    idx = 0
    do i = 1, nbg
        x(1) = gscoo2(i)
        call elrfvf(elrefa_co, x, ff, nno)
        ASSERT(nno .eq. nno_co)
        call elrfdf(elrefa_co, x, df, nno, dim)
        ASSERT(nno .eq. nno_co)
        call elrfd2(elrefa_co, x, dimd, ddf, nno, dim)
        ASSERT(nno .eq. nno_co)
        ASSERT(dim .eq. 1)
        do j = 1, nbg
            idx = idx+1
            gauss_coor(1, idx) = gscoo2(i)
            gauss_coor(2, idx) = gscoo2(j)
            gauss_weight(idx) = gpt_wei(i)*gpt_wei(j)
            do k = 1, nno_co
                ff_co(k, idx) = ff(k)
                do l = 1, dim
                    df_co(k, idx) = df(l, k)
                    ddf_co(k, idx) = ddf(l, l, k)
                end do
            end do
        end do
    end do

    ! construction of local bases at gauss points
    ! vectors "t" and "s" are tangential to  the shell section
    ! the local basis (s, t, n) at gauss point (i) is direct,
    ! with "n" aligned along the thickness direction
    ndim = 3
    do i = 1, nb_gauss
        do l = 1, ndim
            do k = 1, nno_co
                t(l, i) = t(l, i)+df_co(k, i)*coor_pt_co(l, k)
            end do
        end do
        ! tangent t
        magnitude = sqrt(sum(t(:, i)**2))
        ! check before continue
        call provec(s, t(:, i), v3)
        ! If cross product result is zero, skip the subroutine
        if (sqrt(sum(v3(:)**2)) .le. eps*magnitude) then
            return
        end if
        ! use this opportunity to calculate the jacobian
        jac_det(i) = det_jacob(t(:, i), s, ddf_co(:, i), &
                               gauss_coor(2, i), coor_pt_co, nno_co, epai)
        ! normalize
        t(:, i) = t(:, i)/magnitude
        ! tangent s
        call provec(s, t(:, i), n(:, i))
        magnitude = sqrt(sum(n(:, i)**2))
        n(:, i) = n(:, i)/magnitude
    end do

    ! calculation of gauss point coordinates
    ! in the cartesian (global) reference frame
    ndim = 3
    do i = 1, nb_gauss
        do l = 1, ndim
            do k = 1, nno_co
                coor_gp_cartesian(l, i) = coor_gp_cartesian(l, i) &
                                          +ff_co(k, i)*coor_pt_co(l, k)
            end do
            coor_gp_cartesian(l, i) = coor_gp_cartesian(l, i) &
                                      +0.5*epai*gauss_coor(2, i)*n(l, i)
        end do
        ! retrieve the parametric coordinates at the surface of the 3d element
        res = find_parametric_coordinates(coor_gp_cartesian(:, i), &
                                          typma3d, coor_pt_3d)
        ! if the projection exists then get the form
        !functions of the 3d part (at the gauss point i)
        if (res(3) .gt. 0.0d0) then
            skip(i) = .false.
            call elrfvf(elrefa_3d, res(1:2), ff, nno)
            ASSERT(nno .eq. nno_3d)
            do k = 1, nno_3d
                ff_3d(k, i) = ff(k)
            end do
        end if
    end do

end subroutine
