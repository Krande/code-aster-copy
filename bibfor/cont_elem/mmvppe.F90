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

! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
! aslint: disable=W1504
!
subroutine mmvppe(ndim, nne, nnm, nnl, nbdm, &
                  i_reso_geom, l_large_slip, &
                  jeusup, &
                  tau1, tau2, &
                  ffe, ffm, ffl, dffm, ddffm, &
                  jeu, djeu, djeut, &
                  dlagrc, dlagrf, &
                  norm, mprojt, &
                  mprt1n, mprt2n, &
                  mprt11, mprt12, mprt21, mprt22, &
                  kappa, &
                  taujeu1, taujeu2, &
                  dnepmait1, dnepmait2)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/lcgeominit.h"
#include "asterfort/mmdepm.h"
#include "asterfort/mmgeom.h"
#include "asterfort/mmlagm.h"
#include "asterfort/mmmjeu.h"
#include "asterfort/mmcalg.h"
#include "asterfort/mmreac.h"
#include "Contact_type.h"
!
    integer(kind=8), intent(in) :: ndim, nne, nnm, nnl, nbdm
    integer(kind=8), intent(in) :: i_reso_geom
    aster_logical, intent(in) :: l_large_slip
    real(kind=8), intent(in) :: jeusup
    real(kind=8), intent(in) :: tau1(3), tau2(3)
    real(kind=8), intent(in) :: ffe(9), ffm(9), ffl(9), dffm(2, 9), ddffm(3, 9)
    real(kind=8), intent(out) :: jeu
    real(kind=8), intent(out) :: djeut(3), djeu(3), dlagrc, dlagrf(2)
    real(kind=8), intent(out) :: norm(3)
    real(kind=8), intent(out) :: mprojt(3, 3)
    real(kind=8), intent(out) :: mprt1n(3, 3), mprt2n(3, 3)
    real(kind=8), intent(out) :: mprt11(3, 3), mprt12(3, 3), mprt21(3, 3), mprt22(3, 3)
    real(kind=8), intent(out) :: kappa(2, 2)
    real(kind=8), intent(out) :: dnepmait1, dnepmait2, taujeu1, taujeu2
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Compute quantities (for vector)
!
! --------------------------------------------------------------------------------------------------
!
! In  ndim             : dimension of problem (2 or 3)
! In  nne              : number of slave nodes
! In  nnm              : number of master nodes
! In  nnl              : number of nodes with Lagrange multiplicators (contact and friction)
! In  nbdm             : number of components by node for all dof
! In  i_reso_geom      : algorithm for geometry
! In  l_large_slip     : flag for GRAND_GLISSEMENT
! In  jeusup           : gap from DIST_ESCL/DIST_MAIT
! In  tau1             : first tangent at current contact point
! In  tau2             : second tangent at current contact point
! In  ffe              : shape function for slave nodes
! In  ffm              : shape function for master nodes
! In  ffl              : shape function for Lagrange dof
! In  dffm             : first derivative of shape function for master nodes
! In  ddffm            : second derivative of shape function for master nodes
! Out jeu              : normal gap
! Out djeut            : increment of tangent gaps
! Out dlagrc           : increment of contact Lagrange from beginning of time step
! Out dlagrf           : increment of friction Lagrange from beginning of time step
! Out norm             : normal at current contact point
! Out tau1             : first tangent at current contact point
! Out tau2             : second tangent at current contact point
! Out mprojt           : matrix of tangent projection
! Out mprt1n           : projection matrix first tangent/normal
! Out mprt2n           : projection matrix second tangent/normal
! Out mprt11           : projection matrix first tangent/first tangent
!                        tau1*TRANSPOSE(tau1)(matrice 3*3)
! Out mprt12           : projection matrix first tangent/second tangent
!                        tau1*TRANSPOSE(tau2)(matrice 3*3)
! Out mprt21           : Projection matrix second tangent/first tangent
!                        tau2*TRANSPOSE(tau1)(matrice 3*3)
! Out mprt22           : Projection matrix second tangent/second tangent
!                        tau2*TRANSPOSE(tau2)(matrice 3*3)
! Out kappa            : scalar matrix for sliding kinematic
!                        KAPPA(i,j) = INVERSE[tau_i.tau_j-JEU*(ddFFM*geomm)](matrice 2*2)
! Out dnepmait1, dnepmait2, taujeu1, taujeu2
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_disp_incr, jv_disp
    real(kind=8) :: ppe
    real(kind=8) :: ddepmam(9, 3)
    real(kind=8) :: geomm(3), geome(3)
    real(kind=8) :: ddeple(3), ddeplm(3)
    real(kind=8) :: deplme(3), deplmm(3)
    real(kind=8) :: mprojn(3, 3)
    real(kind=8) :: gene11(3, 3), gene21(3, 3), gene22(3, 3)
    real(kind=8) :: a(2, 2), ha(2, 2), h(2, 2), hah(2, 2)
    real(kind=8) :: mprnt1(3, 3), mprnt2(3, 3)
    real(kind=8) :: vech1(3), vech2(3)
    real(kind=8) :: elem_slav_init(9, 3), elem_mast_init(9, 3)
    real(kind=8) :: elem_slav_coor(9, 3), elem_mast_coor(9, 3)
    real(kind=8) :: elem_slav_temp(nne, ndim), elem_mast_temp(nnm, ndim)
    real(kind=8) :: elem_slav_tem2(nne, ndim), elem_mast_tem2(nnm, ndim)
    integer(kind=8) :: i_dime, i_nne, i_nnm
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PDEPL_P', 'L', jv_disp_incr)
    call jevech('PDEPL_M', 'L', jv_disp)
!
! - Initializations
!
    elem_slav_coor = 0.d0
    elem_mast_coor = 0.d0
!
! - Coefficient to update geometry
!
    ppe = 0.d0
    if (i_reso_geom .eq. ALGO_NEWT) then
        ppe = 1.d0
    end if
!
! - Get initial geometry
!
    call lcgeominit(ndim, &
                    nne, nnm, &
                    elem_mast_temp, elem_slav_temp)
    elem_slav_init(:, :) = 0.d0
    elem_mast_init(:, :) = 0.d0
    do i_dime = 1, ndim
        do i_nne = 1, nne
            elem_slav_init(i_nne, i_dime) = elem_slav_temp(i_nne, i_dime)
        end do
        do i_nnm = 1, nnm
            elem_mast_init(i_nnm, i_dime) = elem_mast_temp(i_nnm, i_dime)
        end do
    end do
!
! - Update geometry
!
    call mmreac(ndim, nne, nnm, &
                jv_disp, jv_disp_incr, ppe, &
                elem_slav_temp, elem_mast_temp, &
                elem_slav_tem2, elem_mast_tem2, &
                nbdm_=nbdm, &
                ddepmam_=ddepmam)
    elem_slav_coor(:, :) = 0.d0
    elem_mast_coor(:, :) = 0.d0
    do i_dime = 1, ndim
        do i_nne = 1, nne
            elem_slav_coor(i_nne, i_dime) = elem_slav_tem2(i_nne, i_dime)
        end do
        do i_nnm = 1, nnm
            elem_mast_coor(i_nnm, i_dime) = elem_mast_tem2(i_nnm, i_dime)
        end do
    end do
!
! - Compute local basis
!
    call mmgeom(ndim, &
                nne, nnm, &
                ffe, ffm, &
                elem_slav_coor, elem_mast_coor, &
                tau1, tau2, &
                norm, mprojn, mprojt, &
                geome, geomm)
!
! - Compute increment of Lagrange multipliers
!
    call mmlagm(nbdm, ndim, nnl, jv_disp_incr, ffl, &
                dlagrc, dlagrf)
!
! - Compute increment of displacements
!
    call mmdepm(nbdm, ndim, &
                nne, nnm, &
                jv_disp, jv_disp_incr, &
                ffe, ffm, &
                ddeple, ddeplm, &
                deplme, deplmm)
!
! - Compute gaps
!
    call mmmjeu(ndim, i_reso_geom, jeusup, &
                geome, geomm, &
                ddeple, ddeplm, &
                norm, mprojt, &
                jeu, djeu, djeut)
!
! - Compute projection matrices for second variation of gap
!
    call mmcalg(ndim, l_large_slip, &
                nnm, dffm, ddffm, &
                elem_mast_coor, ddepmam, &
                tau1, tau2, norm, &
                jeu, djeu, &
                gene11, gene21, gene22, &
                kappa, h, &
                vech1, vech2, &
                a, ha, hah, &
                mprt11, mprt12, mprt21, mprt22, &
                mprt1n, mprt2n, mprnt1, mprnt2, &
                taujeu1, taujeu2, &
                dnepmait1, dnepmait2)
!
end subroutine
