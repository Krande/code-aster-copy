! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine te0555(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module, only: behaviourOption
!
    use c_interface_tria_mitc_j
    use iso_c_binding

    implicit none
!
#include "asterfort/elrefe_info.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmdlog.h"
#include "asterfort/nmgpfi.h"
#include "asterfort/nmgrla.h"
#include "asterfort/nmplxd.h"
#include "asterfort/nmtstm.h"
#include "asterfort/rcangm.h"
#include "asterfort/tecach.h"
#include "asterfort/tgveri.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "jeveux.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/dxroep.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/utpvgl.h"
#include "asterfort/dxtpgl.h"
!
!
    character(len=16), intent(in) :: option, nomte
!
! ---------------------------------------------------------
!
! Elementary computation
!
! Elements: PLAT_REIS_MIND
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! ----------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! ----------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
!
    character(len=8) :: typmod(2)
    character(len=4) :: fami
    integer(kind=8) :: sz_tens, ndim, jpres
    integer(kind=8) :: nno, npg, imatuu, lgpg, iret
    integer(kind=8) :: igeom, imate, i, j
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: iinstm, iinstp, ideplm, ideplp, icompo, icarcr
    integer(kind=8) :: ivectu, icontp, ivarip
    integer(kind=8) :: ivarix, jv_mult_comp
    integer(kind=8) :: jtab(7)
    integer(kind=8) :: old_order(21), new_order(21)
    real(kind=8) :: angl_naut(7), pgl(3, 3), xyzl(3, 4)
    aster_logical :: matsym
    character(len=16) :: mult_comp, defo_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: codret
    integer(kind=8) :: jv_codret
!
    real(c_double) :: cst(4), coor(18), w(9), kappa, w_0(21), pres
    real(c_double) :: cdofs_f(9), F_elem_int(21)
    ! Tableau de sortie (18*18 éléments -> matrice de rigidité)
    real(c_double) :: A0(21*21), A1(21*21), A2(21*21), A3(21*21), F_elem(21*21), A_int(21*21)
    real(c_double) :: A5(21*21), A6(21*21), A7(21*21), A8(21*21), A9(21*21), A10(21*21), A_int0(21*21)
    real(c_double) :: A11(21*21), A12(21*21), A13(21*21)
    integer(c_int) :: nw, ncst, ncd, nk, ne0, ne1, ne2, nwinit, np0, np1, quadrature_permutation1(1)
    integer(c_int) :: entities0(1), entities1(1), entities2(1), quadrature_permutation0(1)
!
    real(kind=8) :: b(486), btdb(81, 81), bint(21,21),bint0(21,21) ,bint_perm(21,21),temp_mat(21,21)
    real(kind=8) :: e, nu, temp, epais, rho
    integer(kind=8) :: elas_id, igau, ipoids, nbinco, npg1
    integer(kind=8) :: nnos, ivf, idfde
    character(len=16) :: elas_keyword
    integer(kind=8) :: perm(18), n, k, reorder(21)
! --------------------------------------------------------------------
! - Finite element informations
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Initializations
!
    ndim = 2
    nbinco = ndim*nno
    cst = 0.d0
    coor = 0.d0
    w = 0.d0
    btdb(:, :) = 0.d0
    w_0 = 0.d0
    A_int = 0.d0
    A0 = 0.d0
    bint = 0.d0
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)

    do igau = 1, npg1
! ----- Get elastic parameters (only isotropic elasticity)
!
        call get_elas_id(zi(imate), elas_id, elas_keyword)
        call get_elas_para(fami, zi(imate), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           e_=e, nu_=nu)
! ----- Fill integration weight vector (Divided by 4 for FEniCS)
!
    end do

    call dxroep(rho, epais)

! - Fill material parameters vector
    kappa = 5.0/6.0
    cst(1) = e
    cst(2) = nu
    cst(3) = kappa
    cst(4) = epais
! - Fill material entities vector
! on définit 3 entities car 3 arêtes dans le triangle
    entities0(1) = 0
    entities1(1) = 1
    entities2(1) = 2

! Remplissage du vecteur de coordonnées (3 coordonnées par nœud)
    do i = 0, 5
        coor(3*i+1) = zr(igeom+3*i)
        coor(3*i+2) = zr(igeom+3*i+1)
        coor(3*i+3) = zr(igeom+3*i+2)
    end do

    ! Remplissage des degrés (pas de permut à faire ici)
    cdofs_f(1) = coor(1)
    cdofs_f(2) = coor(2)
    cdofs_f(3) = coor(3)

    cdofs_f(4) = coor(4)
    cdofs_f(5) = coor(5)
    cdofs_f(6) = coor(6)

    cdofs_f(7) = coor(7)
    cdofs_f(8) = coor(8)
    cdofs_f(9) = coor(9)
!
    nwinit = size(w_0)
    ncd = size(cdofs_f)
    ne0 = 1
    ne1 = 1
    ne2 = 1
    np0 = 1
    np1 = 1
    ncst = size(cst)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MATRICE DE RIGIDITÉ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel de la fonction C++ - Matrice de rigidité -> Part1
    call BP4_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A0)
! Appel de la fonction C++ - Matrice de rigidité -> Part2
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A1)
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities1, ne1, cst, ncst, A2)
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities2, ne2, cst, ncst, A3)

! Remplissage de la matrice intermédiaire
    do i = 1, 441
        A_int(i) = A0(i)+A1(i)+A2(i)+A3(i)
    end do
! Remplissage de la matrice K à partir de A_int
    do i = 1, 21
        do j = 1, 21
            bint(i, j) = A_int((j-1)*21+i)
        end do
    end do

![w1, θ_x1, θ_y1,w2, θ_x2, θ_y2,w3, θ_x3, θ_y3,θ_x6, θ_y6,γ_r3, p3,
!θ_x4, θ_y4, γ_r1, p1, θ_x5, θ_y5,γ_r2, p2]

    reorder = (/13, 1, 2, &
                14, 3, 4, &
                15, 5, 6, &
                11, 12, &
                18, 21, &
                7, 8, &
                16, 19, &
                9, 10, &
                17, 20 &
                /)

    do i = 1, 21
        do j = 1, 21
            bint_perm(i, j) = bint(reorder(i), reorder(j))
        end do
    end do

    call writeMatrix('PMATUUR', 21, 21, ASTER_TRUE, bint_perm)

end subroutine
