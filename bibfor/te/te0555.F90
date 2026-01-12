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
    character(len=4) :: fami
    integer(kind=8) :: ndim, nno, igeom, imate, i, j
    integer(kind=8), parameter :: size_init = 21, size_final = 15, size_fenicsx = 21*21

    real(c_double) :: cst(4), coor(18), kappa, cdofs_f(9)
    ! NOTICE: see the size of the arrays in the C file: c_interface_tria_mitc_j
    real(c_double), dimension(size_fenicsx) :: A0, A1, A2, A3, A_int

    integer(c_int) :: ncst, ncd, ne0, ne1, ne2, nwinit
    integer(c_int) :: entities0(1), entities1(1), entities2(1)
!
    integer(kind=8) :: reorder(size_final)
    real(kind=8) :: signs(size_final)
    real(c_double) :: w_0(size_init)
    real(kind=8), dimension(size_init, size_init) :: bint
    real(kind=8), dimension(size_final, size_final) :: bf
    real(kind=8) :: e, nu, epais, rho
    integer(kind=8) :: elas_id, igau, ipoids, npg1
    integer(kind=8) :: nnos, ivf, idfde
    character(len=16) :: elas_keyword

    real(kind=8) :: AA(15, 15), BB(3, 3), CC(15, 3), DD(3, 3), CDinv(15, 3)
    integer(kind=8) :: zz_order(15), gamma_order(3), p_order(3)
    real(kind=8) :: AAcondensed(size_final, size_final)

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
    cst = 0.d0
    coor = 0.d0
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

! Remplissage des coordonées (pas de permut ici)
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
    ncst = size(cst)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MATRICE DE RIGIDITÉ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel de la fonction C++ - Matrice de rigidité -> Part1
    call BP4_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A0)

! Appel de la fonction C++ - Matrice de rigidité -> Part2
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A1)
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities1, ne1, cst, ncst, A2)
    call BP5_tr6_Fortran(w_0, nwinit, cdofs_f, ncd, entities2, ne2, cst, ncst, A3)

! Remplissage de la matrice intermédiaire
    do i = 1, size_fenicsx
        A_int(i) = A0(i)+A1(i)+A2(i)+A3(i)
    end do

! Remplissage de la matrice K à partir de A_int
    do i = 1, size_init
        do j = 1, size_init
            bint(i, j) = A_int((j-1)*size_init+i)
        end do
    end do

! Recupérer les bloques
    ! (z, z) où z = (thetha_x theta_y w)
    AA = bint(1:15, 1:15)

    ! (gamma, gamma)
    BB = bint(16:18, 16:18)

    ! (z, p) où z = (thetha_x theta_y w)
    CC = bint(1:15, 19:21)

    ! (gamma, p)
    DD = bint(16:18, 19:21)

    do j = 1, 3
        CDinv(:, j) = CC(:, j)/DD(j, j)
    end do

    ! Condensation statique
    AAcondensed = AA+matmul(matmul(CDinv, BB), transpose(CDinv))

! Reorganisation du vecteur
    ![w1, θ_y1, -θ_x1, w2, θ_y2, -θ_x2, w3, θ_y3, -θ_x3,
    ! θ_y6, -θ_x6, θ_y4, -θ_x4, θ_y5, -θ_x5]
    !
    reorder = (/ &
              13, 2, 1, &
              14, 4, 3, &
              15, 6, 5, &
              12, 11, &
              8, 7, &
              10, 9 &
              /)
    signs = (/ &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0 &
            /)
!
    do i = 1, size_final
        do j = 1, size_final
            bf(i, j) = signs(i)*signs(j)*AAcondensed(reorder(i), reorder(j))
        end do
    end do

    call writeMatrix('PMATUUR', size_final, size_final, ASTER_TRUE, bf)

end subroutine
