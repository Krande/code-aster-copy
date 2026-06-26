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

subroutine te0029(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module, only: behaviourOption
!
    use c_interface_plaq_mitc_f
    use iso_c_binding

    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dxroep.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/utpvgl.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
!
!
    character(len=16), intent(in) :: option, nomte

!     IN  OPTION : NOM DE L'OPTION A CALCULER
!     IN  NOMTE  : NOM DU TYPE_ELEMENT
!     -----------------------------------------------------------------
!     CALCUL DE PRESSION SUR LES ELEMENTS PLAQ_MITC
!         OPTIONS TRAITEES   ==> CHAR_MECA_PRES_R
!     -----------------------------------------------------------------
    integer(kind=8) :: ndim, nno, nnos, npg, jpoids, jvf, jdfde, jgano, jgeom, jmate
    integer(kind=8) :: i, j, ier, jpres, itemps, ivectu, elas_id
    real(c_double) :: pgl(3, 3), xyzl(3, 4), valpar(4)
    real(c_double) :: e, nu, epais, rho, pres, pr
    character(len=8) :: fami
    character(len=8) :: nompar(4)
    character(len=16) :: elas_keyword

    real(c_double) :: cst(5), coor(27), cdofs_f(12)
    ! NOTICE: see the size of the arrays in the C file: c_interface_plaq_mitc_f

    integer(c_int) :: ncst, ncd, ne0, ne1, ne2, ne3
    integer(c_int) :: entities0(1), entities1(1), entities2(1), entities3(1)

    integer(kind=4), parameter :: size_init = 30, size_final = 22
    real(c_double), dimension(size_init) :: F_elem, F0, F1, F2, F3, F4, w_0
    real(c_double) :: signs(size_final)
    integer(kind=8) :: reorder(size_final)

! --------------------------------------------------------------------
! - Finite element informations
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=jpoids, jvf=jvf, jdfde=jdfde, jgano=jgano)
! - Geometry
!
    call jevech('PGEOMER', 'L', jgeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', jmate)
!
! ----- Get elastic parameters (only isotropic elasticity)
! FIXME: for instance E and NU are supposed to be equal for all quadpoints
!
    call get_elas_id(zi(jmate), elas_id, elas_keyword)
    call get_elas_para(fami, zi(jmate), '+', 1, 1, &
                       elas_id, elas_keyword, &
                       e_=e, nu_=nu)

    call dxroep(rho, epais)
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
!              ------------------------------
        call jevech('PPRESSR', 'L', jpres)
        call dxtpgl(zr(jgeom), pgl)
        call utpvgl(nno, 3, pgl, zr(jgeom), xyzl)
        pres = zr(jpres)
!
! --- CAS DES CHARGEMENTS DE FORME FONCTION
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then

        call jevech('PPRESSF', 'L', jpres)
        if (zk8(jpres) .eq. '&FOZERO') goto 999
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        pres = 0.d0
        do j = 0, nno-1
            valpar(1) = zr(jgeom+3*j)
            valpar(2) = zr(jgeom+3*j+1)
            valpar(3) = zr(jgeom+3*j+2)
            call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                        pr, ier)
            pres = pres+pr
        end do
        pres = pres/nno

    end if
!
!
! Initializations
!
    w_0 = 0.d0
!
! Fill material parameters vector
    cst(1) = e
    cst(2) = nu
    cst(3) = 5.0/6.0
    cst(4) = epais
    cst(5) = pres
    ncst = size(cst)
! Remplissage du vecteur de coordonnées (3 coordonnées par nœud)
    do i = 0, 8
        coor(3*i+1) = zr(jgeom+3*i)
        coor(3*i+2) = zr(jgeom+3*i+1)
        coor(3*i+3) = zr(jgeom+3*i+2)
    end do
! Remplissage des coordonées (N1, N4, N2, N3)
    cdofs_f(1) = coor(1)
    cdofs_f(2) = coor(2)
    cdofs_f(3) = coor(3)
!
    cdofs_f(4) = coor(10)
    cdofs_f(5) = coor(11)
    cdofs_f(6) = coor(12)
!
    cdofs_f(7) = coor(4)
    cdofs_f(8) = coor(5)
    cdofs_f(9) = coor(6)
!
    cdofs_f(10) = coor(7)
    cdofs_f(11) = coor(8)
    cdofs_f(12) = coor(9)
!
    ncd = size(cdofs_f)
!
! On définit 4 entities car 4 arêtes dans le quadrangle
    entities0(1) = 0
    entities1(1) = 1
    entities2(1) = 2
    entities3(1) = 3
    ne0 = 1
    ne1 = 1
    ne2 = 1
    ne3 = 1
!
    call BP1_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities0, ne0, cst, ncst, F0)
    call BP2_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities0, ne0, cst, ncst, F1)
    call BP2_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities1, ne1, cst, ncst, F2)
    call BP2_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities2, ne2, cst, ncst, F3)
    call BP2_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities3, ne3, cst, ncst, F4)

! Remplissage du vecteur
    do i = 1, size_init
        F_elem(i) = F0(i)+F1(i)+F2(i)+F3(i)+F4(i)
    end do
!
! Reorganisation du vecteur
    ![w1, θ_y1, -θ_x1, w2, θ_y2, -θ_x2, w3, θ_y3, -θ_x3, w4, θ_y4, -θ_x4,
    ! θ_y5, -θ_x5, θ_y6, -θ_x6, θ_y7, -θ_x7, θ_y8, -θ_x8, θ_y9, -θ_x9]
    !
    reorder = (/ &
              19, 2, 1, &
              21, 6, 5, &
              22, 8, 7, &
              20, 4, 3, &
              12, 11, &
              16, 15, &
              14, 13, &
              10, 9, &
              18, 17 &
              /)
    signs = (/ &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0 &
            /)
!
! - Set matrix in output field
    call jevech('PVECTUR', 'E', ivectu)
    do i = 1, size_final
        zr(ivectu+i-1) = -signs(i)*F_elem(reorder(i))
    end do
999 continue
!
end subroutine
