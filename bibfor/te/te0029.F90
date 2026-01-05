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
#include "asterfort/dxqfor.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxroep.h"
#include "asterfort/dxtfor.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
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
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, ivectu, imate
    integer(kind=8) :: i, j, ier, jgeom, jpres, itemps, igau
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: rho, epais
    real(kind=8) :: valpar(4), pr
    character(len=8) :: nompar(4)
    character(len=16) :: elas_keyword
    character(len=4) :: fami

    real(c_double) :: cst(5), coor(27), kappa, cdofs_f(12)
    ! NOTICE: see the size of the arrays in the C file: c_interface_plaq_mitc_f

    integer(c_int) :: ncst, ncd, nk, ne0, nwinit
    integer(c_int) :: entities0(1)
    integer(kind=8) :: elas_id

    integer(kind=8), parameter :: size_mat = 30, size_fenicsx = 42*42
    real(c_double) :: F_elem(size_fenicsx), F_elem_int(size_mat), w_0(size_mat)
    real(kind=8) :: signs(size_mat)
    integer(kind=8) :: reorder(size_mat)
    real(c_double) :: pres
    real(kind=8) :: e, nu

! DEB ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
! - Geometry
!
    call jevech('PGEOMER', 'L', jgeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)
!
    do igau = 1, npg
! ----- Get elastic parameters (only isotropic elasticity)
!
        call get_elas_id(zi(imate), elas_id, elas_keyword)
        call get_elas_para(fami, zi(imate), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           e_=e, nu_=nu)
! ----- Fill integration weight vector (Divided by 4 for FEniCS)
!
    end do
!
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
! - Initializations
!
    cst = 0.d0
    coor = 0.d0
    F_elem_int = 0.d0
    F_elem = 0.d0
    w_0 = 0.d0
    nk = 30
!
! - Fill material parameters vector
    kappa = 5.0/6.0
    cst(1) = e
    cst(2) = nu
    cst(3) = kappa
    cst(4) = epais
    cst(5) = pres
! Remplissage du vecteur de coordonnées (3 coordonnées par nœud)
    do i = 0, 8
        coor(3*i+1) = zr(jgeom+3*i)
        coor(3*i+2) = zr(jgeom+3*i+1)
        coor(3*i+3) = zr(jgeom+3*i+2)
    end do
!
    ! Remplissage des degrés de liberté (N1,N2,N4,N3)
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
    nwinit = size(w_0)
    ncd = size(cdofs_f)
    ne0 = 1
    ncst = size(cst)
    entities0(1) = 0
!
    call BP1_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, F_elem)
!
    ![w1, θ_y1, -θ_x1, w2, θ_y2, -θ_x2, w3, θ_y3, -θ_x3, w4, θ_y4, -θ_x4,
    ! θ_y5, -θ_x5, γ_r1, p1, θ_y6, -θ_x6, γ_r2, p2,θ_y7, -θ_x7,
    ! γ_r3, p3, θ_y8, -θ_x8,γ_r4, p4, θ_y9, -θ_x9]
    signs = (/ &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0 &
            /)
    reorder = (/ &
              19, 2, 1, &
              21, 6, 5, &
              22, 8, 7, &
              20, 4, 3, &
              12, 11, &
              24, 28, &
              16, 15, &
              26, 30, &
              14, 13, &
              25, 29, &
              10, 9, &
              23, 27, &
              18, 17 &
              /)

    do i = 1, size_mat
        F_elem_int(i) = -signs(i)*F_elem(reorder(i))
    end do
!
! - Set matrix in output field
    call jevech('PVECTUR', 'E', ivectu)
    do i = 0, nk-1
        zr(ivectu+i) = F_elem_int(i+1)
    end do
999 continue
!
end subroutine
