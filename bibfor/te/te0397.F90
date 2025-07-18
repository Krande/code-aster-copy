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

subroutine te0397(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          COQUE 1D
!                          OPTION : 'CHAR_MECA_PRES_R  '
!                          ELEMENT: MECXSE3
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: nno, nnos, jgano, ndim, nddl, kp, npg, ipoids, ivf, idfdk, igeom
    integer(kind=8) :: ivectu, k, i, l, ipres, ier, iadzi, iazk24, itemps
    real(kind=8) :: valpar(4), poids, r, fx, fy, f3, nx, ny, cour, dfdx(3), pr
    character(len=8) :: nompar(4), nomail, elrefe
    character(len=24) :: valk
! DEB ------------------------------------------------------------------
!
    call elref1(elrefe)
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
    nddl = 3
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
!          ------------------------------
        call jevech('PPRESSR', 'L', ipres)
!
        do kp = 1, npg
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                        cour, poids, nx, ny)
            r = 0.d0
            fx = 0.d0
            fy = 0.d0
            do i = 1, nno
                l = (kp-1)*nno+i
!-----------------------------------------------------
!              LE SIGNE MOINS CORRESPOND A LA CONVENTION :
!                 UNE PRESSION POSITIVE PROVOQUE UN GONFLEMENT
!-----------------------------------------------------
                f3 = -zr(ipres+i-1)
                fx = fx+nx*f3*zr(ivf+l-1)
                fy = fy+ny*f3*zr(ivf+l-1)
                r = r+zr(igeom+2*(i-1))*zr(ivf+l-1)
            end do
            poids = poids*r
            do i = 1, nno
                l = (kp-1)*nno+i
                zr(ivectu+nddl*(i-1)) = zr(ivectu+nddl*(i-1))+fx*zr(ivf+l-1)*poids
                zr(ivectu+nddl*(i-1)+1) = zr(ivectu+nddl*(i-1)+1)+fy*zr(ivf+l-1)*poids
            end do
        end do
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
!              ------------------------------
        call jevech('PPRESSF', 'L', ipres)
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        do i = 0, nno-1
            valpar(1) = zr(igeom+3*i)
            valpar(2) = zr(igeom+3*i+1)
            valpar(3) = zr(igeom+3*i+2)
            call fointe('FM', zk8(ipres), 4, nompar, valpar, &
                        pr, ier)
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk = nomail
                call utmess('F', 'ELEMENTS4_92', sk=valk)
            end if
        end do
!
    end if
!
end subroutine
