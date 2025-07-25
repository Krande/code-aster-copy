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
subroutine te0225(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          COQUE 1D
!                          OPTION : 'CHAR_MECA_TEMP_R'
!                          ELEMENT: MECXSE3
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: i, ip, k, kp, igeom, icaco, ivectt, imate
    integer(kind=8) :: ivf, idfdk, nno, npg, jcoopg, j, nbres, jdfd2
!-----------------------------------------------------------------------
    integer(kind=8) :: ipoids, iret1, iret2, iret3, iret4, jgano, ndim
    integer(kind=8) :: nnos
    real(kind=8) :: tref
!-----------------------------------------------------------------------
    parameter(nbres=3)
    character(len=32) :: phenom
    character(len=16) :: nomres(nbres)
    character(len=4) :: fami
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: valres(nbres), dfdx(3), r, cour, jac, cosa, sina
    real(kind=8) :: tpg1, tpg2, tpg3, tpg, zero, un, deux, x3
    real(kind=8) :: h, epsthe, nu, coef, axis
!
    data zero, un, deux/0.d0, 1.d0, 2.d0/
!     ------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jdfd2=jdfd2, &
                     jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icaco)
    call jevech('PMATERC', 'L', imate)
    call jevech('PVECTUR', 'E', ivectt)
!     TEMPERATURE DE REFERENCE
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, &
                1, tref, iret1)
!
! --- RECUPERATION DE LA NATURE DU MATERIAU DANS PHENOM
!     -------------------------------------------------
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
!
    if (phenom .eq. 'ELAS') then
!
! ==== CALCUL ISOTROPE HOMOGENE =====
!
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nomres(3) = 'ALPHA'
!
        h = zr(icaco)
        axis = zero
        if (nomte .eq. 'MECXSE3 ') axis = un
!
!     ** BOUCLE CONCERNANT LES POINTS DE GAUSS **************
!
        do kp = 1, npg
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                        cour, jac, cosa, sina)
            r = zero
            tpg = zero
            call rcvarc(' ', 'TEMP', '+', fami, kp, &
                        1, tpg2, iret2)
            call rcvarc(' ', 'TEMP', '+', fami, kp, &
                        2, tpg1, iret3)
            call rcvarc(' ', 'TEMP', '+', fami, kp, &
                        3, tpg3, iret4)
            do i = 1, nno
                r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
            end do
            if (nomte .eq. 'MECXSE3 ') jac = jac*r
!
!---- UTILISATION DE 4 POINTS DE GAUSS DANS L'EPAISSEUR
!---- COMME POUR LA LONGUEUR
!
            do ip = 1, npg
                x3 = zr(jcoopg+ip-1)
                tpg = tpg1*(un-x3**2)+x3*(tpg3*(un+x3)-tpg2*(un-x3))/deux
                call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                            ' ', 'ELAS', 1, 'TEMP', [tpg], &
                            2, nomres, valres, icodre, 1)
                call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                            ' ', 'ELAS', 1, 'TEMP', [tpg], &
                            1, nomres(3), valres(3), icodre(3), 0)
                if (((iret1+iret2+iret3+iret4) .ge. 1) .and. (icodre(3) .eq. 0)) then
                    call utmess('F', 'CALCULEL_15')
                else if (icodre(3) .ne. 0) then
                    epsthe = 0.d0
                else
                    epsthe = (tpg-tref)*valres(3)
                end if
                nu = valres(2)
                coef = valres(1)*jac*epsthe*zr(ipoids+ip-1)*(h/deux)
                coef = coef/(un-nu)
!
                do i = 1, nno
                    j = 3*(i-1)
                    zr(ivectt+j) = zr(ivectt+j)+coef*(axis*zr(ivf+k+i-1)/r-dfdx(i)*sina)
                    zr(ivectt+j+1) = zr(ivectt+j+1)+coef*dfdx(i)*cosa
                    zr(ivectt+j+2) = zr(ivectt+j+2)-coef*x3*h/deux*(axis*zr(ivf+k+i-1)*sina/r&
                                     &-dfdx(i))
                end do
            end do
        end do
    else
!  ==== CALCUL ANISOTROPE  =====
        call utmess('F', 'ELEMENTS3_49')
    end if
end subroutine
