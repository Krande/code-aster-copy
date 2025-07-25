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
subroutine te0553(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/vff2dn.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN_ABSO
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          AMOR_MECA
!          RIGI_MECA
!          FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) ::  fami, poum
    integer(kind=8) :: icodre(5), kpg
    real(kind=8) :: poids, nx, ny, valres(5), e, nu, lambda, mu
    real(kind=8) :: rhocp, rhocs, l0, usl0, depla(6), coef_amor
    real(kind=8) :: rho, taux, tauy, nux, nuy, scal, vnx, vny, vtx, vty
    real(kind=8) :: vituni(2, 2), vect(3, 2, 6), matr(6, 6), jac
    integer(kind=8) :: nno, npg, ipoids, ivf, idfde, igeom, jvDisp
    integer(kind=8) :: ldec, i, l, mater, ndim2
    character(len=8) :: nompar(2)
    integer(kind=8) :: imate, j, ll, ndim
    character(len=16), parameter :: nomres(5) = (/'E        ', 'NU       ', &
                                                  'RHO      ', &
                                                  'COEF_AMOR', 'LONG_CARA'/)
    integer(kind=8) :: jvDispm, jvDispp, jvVect
    integer(kind=8) :: nnos
    aster_logical :: lDamp, lMatr, lVect
    real(kind=8) :: xygau(2)
    integer(kind=8) :: idecpg, idecno
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
    mater = zi(imate)
    fami = 'RIGI'
    poum = '+'
    ASSERT(ndim .ne. 2)
    ndim2 = ndim+1
!
    nompar(1) = 'X'
    nompar(2) = 'Y'
!
    lDamp = option .eq. 'AMOR_MECA'
    lVect = L_VECT(option)
    lMatr = L_MATR(option)
!
!     VITESSE UNITAIRE DANS LES 3 DIRECTIONS
!
    vituni(1, 1) = 1.d0
    vituni(1, 2) = 0.d0
    vituni(2, 1) = 0.d0
    vituni(2, 2) = 1.d0
!
    vect = 0.d0
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kpg = 1, npg
!
! - Get material properties
!
        idecpg = nno*(kpg-1)-1
        ! ----- Coordinates for current Gauss point
        xygau(:) = 0.d0
        do i = 1, nno
            idecno = ndim2*(i-1)-1
            do j = 1, ndim2
                xygau(j) = xygau(j)+zr(ivf+i+idecpg)*zr(igeom+j+idecno)
            end do
        end do
!
        call rcvalb(fami, kpg, 1, poum, mater, &
                    ' ', 'ELAS', 2, nompar, xygau, &
                    4, nomres, valres, icodre, 1)
!       appel LONG_CARA en iarret = 0
        call rcvalb(fami, kpg, 1, poum, mater, &
                    ' ', 'ELAS', 2, nompar, xygau, &
                    1, nomres(5), valres(5), icodre(5), 0)
!
        e = valres(1)
        nu = valres(2)
        rho = valres(3)
        coef_amor = valres(4)
!
        usl0 = 0.d0
        if (icodre(5) .eq. 0) then
            l0 = valres(5)
            usl0 = 1.d0/l0
        end if
        lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
        mu = e/2.d0/(1.d0+nu)
        if (lDamp) then
            rhocp = coef_amor*sqrt((lambda+2.d0*mu)*rho)
            rhocs = coef_amor*sqrt(mu*rho)
        else
            rhocp = (lambda+2.d0*mu)*usl0
            rhocs = mu*usl0
        end if
!
        ldec = (kpg-1)*nno
        call vff2dn(ndim, nno, kpg, ipoids, idfde, &
                    zr(igeom), nx, ny, poids)
        jac = sqrt(nx*nx+ny*ny)
!
!        --- CALCUL DE LA NORMALE UNITAIRE ---
!
        nux = nx/jac
        nuy = ny/jac
!
!        --- CALCUL DE V.N ---
!
        scal = 0.d0
        do i = 1, nno
            do j = 1, 2
                scal = nux*zr(ivf+ldec+i-1)*vituni(j, 1)
                scal = scal+nuy*zr(ivf+ldec+i-1)*vituni(j, 2)
!
!        --- CALCUL DE LA VITESSE NORMALE ET DE LA VITESSE TANGENCIELLE
!
                vnx = nux*scal
                vny = nuy*scal
                vtx = zr(ivf+ldec+i-1)*vituni(j, 1)
                vty = zr(ivf+ldec+i-1)*vituni(j, 2)
                vtx = vtx-vnx
                vty = vty-vny
!
!        --- CALCUL DU VECTEUR CONTRAINTE
!
                taux = rhocp*vnx+rhocs*vtx
                tauy = rhocp*vny+rhocs*vty
!
!        --- CALCUL DU VECTEUR ELEMENTAIRE
!
                do l = 1, nno
                    ll = 2*l-1
                    vect(i, j, ll) = vect(i, j, ll)+taux*zr(ivf+ldec+l-1)*poids
                    vect(i, j, ll+1) = vect(i, j, ll+1)+tauy*zr(ivf+ldec+l-1)*poids
                end do
            end do
        end do
    end do
!
    do i = 1, nno
        do j = 1, 2
            do ldec = 1, 2*nno
                matr(2*(i-1)+j, ldec) = vect(i, j, ldec)
            end do
        end do
    end do
!
! - Get output fields
!
    if (lVect .or. option .Eq. "FORC_NODA") then
        call jevech('PVECTUR', 'E', jvVect)
    end if
!
!       --- PASSAGE AU STOCKAGE TRIANGULAIRE
!
    if (lMatr .or. lDamp) then
        call writeMatrix('PMATUUR', 2*nno, 2*nno, ASTER_TRUE, matr)
    end if
    if (lVect .or. option .Eq. "FORC_NODA") then
        if (option .Eq. "FORC_NODA") then
            call jevech('PDEPLAR', 'L', jvDisp)
            do i = 1, 2*nno
                depla(i) = zr(jvDisp+i-1)
                zr(jvVect+i-1) = 0.d0
            end do
        else
            call jevech('PDEPLMR', 'L', jvDispm)
            call jevech('PDEPLPR', 'L', jvDispp)
            do i = 1, 2*nno
                depla(i) = zr(jvDispm+i-1)+zr(jvDispp+i-1)
                zr(jvVect+i-1) = 0.d0
            end do
        end if
        do i = 1, 2*nno
            do j = 1, 2*nno
                zr(jvVect+i-1) = zr(jvVect+i-1)+matr(i, j)*depla(j)
            end do
        end do
    end if
!
end subroutine
