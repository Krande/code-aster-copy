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
subroutine te0105(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/cq3d2d.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_FLUN_F'
!                          CAS COQUE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
    integer(kind=8) :: nbres
    parameter(nbres=4)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: pc(3), valpar(nbres), zero, un
    real(kind=8) :: coor2d(18), x, y, z, theta, fplnp1, fpln, fmonp1, fmon
    real(kind=8) :: dfdx(9), dfdy(9), poids, flux, cour, cosa, sina
    real(kind=8) :: matnp(9), fluxp1, coef, long, rp1, rp2, rp3, poid
    integer(kind=8) :: nno, kp, npg1, gi, pi, ivectt, iflux, itemps, nnos
    integer(kind=8) :: ipoids, ivf, idfde, igeom, ndim, jgano, i, k, ier
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    zero = 0.d0
    un = 1.d0
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PFLUXNF', 'L', iflux)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PVECTTR', 'E', ivectt)
!
    theta = zr(itemps+2)
!
    if (nomte .ne. 'THCPSE3 ' .and. nomte .ne. 'THCASE3 ' .and. nomte .ne. 'THCOSE3 ' .and. &
        nomte .ne. 'THCOSE2 ') then
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
!
        call cq3d2d(nno, zr(igeom), 1.d0, zero, coor2d)
!
        do kp = 1, npg1
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, coor2d, &
                        poids, dfdx, dfdy)
            x = zero
            y = zero
            z = zero
            do i = 1, nno
                x = x+zr(igeom+3*(i-1))*zr(ivf+k+i-1)
                y = y+zr(igeom+3*(i-1)+1)*zr(ivf+k+i-1)
                z = z+zr(igeom+3*(i-1)+2)*zr(ivf+k+i-1)
            end do
            valpar(1) = x
            valpar(2) = y
            valpar(3) = z
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(iflux), 4, nompar, valpar, &
                        fmonp1, ier)
            call fointe('FM', zk8(iflux+1), 4, nompar, valpar, &
                        fplnp1, ier)
            valpar(4) = zr(itemps)-zr(itemps+1)
            call fointe('FM', zk8(iflux), 4, nompar, valpar, &
                        fmon, ier)
            call fointe('FM', zk8(iflux+1), 4, nompar, valpar, &
                        fpln, ier)
            if (theta < -0.5) then
                pc(1) = zero
                pc(2) = fmonp1
                pc(3) = fplnp1
            else
                pc(1) = zero
                pc(2) = theta*(fmonp1)+(un-theta)*(fmon)
                pc(3) = theta*(fplnp1)+(un-theta)*(fpln)
            end if
            do gi = 1, nno
                do pi = 1, 3
                    i = 3*(gi-1)+pi-1+ivectt
                    zr(i) = zr(i)+pc(pi)*zr(ivf+k+gi-1)*poids
                end do
            end do
        end do
!
    else if (nomte .eq. 'THCPSE3' .or. nomte .eq. 'THCASE3') &
        then
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'INST'
!
        do kp = 1, npg1
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfde+k), zr(igeom), dfdx, &
                        cour, poids, cosa, sina)
            x = zero
            y = zero
            do i = 1, nno
                x = x+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
                y = y+zr(igeom+2*(i-1)+1)*zr(ivf+k+i-1)
            end do
!
            if (nomte .eq. 'THCASE3') poids = poids*x
!
            valpar(1) = x
            valpar(2) = y
            valpar(3) = zr(itemps)
            call fointe('FM', zk8(iflux), 3, nompar, valpar, &
                        fmonp1, ier)
            call fointe('FM', zk8(iflux+1), 3, nompar, valpar, &
                        fplnp1, ier)
            valpar(3) = zr(itemps)-zr(itemps+1)
            call fointe('FM', zk8(iflux), 3, nompar, valpar, &
                        fmon, ier)
            call fointe('FM', zk8(iflux+1), 3, nompar, valpar, &
                        fpln, ier)
            if (theta < -0.5) then
                pc(1) = zero
                pc(2) = fmonp1
                pc(3) = fplnp1
            else
                pc(1) = zero
                pc(2) = theta*(fmonp1)+(un-theta)*(fmon)
                pc(3) = theta*(fplnp1)+(un-theta)*(fpln)
            end if
            do gi = 1, nno
                do pi = 1, 3
                    i = 3*(gi-1)+pi-1+ivectt
                    zr(i) = zr(i)+pc(pi)*zr(ivf+k+gi-1)*poids
                end do
            end do
        end do
!
    else if (nomte .eq. 'THCOSE3' .or. nomte .eq. 'THCOSE2') &
        then
!
        long = ( &
               zr(igeom+3)-zr(igeom))**2+(zr(igeom+3+1)-zr(igeom+1))**2+(zr(igeom+3+2)-zr(ig&
               &eom+2) &
               )**2
        long = sqrt(long)/2.d0
!       EP  =EP/2.d0
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
!
        rp1 = 1.33333333333333d0
        rp2 = 0.33333333333333d0
        rp3 = 0.33333333333333d0
!
        do kp = 1, npg1
            k = (kp-1)*nno
!
            x = zero
            y = zero
            z = zero
            do i = 1, nno
                x = x+zr(igeom+3*(i-1))*zr(ivf+k+i-1)
                y = y+zr(igeom+3*(i-1)+1)*zr(ivf+k+i-1)
                z = y+zr(igeom+3*(i-1)+2)*zr(ivf+k+i-1)
            end do
!
            valpar(1) = x
            valpar(2) = y
            valpar(3) = z
!
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(iflux), 4, nompar, valpar, &
                        fluxp1, ier)
            valpar(4) = zr(itemps)-zr(itemps+1)
            call fointe('FM', zk8(iflux), 4, nompar, valpar, &
                        flux, ier)
!
!      IMPORTANT: FLUXP1 OU FLUX = FLUX * EPAISSEUR
!
            if (theta < -0.5) then
                coef = fluxp1
            else
                coef = (theta*fluxp1+(un-theta)*flux)/2.d0
            end if
!
            poid = zr(ipoids-1+kp)
!
            matnp(1) = rp1*poid*zr(ivf-1+k+1)
            matnp(2) = rp2*poid*zr(ivf-1+k+1)
            matnp(3) = rp3*poid*zr(ivf-1+k+1)
!
            matnp(4) = rp1*poid*zr(ivf-1+k+2)
            matnp(5) = rp2*poid*zr(ivf-1+k+2)
            matnp(6) = rp3*poid*zr(ivf-1+k+2)
!
            if (nomte .eq. 'THCOSE3') then
                matnp(7) = rp1*poid*zr(ivf-1+k+3)
                matnp(8) = rp2*poid*zr(ivf-1+k+3)
                matnp(9) = rp3*poid*zr(ivf-1+k+3)
            end if
!
            do i = 1, 3*nno
                zr(ivectt-1+i) = zr(ivectt-1+i)+coef*long*matnp(i)
            end do
        end do
!
    end if
!
end subroutine
