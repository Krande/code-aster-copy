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
subroutine te0108(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/cq3d2d.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_ECHA_R'
!                          CAS COQUE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
    real(kind=8) :: pc(3), fpl, fmo, coor2d(18), zero, coefmo, coefpl
    real(kind=8) :: dfdx(9), dfdy(9), poids, cour, cosa, sina
    real(kind=8) :: r, matnp(9), long, lamb, text, textmo, textpl, rp1, rp2, rp3
    integer(kind=8) :: i, k, nno, kp, npg1, gi, pi, ivectt, icoef
    integer(kind=8) :: itemps, itex, nnos, ndim, jgano
    integer(kind=8) :: ipoids, ivf, idfde, igeom
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    zero = 0.d0
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTTR', 'E', ivectt)
    call jevech('PCOEFHR', 'L', icoef)
    call jevech('PT_EXTR', 'L', itex)
    call jevech('PINSTR', 'L', itemps)
!
    if (nomte .ne. 'THCOSE3 ' .and. nomte .ne. 'THCOSE2 ') then
        coefmo = zr(icoef)
        coefpl = zr(icoef+1)
        textmo = zr(itex+1)
        textpl = zr(itex+2)
!
        fmo = coefmo*textmo
        fpl = coefpl*textpl
        pc(1) = zero
        pc(2) = fmo
        pc(3) = fpl
    end if
!
    if (nomte .ne. 'THCPSE3 ' .and. nomte .ne. 'THCASE3 ' .and. nomte .ne. 'THCOSE3 ' .and. &
        nomte .ne. 'THCOSE2 ') then
!
        call cq3d2d(nno, zr(igeom), 1.d0, zero, coor2d)
!
        do kp = 1, npg1
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, coor2d, &
                        poids, dfdx, dfdy)
            do gi = 1, nno
                do pi = 1, 3
                    i = 3*(gi-1)+pi-1+ivectt
                    zr(i) = zr(i)+poids*zr(ivf+k+gi-1)*pc(pi)
                end do
            end do
        end do
!
    else if (nomte .eq. 'THCPSE3' .or. nomte .eq. 'THCASE3') &
        then
!
        do kp = 1, npg1
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfde+k), zr(igeom), dfdx, &
                        cour, poids, cosa, sina)
            do pi = 1, 3
!
                if (nomte .eq. 'THCASE3') then
                    r = zero
                    do i = 1, nno
                        r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
                    end do
                    poids = poids*r
                end if
!
            end do
            do gi = 1, nno
                do pi = 1, 3
                    i = 3*(gi-1)+pi-1+ivectt
                    zr(i) = zr(i)+poids*zr(ivf+k+gi-1)*pc(pi)
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
!
!      IMPORTANT: LAMB = CONV * EPAISSEUR
!
        lamb = zr(icoef)/2.d0
        text = zr(itex)
!
        rp1 = 1.33333333333333d0
        rp2 = 0.33333333333333d0
        rp3 = 0.33333333333333d0
!
        do kp = 1, npg1
!
            k = (kp-1)*nno
!
            poids = zr(ipoids-1+kp)
!
            matnp(1) = rp1*poids*zr(ivf-1+k+1)
            matnp(2) = rp2*poids*zr(ivf-1+k+1)
            matnp(3) = rp3*poids*zr(ivf-1+k+1)
!
            matnp(4) = rp1*poids*zr(ivf-1+k+2)
            matnp(5) = rp2*poids*zr(ivf-1+k+2)
            matnp(6) = rp3*poids*zr(ivf-1+k+2)
!
            if (nomte .eq. 'THCOSE3') then
                matnp(7) = rp1*poids*zr(ivf-1+k+3)
                matnp(8) = rp2*poids*zr(ivf-1+k+3)
                matnp(9) = rp3*poids*zr(ivf-1+k+3)
            end if
!
            do i = 1, 3*nno
                zr(ivectt-1+i) = zr(ivectt-1+i)+long*matnp(i)*lamb*text
            end do
        end do
!
    end if
!
end subroutine
