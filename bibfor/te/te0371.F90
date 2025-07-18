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
subroutine te0371(option, nomte)
!.......................................................................
    implicit none
!
!     BUT: CALCUL DES MATRICES DE RIGIDITE  ELEMENTAIRES EN MECANIQUE
!          ELEMENTS 2D DE COUPLAGE PESANTEUR-SURFACE LIBRE D'UN FLUIDE
!
!          OPTION : 'MASS_MECA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
!
    integer(kind=8) :: icodre(1)
    character(len=8) :: fami, poum
    character(len=16) :: nomte, option
    integer(kind=8) :: igeom, imate, kpg, spt
    integer(kind=8) :: i, j, k, l, ik, ijkl, ldec, kdec, ino, jno
    integer(kind=8) :: ndim, nno, npg2, ipg, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfrde, imatuu
    real(kind=8) :: a(2, 2, 27, 27), rho(1)
    real(kind=8) :: poids, jac
    real(kind=8) :: dxde, dxdk, dyde, dydk
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfrde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PMATUUR', 'E', imatuu)
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'FLUIDE', 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
!
!     INITIALISATION DE LA MATRICE
!
    do k = 1, 2
        do l = 1, 2
            do i = 1, nno
                do j = 1, i
                    a(k, l, i, j) = 0.d0
                end do
            end do
        end do
    end do
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg2
!
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        dxde = 0.d0
        dxdk = 0.d0
        dyde = 0.d0
        dydk = 0.d0
        do i = 1, nno
            dxde = dxde+zr(igeom+3*(i-1))*zr(idfrde+kdec+(i-1)*ndim)
            dxdk = dxdk+zr(igeom+3*(i-1))*zr(idfrde+kdec+(i-1)*ndim+1)
            dyde = dyde+zr(igeom+3*(i-1)+1)*zr(idfrde+kdec+(i-1)*ndim)
            dydk = dydk+zr(igeom+3*(i-1)+1)*zr(idfrde+kdec+(i-1)*ndim+1)
        end do
        jac = dxde*dydk-dxdk*dyde
        poids = abs(jac)*zr(ipoids+ipg-1)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       CALCUL DU TERME RHO * PHI * Z      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
        do ino = 1, nno
            do jno = 1, ino
                a(1, 2, ino, jno) = a(1, 2, ino, jno)+poids*rho(1)*zr(ivf+ldec+ino-1)*zr(ivf+ld&
                                 &ec+jno-1)
            end do
        end do
    end do
!
    do ino = 1, nno
        do jno = 1, ino
            a(2, 1, ino, jno) = a(1, 2, ino, jno)
        end do
    end do
!
! PASSAGE DU STOCKAGE RECTANGULAIRE AU STOCKAGE TRIANGULAIRE
!
    ijkl = 0
    ik = 0
    do k = 1, 2
        do l = 1, 2
            do i = 1, nno
                ik = ((2*i+k-3)*(2*i+k-2))/2
                do j = 1, i
                    ijkl = ik+2*(j-1)+l
                    zr(imatuu+ijkl-1) = a(k, l, i, j)
                end do
            end do
        end do
    end do
!
end subroutine
