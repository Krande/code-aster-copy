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
subroutine lcoptg(nmat, mater, nr, nvi, drdy, &
                  sigeps, dsde, iret)
!     CALCUL DU JACOBIEN DU SYSTEME NL A RESOUDRE = DRDY(DY)
!     POUR LE MODELE BETON_BURGER
!     IN  NR     :  DIMENSION JACOBIEN
!         NMAT   :  DIMENSION MATER
!         MATER  :  COEFFICIENTS MATERIAU
!         NR     :  DIMENSION MATRICE JACOBIENNE
!         DRDY   :  MATRICE JACOBIENNE
!         SIGEPS :  =1 si DY=(DSIG, DZ) et =0 si DY=(DEPSE, DZ)
!     OUT DSDE   :  MATRICE TANGENTE EN VITESSE
!     ----------------------------------------------------------------
! aslint: disable=W1306
    implicit none
!     ----------------------------------------------------------------
#include "asterc/r8prem.h"
#include "asterfort/lcopli.h"
#include "asterfort/mgauss.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
    integer(kind=8) :: nmat, nr, iret, ndt, ndi, i, j, k, nvi, norm, sigeps
    real(kind=8) :: hook(6, 6), drdy(nr, nr), dsde(6, 6), mater(nmat, 2)
    real(kind=8) :: y0(ndt, ndt), y1(ndt, nvi), y2(nvi, ndt), y3(nvi, nvi)
    real(kind=8) :: y4(ndt, ndt), y5(ndt, ndt), det, maxi, mini
    real(kind=8) :: dsdeb(ndt, ndt)
    character(len=4) :: cargau
! === =================================================================
! --- RECHERCHE DU MAXIMUM DE DRDY
! === =================================================================
!
    norm = 0
    if (norm .eq. 0) goto 30
!
    maxi = 0.d0
    do i = 1, nr
        do j = 1, nr
            if (abs(drdy(i, j)) .gt. maxi) maxi = abs(drdy(i, j))
        end do
    end do
!
! === =================================================================
! --- DIMENSIONNEMENT A R8PREM
! === =================================================================
    mini = r8prem()*maxi
    do i = 1, nr
        do j = 1, nr
            if (abs(drdy(i, j)) .lt. mini) drdy(i, j) = 0.d0
        end do
    end do
!
30  continue
!
! === =================================================================
! --- SEPARATION DES TERMES DU JACOBIEN
! === =================================================================
    do i = 1, ndt
        do j = 1, ndt
            y0(i, j) = drdy(i, j)
        end do
    end do
    do i = 1, ndt
        do j = 1, nvi
            y1(i, j) = drdy(i, j+ndt)
        end do
    end do
    do i = 1, nvi
        do j = 1, ndt
            y2(i, j) = drdy(i+ndt, j)
        end do
    end do
    do i = 1, nvi
        do j = 1, nvi
            y3(i, j) = drdy(i+ndt, j+ndt)
        end do
    end do
!
! === =================================================================
! --- CONSTRUCTION TENSEUR RIGIDITE ELASTIQUE A T+DT
! === =================================================================
    call lcopli('ISOTROPE', '3D      ', mater, hook)
!
!     CHOIX DES PARAMETRES DE LANCEMENT DE MGAUSS
!     METHODE 'S' : SURE
    cargau = 'NCSP'
!     METHODE 'W' : RATEAU
!      CARGAU = 'NCWP'
! === =================================================================
! --- CONSTRUCTION TENSEUR CONSTITUTIF TANGENT DSDE
! === =================================================================
!     Y2=INVERSE(Y3)*Y2
    call mgauss(cargau, y3, y2, nvi, nvi, &
                ndt, det, iret)
    if (iret .gt. 1) then
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)
        goto 999
    end if
! --- PRODUIT DU TERME (Y3)^-1 * Y2 = Y4
    call promat(y1, ndt, ndt, nvi, y2, &
                nvi, nvi, ndt, y4)
!
! --- DIFFERENCE DE MATRICE (DR1DY1 - Y4) = Y5
    do i = 1, ndt
        do j = 1, ndt
            y5(i, j) = y0(i, j)-y4(i, j)
        end do
    end do
!
! --- INVERSION DU TERME Y5
    call r8inir(ndt*ndt, 0.d0, dsdeb, 1)
    do i = 1, ndt
        dsdeb(i, i) = 1.d0
    end do
    call mgauss(cargau, y5, dsdeb, ndt, ndt, &
                ndt, det, iret)
!
    if (iret .gt. 1) then
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)
    else
        call r8inir(36, 0.d0, dsde, 1)
        if (sigeps .eq. 1) then
            do i = 1, ndt
                do j = 1, ndt
                    do k = 1, ndt
                        dsde(i, j) = dsdeb(i, j)
                    end do
                end do
            end do
        else
            do i = 1, ndt
                do j = 1, ndt
                    do k = 1, ndt
                        dsde(i, j) = dsde(i, j)+hook(i, k)*dsdeb(k, j)
                    end do
                end do
            end do
        end if
!
    end if
!
999 continue
end subroutine
