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
subroutine projsg(x3dca, x3d1, x3d2, normal, x3dp, &
                  xbar, iproj, excent)
    implicit none
!  DESCRIPTION : PROJECTION DU NOEUD CABLE X3DCA(3) SUR UN SEGMENT
!  -----------   DE NOEUDS EXTREMITES X3D1(3) ET X3D2(3)
!                DEUX ETAPES :
!                  - PROJECTION DU NOEUD X3DCA(3) SUR LA DROITE PASSANT
!                    PAR LES NOEUDS X3D1(3) ET X3D2(3)
!                  - TEST D'APPARTENANCE DU POINT PROJETE X3DP(3) AU
!                    SEGMENT D'EXTREMITES X3D1(3) ET X3D2(3), PAR
!                    CALCUL DES COORDONNEES BARYCENTRIQUES
!
!                APPELANT : PROJKB
!
!  IN     : X3DCA  : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU NOEUD CABLE CONSIDERE
!  IN     : X3D1   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU PREMIER NOEUD EXTREMITE DU SEGMENT
!  IN     : X3D2   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU SECOND NOEUD EXTREMITE DU SEGMENT
!  OUT    : NORMAL : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DANS LE REPERE GLOBAL DU VECTEUR NORMAL
!                    A LA DROITE SUR LAQUELLE LE NOEUD CABLE EST PROJETE
!  OUT    : X3DP   : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU POINT PROJETE
!  OUT    : XBAR   : REAL*8 , VECTEUR DE DIMENSION 3
!                    SI PROJECTION REUSSIE : COORDONNEES BARYCENTRIQUES
!                    DU POINT PROJETE (BARYCENTRE DES EXTREMITES DU
!                    SEGMENT)
!  OUT    : IPROJ  : INTEGER , SCALAIRE
!                    INDICE DE PROJECTION
!                    IPROJ = -1  PROJECTION NON REUSSIE
!                    IPROJ =  0  LE POINT PROJETE EST A L'INTERIEUR
!                                DU SEGMENT
!                    IPROJ =  2  LE POINT PROJETE COINCIDE AVEC UN DES
!                                NOEUDS EXTREMITES
!  OUT    : EXCENT : REAL*8 , SCALAIRE
!                    DISTANCE DU NOEUD CABLE A LA DROITE SUR LAQUELLE
!                    IL EST PROJETE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterc/matfpe.h"
#include "asterc/r8prem.h"
#include "asterfort/r8inir.h"
#include "asterfort/tstbar.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
    real(kind=8) :: x3dca(*), x3d1(*), x3d2(*), normal(*), x3dp(*), xbar(*)
    real(kind=8) :: excent
    integer(kind=8) :: iproj
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: alpha1, alpha2, beta1, beta2, dx12, dx3d(3), dy12, dz12
    real(kind=8) :: epsg, n1n1, n1n2, n2n2, nrm2, plan1(4), plan2(4)
    blas_int :: b_incx, b_incy, b_n
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call matfpe(-1)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   PROJECTION DU NOEUD CABLE SUR LA DROITE PASSANT PAR LES DEUX
!     NOEUDS EXTREMITES DU SEGMENT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    epsg = 1.0d+08*r8prem()
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    nrm2 = dble(max(dnrm2(b_n, x3d1(1), b_incx), dnrm2(b_n, x3d2(1), b_incx)))
    if (nrm2 .eq. 0.0d0) then
        iproj = -1
        goto 9999
    end if
!
! 1.1 DETERMINATION DES DEUX PLANS DONT L'INTERSECTION DEFINIT LA DROITE
! --- PASSANT PAR LES DEUX NOEUDS EXTREMITES DU SEGMENT
!
    dx12 = x3d2(1)-x3d1(1)
    dy12 = x3d2(2)-x3d1(2)
    dz12 = x3d2(3)-x3d1(3)
!
!.... CAS OU LE PREMIER PLAN A POUR EQUATION X = CSTE
!
    if (dble(abs(dx12))/nrm2 .lt. epsg) then
!
        plan1(1) = 1.0d0
        plan1(2) = 0.0d0
        plan1(3) = 0.0d0
        plan1(4) = -x3d1(1)
!
!....... CAS OU LE SECOND PLAN A POUR EQUATION Y = CSTE
!
        if (dble(abs(dy12))/nrm2 .lt. epsg) then
!
            plan2(1) = 0.0d0
            plan2(2) = 1.0d0
            plan2(3) = 0.0d0
            plan2(4) = -x3d1(2)
!
!....... CAS OU LE SECOND PLAN A POUR EQUATION Z = CSTE
!
        else if (dble(abs(dz12))/nrm2 .lt. epsg) then
!
            plan2(1) = 0.0d0
            plan2(2) = 0.0d0
            plan2(3) = 1.0d0
            plan2(4) = -x3d1(3)
!
!....... AUTRES CAS
!
        else
!
            plan2(1) = 0.0d0
            plan2(2) = dz12
            plan2(3) = -dy12
            plan2(4) = -x3d1(2)*dz12+x3d1(3)*dy12
!
        end if
!
!.... CAS OU LE PREMIER PLAN A POUR EQUATION Y = CSTE
!
    else if (dble(abs(dy12))/nrm2 .lt. epsg) then
!
        plan1(1) = 0.0d0
        plan1(2) = 1.0d0
        plan1(3) = 0.0d0
        plan1(4) = -x3d1(2)
!
!....... CAS OU LE SECOND PLAN A POUR EQUATION Z = CSTE
!
        if (dble(abs(dz12))/nrm2 .lt. epsg) then
!
            plan2(1) = 0.0d0
            plan2(2) = 0.0d0
            plan2(3) = 1.0d0
            plan2(4) = -x3d1(3)
!
!....... AUTRES CAS
!
        else
!
            plan2(1) = -dz12
            plan2(2) = 0.0d0
            plan2(3) = dx12
            plan2(4) = x3d1(1)*dz12-x3d1(3)*dx12
!
        end if
!
!.... CAS OU LE PREMIER PLAN A POUR EQUATION Z = CSTE
!
    else if (dble(abs(dz12))/nrm2 .lt. epsg) then
!
        plan1(1) = 0.0d0
        plan1(2) = 0.0d0
        plan1(3) = 1.0d0
        plan1(4) = -x3d1(3)
!
        plan2(1) = dy12
        plan2(2) = -dx12
        plan2(3) = 0.0d0
        plan2(4) = -x3d1(1)*dy12+x3d1(2)*dx12
!
!.... CAS GENERAL
!
    else
!
        plan1(1) = 0.0d0
        plan1(2) = dz12
        plan1(3) = -dy12
        plan1(4) = -x3d1(2)*dz12+x3d1(3)*dy12
!
        plan2(1) = -dz12
        plan2(2) = 0.0d0
        plan2(3) = dx12
        plan2(4) = x3d1(1)*dz12-x3d1(3)*dx12
!
    end if
!
! 1.2 EXCENTRICITE ET COORDONNEES DU POINT PROJETE
! ---
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    n1n1 = ddot(b_n, plan1(1), b_incx, plan1(1), b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    n1n2 = ddot(b_n, plan1(1), b_incx, plan2(1), b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    n2n2 = ddot(b_n, plan2(1), b_incx, plan2(1), b_incy)
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    alpha1 = ddot(b_n, plan1(1), b_incx, x3dca(1), b_incy)+plan1(4)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    alpha2 = ddot(b_n, plan2(1), b_incx, x3dca(1), b_incy)+plan2(4)
!
    beta1 = -n2n2*alpha1+n1n2*alpha2
    beta2 = n1n2*alpha1-n1n1*alpha2
!
    excent = ( &
             plan1(1)*beta1+plan2(1)*beta2)*(plan1(1)*beta1+plan2(1)*beta2)+(plan1(2)*beta1+plan2&
             &(2)*beta2)*(plan1(2)*beta1+plan2(2)*beta2)+(plan1(3)*beta1+plan2(3)*beta2)*(plan1(3&
             &)*beta1+plan2(3)*beta2 &
             )
    excent = dble(sqrt(excent))/(n1n1*n2n2-n1n2*n1n2)
    dx3d(1) = x3dca(1)-x3d1(1)
    dx3d(2) = x3dca(2)-x3d1(2)
    dx3d(3) = x3dca(3)-x3d1(3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    nrm2 = dnrm2(b_n, dx3d(1), b_incx)
    dx3d(1) = x3dca(1)-x3d2(1)
    dx3d(2) = x3dca(2)-x3d2(2)
    dx3d(3) = x3dca(3)-x3d2(3)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    nrm2 = dble(max(nrm2, dnrm2(b_n, dx3d(1), b_incx)))
    if (nrm2 .eq. 0.0d0) then
        iproj = -1
        goto 9999
    end if
    if (dble(abs(excent))/nrm2 .lt. epsg) excent = 0.0d0
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, x3dca(1), b_incx, x3dp(1), b_incy)
    if (excent .gt. 0.0d0) then
        normal(1) = (plan1(1)*beta1+plan2(1)*beta2)/(n1n1*n2n2-n1n2*n1n2)
        normal(2) = (plan1(2)*beta1+plan2(2)*beta2)/(n1n1*n2n2-n1n2*n1n2)
        normal(3) = (plan1(3)*beta1+plan2(3)*beta2)/(n1n1*n2n2-n1n2*n1n2)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.0d0, normal(1), b_incx, x3dp(1), &
                   b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, -1.0d0/excent, normal(1), b_incx)
    else
        call r8inir(3, 0.0d0, normal(1), 1)
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   TEST D'APPARTENANCE DU POINT PROJETE AU SEGMENT, PAR DETERMINATION
!     DES COORDONNEES BARYCENTRIQUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    call tstbar(2, x3d1(1), x3d2(1), x3d2(1), x3d2(1), &
                x3dp(1), xbar(1), iproj)
!
9999 continue
!
! --- FIN DE PROJSG.
    call matfpe(1)
!
end subroutine
