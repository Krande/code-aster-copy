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
subroutine dxiner(nnoe, xyzg1, rho, epais, mass, &
                  cdg, inerti)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
    integer(kind=8) :: nnoe
    real(kind=8) :: xyzg1(3, *), rho, epais, mass, cdg(*), inerti(*)
!     CALCULE LE CDG ET LA MASSE D'UNE MAILLE TRIA ET QUAD
!
!     ------------------------------------------------------------------
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), zero
    real(kind=8) :: pgl(3, 3), xyzl1(3, 4)
    real(kind=8) :: xyzg(3, 8), xyzl(3, 8)
    real(kind=8) :: igxx, igyy, igxy, matine(6), igzz
    real(kind=8) :: inert0(6)
!
! --- INITIALISATIONS :
!     ---------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idec, idfdx, idfdy, ino, ipg, ipoids
    integer(kind=8) :: ivf, j, jdec, jgano, jno, k, kdec
    integer(kind=8) :: ldec, ndim, nno, nnos, npg1
    real(kind=8) :: aire, axg, axggau, axl, axlgau, axx, axxgau
    real(kind=8) :: axy, axygau, ayg, ayggau, ayl, aylgau, ayy
    real(kind=8) :: ayygau, azg, azggau, douze, roep, s1
    real(kind=8) :: sigau, un, undemi, xg, xgau, xl, yg
    real(kind=8) :: ygau, yl, zg
!-----------------------------------------------------------------------
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    douze = 12.0d0
!
    aire = zero
    axg = zero
    ayg = zero
    azg = zero
    axl = zero
    ayl = zero
    axx = zero
    ayy = zero
    axy = zero
!
! --- RECUPERATION DES DONNEES RELATIVES A L'INTEGRATION DES ELEMENTS
! --- DE TYPE 'FACE6' ET 'FACE8' :
!     -------------------------
    call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    roep = rho*epais
!
! --- DETERMINATION DE LA MATRICE DE PASSAGE DU REPERE GLOBAL
! --- AU REPERE LOCAL :
!     ---------------
    if (nnoe .eq. 3) then
        call dxtpgl(xyzg1, pgl)
    else if (nnoe .eq. 4) then
        call dxqpgl(xyzg1, pgl)
    end if
!
! --- DETERMINATION DES COORDONNEES DES NOEUDS DANS LE REPERE LOCAL :
!     -------------------------------------------------------------
    call utpvgl(nnoe, 3, pgl, xyzg1, xyzl1)
!
! --- AFFECTATION DES COORDONNEES DES NOEUDS DE L'ELEMENT FACE6
! --- OU FACE8 CORRESPONDANT A L'ELEMENT DE PLAQUE COURANT,
! --- XYZG DESIGNENT LES COORDONNNEES DES CONNECTIVITES DANS
! --- LE REPERE GLOBAL, XYZL DESIGNENT LES COORDONNEES DES
! --- CONNECTIVITES DANS LE REPERE LOCAL A L'ELEMENT DE PLAQUE:
!     -------------------------------------------------------
!
! --- NOEUDS SOMMETS :
!     --------------
    do ino = 1, nnoe
        do k = 1, 3
            xyzl(k, ino) = xyzl1(k, ino)
            xyzg(k, ino) = xyzg1(k, ino)
        end do
    end do
!
! --- NOEUDS MILIEUX DES COTES , ON PREND COMME COORDONNEES
! --- LA DEMI-SOMME DES COORDONNEES DES NOEUDS SOMMETS PUISQU'IL
! --- S'AGIT D'ELEMENTS DE PLAQUE :
!     ---------------------------
    do ino = 1, nnoe-1
        do k = 1, 3
            xyzl(k, nnoe+ino) = undemi*(xyzl1(k, ino)+xyzl1(k, ino+1))
            xyzg(k, nnoe+ino) = undemi*(xyzg1(k, ino)+xyzg1(k, ino+1))
        end do
    end do
!
    do k = 1, 3
        xyzl(k, nnoe+nnoe) = undemi*(xyzl1(k, 1)+xyzl1(k, nnoe))
        xyzg(k, nnoe+nnoe) = undemi*(xyzg1(k, 1)+xyzg1(k, nnoe))
    end do
!
! --- CALCUL DES PRODUITS VECTORIELS OMI X OMJ :
!     ----------------------------------------
    do ino = 1, nno
        do jno = 1, nno
            sx(ino, jno) = xyzg(2, ino)*xyzg(3, jno)-xyzg(3, ino)*xyzg(2, jno)
            sy(ino, jno) = xyzg(3, ino)*xyzg(1, jno)-xyzg(1, ino)*xyzg(3, jno)
            sz(ino, jno) = xyzg(1, ino)*xyzg(2, jno)-xyzg(2, ino)*xyzg(1, jno)
        end do
    end do
!
! --- BOUCLE SUR LES POINTS DE GAUSS :
!     ------------------------------
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
!
        nx = zero
        ny = zero
        nz = zero
!
! ---   CALCUL DE LA NORMALE AU POINT DE GAUSS IPG :
!       ------------------------------------------
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
!
                nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
!
            end do
        end do
!
! ---   LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE :
!       ---------------------------------------------
        jac = sqrt(nx*nx+ny*ny+nz*nz)
!
        sigau = zr(ipoids+ipg-1)*jac
!
! ---   CALCUL DE AX, AY, AZ = SOMME(X.DS, Y.DS, Z.DS)
! ---   DANS LE REPERE GLOBAL ET DANS LE REPERE LOCAL :
!       ---------------------------------------------
        axggau = zero
        ayggau = zero
        azggau = zero
!
        axlgau = zero
        aylgau = zero
!
        do ino = 1, nno
!
            axggau = axggau+zr(ivf+ldec+ino-1)*xyzg(1, ino)
            ayggau = ayggau+zr(ivf+ldec+ino-1)*xyzg(2, ino)
            azggau = azggau+zr(ivf+ldec+ino-1)*xyzg(3, ino)
!
            axlgau = axlgau+zr(ivf+ldec+ino-1)*xyzl(1, ino)
            aylgau = aylgau+zr(ivf+ldec+ino-1)*xyzl(2, ino)
!
        end do
!
! ---     CALCUL DE  AXX, AYY, AZZ, AXY
! ---     = SOMME(X*X.DS, Y*Y.DS, Z*Z.DS, X*Y.DS) DANS LE REPERE LOCAL:
!         ------------------------------------------------------------
        xgau = zero
        ygau = zero
!
        do ino = 1, nno
!
            xgau = xgau+zr(ivf+ldec+ino-1)*xyzl(1, ino)
            ygau = ygau+zr(ivf+ldec+ino-1)*xyzl(2, ino)
        end do
!
        axxgau = xgau*xgau
        ayygau = ygau*ygau
        axygau = xgau*ygau
!
! ---      CALCUL DE LA SURFACE :
        aire = aire+sigau
! ---      AX :
        axg = axg+axggau*sigau
        axl = axl+axlgau*sigau
! ---      AY :
        ayg = ayg+ayggau*sigau
        ayl = ayl+aylgau*sigau
! ---      AZ :
        azg = azg+azggau*sigau
! ---      AXX :
        axx = axx+axxgau*sigau
! ---      AYY :
        ayy = ayy+ayygau*sigau
! ---      AXY :
        axy = axy+axygau*sigau
    end do
!
    if (abs(aire) .lt. r8prem()) then
        call utmess('F', 'ELEMENTS_48')
    end if
!
    s1 = un/aire
!
! --- COORDONNEES DU CENTRE GEOMETRIQUE G DE L'ELEMENT
! --- DANS LE REPERE GLOBAL.
! --- XG = AX/S, YG = AY/S, ZG = AZ/S :
!     -------------------------------
    xg = s1*axg
    yg = s1*ayg
    zg = s1*azg
!
! --- COORDONNEES DU CENTRE GEOMETRIQUE G DE L'ELEMENT
! --- DANS LE REPERE LOCAL :
!     --------------------
    xl = s1*axl
    yl = s1*ayl
!
! --- CALCUL DES TERMES DU TENSEUR D'INERTIE AU CENTRE GEOMETRIQUE
! --- DE L'ELEMENT ET DANS LE REPERE LOCAL :
!     -----------------------------------
! ---        IGXX = EPAIS*AYY + SOMME (Z**2.DV) -V*(YL**2 + ZL**2)
! ---   SOIT IGXX = EPAIS*AYY + S*(EPAIS**3)/12 -V*(YL**2)
! ---   EN FAIT ON CALCULE IGXX/EPAIS, LA MULTIPLICATION PAR
! ---   L'EPAISSEUR SERA FAITE EN FIN DE ROUTINE LORS DE LA
! ---   MULTIPLICATION PAR ROEP :
!       -----------------------
    igxx = ayy+aire*epais*epais/douze-aire*yl*yl
!
! ---        IGYY = EPAIS*AXX + SOMME (Z**2.DV) -V*(XL**2 + ZL**2) :
!            ----------------------------------------------------
    igyy = axx+aire*epais*epais/douze-aire*xl*xl
!
! ---        IGXY = EPAIS*AXY - V*XL*YL
!            ------------------------
    igxy = axy-aire*xl*yl
!
! ---        IGZZ = EPAIS*(AXX+AYY)  - V*(XL**2 + YL**2) :
!            -------------------------------------------
    igzz = axx+ayy-aire*(xl*xl+yl*yl)
!
! --- AFFECTATION DES TERMES DU TENSEUR D'INERTIE LOCAL :
!     -------------------------------------------------
!     MULTIPLICATION PAR MOINS DES TERMES EXTRA_DIAGONAUX
    matine(1) = roep*igxx
    matine(2) = -roep*igxy
    matine(3) = roep*igyy
    matine(4) = zero
    matine(5) = zero
    matine(6) = roep*igzz
!
! --- PASSAGE DU TENSEUR D'INERTIE DANS LE REPERE GLOBAL :
!     --------------------------------------------------
    call utpslg(1, 3, pgl, matine, inert0)
!
!     REMULTIPLICATION PAR MOINS DES TERMES EXTRA_DIAGONAUX
    inerti(1) = inert0(1)
    inerti(2) = inert0(3)
    inerti(3) = inert0(6)
    inerti(4) = -inert0(2)
    inerti(5) = -inert0(4)
    inerti(6) = -inert0(5)
!
! --- AFFECTATION DU VECTEUR DES COORDONNEES DU CENTRE DE GRAVITE
! --- DANS LE REPERE GLOBAL :
!     ---------------------
    cdg(1) = xg
    cdg(2) = yg
    cdg(3) = zg
!
! --- CALCUL DE LA MASSE DE L'ELEMENT :
!     -------------------------------
    mass = roep*aire
!
end subroutine
