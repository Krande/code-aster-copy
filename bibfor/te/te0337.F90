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
subroutine te0337(option, nomte)
!.......................................................................
!
!     BUT:
!       1) POUR L'OPTION : 'CARA_SECT_POUT3 ' :
!          CALCUL DU CHAMP ELEMENTAIRE A 10 COMPOSANTES :
!          SOMME/S_ELEMENT(DS,X.DS,Y.DS,Z.DS,X*X.DS,Y*Y.DS,Z*Z.DS,
!                             X*Y.DS,X*Z.DS,Y*Z.DS)
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          CES 10 QUANTITES GEOMETRIQUES SONT NOTEES :
!          A1 = S,AX,AY,AZ,AXX,AYY,AZZ,AXY,AXZ,AYZ
!
!       2) POUR L'OPTION : 'CARA_SECT_POUT4 ' :
!          CALCUL DES 2 VECTEURS DEFINIS AUX NOEUDS DES ELEMENTS
!          AYANT POURS VALEURS AU NOEUD I DE L'ELEMENT:
!          POUR LE PREMIER VECTEUR
!               SOMME/S_ELEMENT(NI.DS,0,0)
!          POUR LE SECOND VECTEUR
!               SOMME/S_ELEMENT(X*NI.DS,Y*NI.DS,Z*NI.DS)
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          AVEC X = XM - XG = NJ*XJ - XG
!               Y = YM - YG = NJ*YJ - YG
!               Z = ZM - ZG = NJ*ZJ - ZG
!          OU (XG,YG,ZG) SONT LES COORDONNEES DU CENTRE GEOMETRIQUE
!                        DU LIGREL DES MAILLES DE SURFACE TRAITE
!
!       3) POUR L'OPTION : 'CARA_SECT_POUT5 ' : 3D_TUYAU
!          CALCUL DU VECTEUR DEFINI AUX NOEUDS DES ELEMENTS
!          AYANT POUR VALEURS AU NOEUD I DE L'ELEMENT:
!          SOMME/S_ELEMENT(NI.COS(M.PHI).P.DS)
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    implicit none
#include "jeveux.h"
#include "asterfort/angvxy.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "blas/ddot.h"
    character(len=16) :: nomte, option
    real(kind=8) :: jacpoi, e1(3), e2(3)
    real(kind=8) :: jac, nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), zero
    real(kind=8) :: gp0(3), gpg(3), xpg(3), vsin(3)
    real(kind=8) :: norgp0, norgpg, angl(3), pgl(3, 3)
    real(kind=8) :: cosphi, sinphi, phi, cosmfi, sinmfi, phi0, ray
    real(kind=8) :: sigau, axgau, aygau, azgau, xgau, ygau, zgau
    real(kind=8) :: axxgau, ayygau, azzgau, axygau, axzgau, ayzgau
    integer(kind=8) :: ipoids, ivf, idfdx, idfdy, igeom, m
    integer(kind=8) :: ndim, nno, ipg, npg1, inumod, nnos, jgano
    integer(kind=8) :: idec, jdec, kdec, ldec
    integer(kind=8) :: ino, jno, ii
    integer(kind=8) :: i, j, isect, iorig, iorifi, iaxe
    integer(kind=8) :: ivect1, ivect2, ivect3, ivect4, ivect5, ivect6
    blas_int :: b_incx, b_incy, b_n
!
!
!
!
    zero = 0.0d0
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
!
    if (option .eq. 'CARA_SECT_POUT3') then
        call jevech('PCASECT', 'E', isect)
!
        do i = 1, 10
            zr(isect+i-1) = 0.0d0
        end do
!
    else if (option .eq. 'CARA_SECT_POU3R') then
        call jevech('PORIGIN', 'L', iorig)
        call jevech('PRAYONM', 'E', isect)
        zr(isect+1-1) = zero
!
    else if (option .eq. 'CARA_SECT_POUT4') then
        call jevech('PORIGIN', 'L', iorig)
        call jevech('PVECTU1', 'E', ivect1)
        call jevech('PVECTU2', 'E', ivect2)
!
        do i = 1, 3*nno
            zr(ivect1+i-1) = zero
            zr(ivect2+i-1) = zero
        end do
!
    else if (option .eq. 'CARA_SECT_POUT5') then
        call jevech('PORIGIN', 'L', iorig)
        call jevech('PORIGFI', 'L', iorifi)
        call jevech('PNUMMOD', 'L', inumod)
        call jevech('PCAORIE', 'L', iaxe)
        call jevech('PVECTU1', 'E', ivect1)
        call jevech('PVECTU2', 'E', ivect2)
        call jevech('PVECTU3', 'E', ivect3)
        call jevech('PVECTU4', 'E', ivect4)
        call jevech('PVECTU5', 'E', ivect5)
        call jevech('PVECTU6', 'E', ivect6)
!
        e1(1) = zr(iaxe+1-1)
        e1(2) = zr(iaxe+2-1)
        e1(3) = zr(iaxe+3-1)
!
!         COORDONNES DU POINT P TEL QUE GP EST L'ORIGINE
!         DE L'ANGLE PHI
!
        do i = 1, 3
            gp0(i) = zr(iorifi-1+i)-zr(iorig-1+i)
        end do
        call normev(gp0, norgp0)
!
!         NUMERO DE MODE DE FOURIER
!
        m = zi(inumod)
        do i = 1, 3*nno
            zr(ivect1+i-1) = zero
            zr(ivect2+i-1) = zero
            zr(ivect3+i-1) = zero
            zr(ivect4+i-1) = zero
            zr(ivect5+i-1) = zero
            zr(ivect6+i-1) = zero
        end do
    end if
!
!--- CALCUL DES PRODUITS VECTORIELS OMI X OMJ
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
!    ---------------------------
!--- - OPTION : CARA_SECT_POUT3-
!    ---------------------------
!
    if (option .eq. 'CARA_SECT_POUT3') then
!
!---  BOUCLE SUR LES POINTS DE GAUSS
!     ------------------------------
!
        do ipg = 1, npg1
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
!
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
!
!---    CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!       ------------------------------------------
!
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
!---  LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE
!     ---------------------------------------------
!
            jac = sqrt(nx*nx+ny*ny+nz*nz)
!
            sigau = zr(ipoids+ipg-1)*jac
!
!---  CALCUL DE AX, AY, AZ = SOMME(X.DS, Y.DS, Z.DS)
!     ----------------------------------------------
!
            axgau = zero
            aygau = zero
            azgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
!
                axgau = axgau+zr(ivf+ldec+ino-1)*zr(i+1)
                aygau = aygau+zr(ivf+ldec+ino-1)*zr(i+2)
                azgau = azgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
!---   CALCUL DE  AXX, AYY, AZZ, AXY, AXZ, AYZ
!---   = SOMME(X*X.DS, Y*Y.DS, Z*Z.DS, X*Y.DS, X*Z.DS, Y*Z.DS)
!      -------------------------------------------------------
!
            xgau = zero
            ygau = zero
            zgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
!
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
                zgau = zgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
            axxgau = xgau*xgau
            ayygau = ygau*ygau
            azzgau = zgau*zgau
            axygau = xgau*ygau
            axzgau = xgau*zgau
            ayzgau = ygau*zgau
!
!---  CALCUL DE A1 = S
            zr(isect+1-1) = zr(isect+1-1)+sigau
!---  AX
            zr(isect+2-1) = zr(isect+2-1)+axgau*sigau
!---  AY
            zr(isect+3-1) = zr(isect+3-1)+aygau*sigau
!---  AZ
            zr(isect+4-1) = zr(isect+4-1)+azgau*sigau
!---  AXX
            zr(isect+5-1) = zr(isect+5-1)+axxgau*sigau
!---  AYY
            zr(isect+6-1) = zr(isect+6-1)+ayygau*sigau
!---  AZZ
            zr(isect+7-1) = zr(isect+7-1)+azzgau*sigau
!---  AXY
            zr(isect+8-1) = zr(isect+8-1)+axygau*sigau
!---  AXZ
            zr(isect+9-1) = zr(isect+9-1)+axzgau*sigau
!---  AYZ
            zr(isect+10-1) = zr(isect+10-1)+ayzgau*sigau
!
        end do
!--- FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
!---  ET FIN DE L'OPTION 'CARA_SECT_POUT3'
!
!    ---------------------------
!--- - OPTION : CARA_SECT_POU3R-
!    ---------------------------
!
    else if (option .eq. 'CARA_SECT_POU3R') then
!
!
!---  BOUCLE SUR LES POINTS DE GAUSS
!     ------------------------------
!
        do ipg = 1, npg1
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
!
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
!
!---    CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!       ------------------------------------------
!
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
!---  LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE
!     ---------------------------------------------
!
            jac = sqrt(nx*nx+ny*ny+nz*nz)
            sigau = zr(ipoids+ipg-1)*jac
            xgau = zero
            ygau = zero
            zgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
!
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
                zgau = zgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
            ray = (xgau-zr(iorig+1-1))**2.d0+(ygau-zr(iorig+2-1))**2.d0+(zgau-zr(iorig+3-1) &
                                                                         )**2.d0
            ray = sqrt(ray)
!
!---  CALCUL DE A1 = RAYON
            zr(isect+1-1) = zr(isect+1-1)+ray*sigau
!
        end do
!--- FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
!---  ET FIN DE L'OPTION 'CARA_SECT_POU3R'
!
!    ---------------------------
!--- - OPTION : CARA_SECT_POUT4-
!    ---------------------------
!
    else if (option .eq. 'CARA_SECT_POUT4') then
!
!---  BOUCLE SUR LES POINTS DE GAUSS
!     ------------------------------
!
        do ipg = 1, npg1
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
!
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
!
!---    CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!       ------------------------------------------
!
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
!---  LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE
!     ---------------------------------------------
!
            jac = sqrt(nx*nx+ny*ny+nz*nz)
!
            sigau = zr(ipoids+ipg-1)*jac
!
!---  CALCUL DE VECT1(I) = SOMME(NI.DS, 0, 0)
!     ---------------------------------------
!
            do ino = 1, nno
!
                zr(ivect1+3*(ino-1)+1-1) = zr(ivect1+3*(ino-1)+1-1)+zr(ivf+ldec+ino-1)*sigau
            end do
!
!---  CALCUL DE VECT2(I) = SOMME(X*NI.DS, Y*NI.DS, Z*NI.DS)
!     -----------------------------------------------------
!
            xgau = zero
            ygau = zero
            zgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
!
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
                zgau = zgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
            do ino = 1, nno
!
                zr(ivect2+3*(ino-1)+1-1) = zr( &
                                           ivect2+3*(ino-1)+1-1)+zr(ivf+ldec+ino-1)*(xgau-zr(iori&
                                           &g+1-1) &
                                           )*sigau
!
                zr(ivect2+3*(ino-1)+2-1) = zr( &
                                           ivect2+3*(ino-1)+2-1)+zr(ivf+ldec+ino-1)*(ygau-zr(iori&
                                           &g+2-1) &
                                           )*sigau
!
                zr(ivect2+3*(ino-1)+3-1) = zr( &
                                           ivect2+3*(ino-1)+3-1)+zr(ivf+ldec+ino-1)*(zgau-zr(iori&
                                           &g+3-1) &
                                           )*sigau
            end do
!
        end do
!---  FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
!---   ET FIN DE L'OPTION 'CARA_SECT_POUT4'
!
!     ---------------------------
! --- - OPTION : CARA_SECT_POUT5-
!     ---------------------------
!
    else if (option .eq. 'CARA_SECT_POUT5') then
!
!
        do ipg = 1, npg1
            kdec = (ipg-1)*nno*ndim
            ldec = (ipg-1)*nno
!
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
!
!---    CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
!       ------------------------------------------
!
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
!---  LE JACOBIEN EST EGAL A LA NORME DE LA NORMALE
!     ---------------------------------------------
!
            jac = sqrt(nx*nx+ny*ny+nz*nz)
            jacpoi = jac*zr(ipoids+ipg-1)
!
! ---   COORDONNEES DU POINT D'INTEGRATION COURANT :
!       ------------------------------------------
            do ii = 1, 3
                xpg(ii) = zero
            end do
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
                xpg(1) = xpg(1)+zr(ivf+ldec+ino-1)*zr(i+1)
                xpg(2) = xpg(2)+zr(ivf+ldec+ino-1)*zr(i+2)
                xpg(3) = xpg(3)+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
!  CALCUL DU VECTEUR G-PG ET DE L'ANGLE PHI ENTRE G-P0 ET G-PG
!
            do i = 1, 3
                gpg(i) = xpg(i)-zr(iorig-1+i)
            end do
            call normev(gpg, norgpg)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            cosphi = ddot(b_n, gp0, b_incx, gpg, b_incy)
!PM          CALL PROVEC(GP0,GPG,VSIN)
            call provec(gpg, gp0, vsin)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinphi = ddot(b_n, e1, b_incx, vsin, b_incy)
            phi0 = atan2(sinphi, cosphi)
!JMP          PHI=-PHI0
            phi = phi0
            cosmfi = cos(m*phi)
            sinmfi = sin(m*phi)
!
!  CALCUL DE PGL MATRICE DE PASSAGE DE X,Y,Z GLOBAL A E1,E2,E3
!
            call provec(gpg, e1, e2)
            call angvxy(e1, e2, angl)
            call matrot(angl, pgl)
!
            do ino = 1, nno
                do ii = 1, 3
!
! CALCUL DE VECT1(I) : TERMES EN UMI(COS(M.PHI)) ET UMO (SIN(M.PHI))
!
                    zr(ivect1+3*(ino-1)+ii-1) = &
                        zr(ivect1+3*(ino-1)+ii-1)+cosmfi*pgl(1, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect4+3*(ino-1)+ii-1) = &
                        zr(ivect4+3*(ino-1)+ii-1)+sinmfi*pgl(1, ii)*zr(ivf+ldec+ino-1)*jacpoi
!
! CALCUL DE VECT2(I) : TERMES EN VMI(COS(M.PHI)) ET VMO (SIN(M.PHI))
!
                    zr(ivect2+3*(ino-1)+ii-1) = &
                        zr(ivect2+3*(ino-1)+ii-1)+cosmfi*pgl(2, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect5+3*(ino-1)+ii-1) = &
                        zr(ivect5+3*(ino-1)+ii-1)+sinmfi*pgl(2, ii)*zr(ivf+ldec+ino-1)*jacpoi
!
! CALCUL DE VECT3(I) : TERMES EN WMI(COS(M.PHI)) ET WMO (SIN(M.PHI))
!
                    zr(ivect3+3*(ino-1)+ii-1) = &
                        zr(ivect3+3*(ino-1)+ii-1)+cosmfi*pgl(3, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect6+3*(ino-1)+ii-1) = &
                        zr(ivect6+3*(ino-1)+ii-1)+sinmfi*pgl(3, ii)*zr(ivf+ldec+ino-1)*jacpoi
                end do
            end do
        end do
! ---  FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
! ---  ET FIN DE L'OPTION 'CARA_SECT_POUT5'
    end if
end subroutine
