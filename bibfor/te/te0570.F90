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
subroutine te0570(option, nomte)
!.......................................................................
!
!     BUT:
!       1) POUR L'OPTION : 'CARA_SECT_POUT3 ' :
!          CALCUL DU CHAMP ELEMENTAIRE A 10 COMPOSANTES :
!          SOMME/S_ELEMENT(DS,X.DS,Y.DS,Z.DS,X*X.DS,Y*Y.DS,Z*Z.DS,
!                             X*Y.DS,X*Z.DS,Y*Z.DS)
!          SUR LES ELEMENTS DE BORD DE COQUE :
!          MEBODKT, MEBODST, MEBOQ4G, MEBOCQ3
!
!          CES 10 QUANTITES GEOMETRIQUES SONT NOTEES :
!          A1 = S,AX,AY,AZ,AXX,AYY,AZZ,AXY,AXZ,AYZ
!
!       2) POUR L'OPTION : 'CARA_SECT_POUT4 ' :
!          CALCUL DU VECTEUR DEFINIS AUX NOEUDS DES ELEMENTS
!          AYANT POURS VALEURS AU NOEUD I DE L'ELEMENT:
!          SOMME/S_ELEMENT(X*NI.DS,Y*NI.DS,Z*NI.DS,NI.DS,NI.DS*H3/12,0)
!
!          SUR LES ELEMENTS DE BORD DE COQUE :
!          MEBODKT, MEBODST, MEBOQ4G, MEBOCQ3
!
!          AVEC X = XM - XG = NJ*XJ - XG
!               Y = YM - YG = NJ*YJ - YG
!               Z = ZM - ZG = NJ*ZJ - ZG
!          OU (XG,YG,ZG) SONT LES COORDONNEES DU CENTRE GEOMETRIQUE
!                        DU LIGREL DES MAILLES DE BORD DE COQUE TRAITE
!
!       3) POUR L'OPTION : 'CARA_SECT_POUT5 ' : COQ_TUYAU
!          CALCUL DU VECTEUR DEFINI AUX NOEUDS DES ELEMENTS
!          AYANT POUR VALEURS AU NOEUD I DE L'ELEMENT:
!          SOMME/S_ELEMENT(NI.COS(M.PHI).P.DS)
!          SUR LES ELEMENTS DE BORD DE COQUE :
!          MEBODKT
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/angvxy.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
    character(len=8) :: elrefe
    character(len=16) :: nomte, option
    real(kind=8) :: jac, jacpoi, jacpo2, zero, e1(3), e2(3), e3(3), gn1(3)
    real(kind=8) :: xg, yg, zg, gp0(3), gpg(3), xpg(3), xn1(3), vsin(3)
    real(kind=8) :: norgp0, norgpg, angl(3), pgl(3, 3), pi, rayon
    real(kind=8) :: epais, coef, dxdk, dydk, dzdk, axgau, aygau
    real(kind=8) :: azgau, xgau, ygau, zgau, axxgau, ayygau, azzgau, sinphi
    real(kind=8) :: axygau, axzgau, ayzgau, e3xx, e3xy, e3xz, e3yy, e3yz
    real(kind=8) :: e3zz, cosphi, phi, cosmfi, sinmfi, phi0
    integer(kind=8) :: nno, nnos, jgano, ndim, ipg, npg, idfdk, icoqu, iopt
    integer(kind=8) :: ldec, iorifi, m, isect, i, iorig, ivect1, inumod, iret
    integer(kind=8) :: ivect2, ivect3, iaxe, ino, ii, ipoids, ivf, igeom, itabm(8)
    blas_int :: b_incx, b_incy, b_n
!
!
!
    call elref1(elrefe)
!
    zero = 0.0d0
    iopt = 0
    pi = r8pi()
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
!
!
! --- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!     ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! --- RECUPERATION L'EPAISSEUR DE L'ELEMENT :
!     -------------------------------------
    call tecach('OON', 'PCACOQU', 'L', iret, nval=8, &
                itab=itabm)
!
    if (.not. zl(itabm(8))) then
        call utmess('F', 'ELEMENTS4_32')
    end if
!
    icoqu = itabm(1)
    epais = zr(icoqu+1-1)
    coef = epais*epais*epais/12.0d0
!
    if (epais .le. r8prem()) then
        call utmess('F', 'ELEMENTS4_33')
    end if
!
    if (option .eq. 'CARA_SECT_POUT3') then
        call jevech('PCASECT', 'E', isect)
        iopt = 3
        do i = 1, 10
            zr(isect+i-1) = zero
        end do
!
    else if (option .eq. 'CARA_SECT_POUT4') then
        call jevech('PORIGIN', 'L', iorig)
        call jevech('PVECTU1', 'E', ivect1)
        call jevech('PVECTU2', 'E', ivect2)
        iopt = 4
        xg = zr(iorig+1-1)
        yg = zr(iorig+2-1)
        zg = zr(iorig+3-1)
!
        do i = 1, 6*nno
            zr(ivect1+i-1) = zero
        end do
!
    else if (option .eq. 'CARA_SECT_POUT5') then
        call jevech('PORIGIN', 'L', iorig)
        call jevech('PORIGFI', 'L', iorifi)
        call jevech('PNUMMOD', 'L', inumod)
        call jevech('PVECTU1', 'E', ivect1)
        call jevech('PVECTU2', 'E', ivect2)
        call jevech('PVECTU3', 'E', ivect3)
        iopt = 5
        xg = zr(iorig+1-1)
        yg = zr(iorig+2-1)
        zg = zr(iorig+3-1)
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
        do i = 1, 6*nno
            zr(ivect1+i-1) = zero
            zr(ivect2+i-1) = zero
            zr(ivect3+i-1) = zero
        end do
    end if
!
    call jevech('PCAORIE', 'L', iaxe)
    e1(1) = zr(iaxe+1-1)
    e1(2) = zr(iaxe+2-1)
    e1(3) = zr(iaxe+3-1)
!
!     ---------------------------
! --- - OPTION : CARA_SECT_POUT3-
!     ---------------------------
!
    if (iopt .eq. 3) then
!
! --- BOUCLE SUR LES POINTS DE GAUSS :
!     ------------------------------
        do ipg = 1, npg
!
            ldec = (ipg-1)*nno
!
            dxdk = zero
            dydk = zero
            dzdk = zero
!
! ---   DERIVEES DES FONCTION DE FORME SUR L'ELEMENT REEL :
!       -------------------------------------------------
            do i = 1, nno
                dxdk = dxdk+zr(igeom+3*(i-1)+1-1)*zr(idfdk+ldec+i-1)
                dydk = dydk+zr(igeom+3*(i-1)+2-1)*zr(idfdk+ldec+i-1)
                dzdk = dzdk+zr(igeom+3*(i-1)+3-1)*zr(idfdk+ldec+i-1)
            end do
!
! ---   JACOBIEN :
!       --------
            jac = sqrt(dxdk*dxdk+dydk*dydk+dzdk*dzdk)
            if (jac .le. r8prem()) then
                call utmess('F', 'ELEMENTS4_34')
            end if
            jacpoi = jac*zr(ipoids+ipg-1)*epais
            jacpo2 = jac*zr(ipoids+ipg-1)*coef
!
! ---   CALCUL DU VECTEUR E2 TANGENT A LA FIBRE MOYENNE AU POINT
! ---   D'INTEGRATION COURANT :
!
            e2(1) = dxdk/jac
            e2(2) = dydk/jac
            e2(3) = dzdk/jac
!
! ---   CALCUL DU VECTEUR E3 NORMAL A E1 ET E2
!
            e3(1) = e1(2)*e2(3)-e1(3)*e2(2)
            e3(2) = e1(3)*e2(1)-e1(1)*e2(3)
            e3(3) = e1(1)*e2(2)-e1(2)*e2(1)
!
! ---   CALCUL DE AX, AY, AZ = SOMME(X.DS, Y.DS, Z.DS) :
!       ----------------------------------------------
            axgau = zero
            aygau = zero
            azgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
                axgau = axgau+zr(ivf+ldec+ino-1)*zr(i+1)
                aygau = aygau+zr(ivf+ldec+ino-1)*zr(i+2)
                azgau = azgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
! ---   CALCUL DE  AXX, AYY, AZZ, AXY, AXZ, AYZ
! ---   = SOMME(X*X.DS, Y*Y.DS, Z*Z.DS, X*Y.DS, X*Z.DS, Y*Z.DS) :
!       -------------------------------------------------------
            xgau = zero
            ygau = zero
            zgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
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
! ---   CALCUL DES TERMES EN E3*E3
!
            e3xx = e3(1)*e3(1)*jacpo2
            e3xy = e3(1)*e3(2)*jacpo2
            e3xz = e3(1)*e3(3)*jacpo2
            e3yy = e3(2)*e3(2)*jacpo2
            e3yz = e3(2)*e3(3)*jacpo2
            e3zz = e3(3)*e3(3)*jacpo2
!
!---  CALCUL DE A1 = S
            zr(isect+1-1) = zr(isect+1-1)+jacpoi
!---  AX
            zr(isect+2-1) = zr(isect+2-1)+axgau*jacpoi
!---  AY
            zr(isect+3-1) = zr(isect+3-1)+aygau*jacpoi
!---  AZ
            zr(isect+4-1) = zr(isect+4-1)+azgau*jacpoi
!---  AXX
            zr(isect+5-1) = zr(isect+5-1)+axxgau*jacpoi+e3xx
!---  AYY
            zr(isect+6-1) = zr(isect+6-1)+ayygau*jacpoi+e3yy
!---  AZZ
            zr(isect+7-1) = zr(isect+7-1)+azzgau*jacpoi+e3zz
!---  AXY
            zr(isect+8-1) = zr(isect+8-1)+axygau*jacpoi+e3xy
!---  AXZ
            zr(isect+9-1) = zr(isect+9-1)+axzgau*jacpoi+e3xz
!---  AYZ
            zr(isect+10-1) = zr(isect+10-1)+ayzgau*jacpoi+e3yz
!
        end do
! --- FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
! --- ET FIN DE L'OPTION 'CARA_SECT_POUT3'
!
!     ---------------------------
! --- - OPTION : CARA_SECT_POUT4-
!     ---------------------------
!
    else if (iopt .eq. 4) then
!
! --- BOUCLE SUR LES POINTS DE GAUSS :
!     ------------------------------
!
        do ipg = 1, npg
!
            ldec = (ipg-1)*nno
!
            dxdk = zero
            dydk = zero
            dzdk = zero
!
! ---   DERIVEES DES FONCTION DE FORME SUR L'ELEMENT REEL :
!       -------------------------------------------------
            do i = 1, nno
                dxdk = dxdk+zr(igeom+3*(i-1)+1-1)*zr(idfdk+ldec+i-1)
                dydk = dydk+zr(igeom+3*(i-1)+2-1)*zr(idfdk+ldec+i-1)
                dzdk = dzdk+zr(igeom+3*(i-1)+3-1)*zr(idfdk+ldec+i-1)
            end do
!
! ---   JACOBIEN :
!       --------
            jac = sqrt(dxdk*dxdk+dydk*dydk+dzdk*dzdk)
            if (jac .le. r8prem()) then
                call utmess('F', 'ELEMENTS4_34')
            end if
            jacpoi = jac*zr(ipoids+ipg-1)*epais
            jacpo2 = jac*zr(ipoids+ipg-1)*coef
!
! ---   CALCUL DU VECTEUR E2 TANGENT A LA FIBRE MOYENNE AU POINT
! ---   D'INTEGRATION COURANT :
!       ---------------------
            e2(1) = dxdk/jac
            e2(2) = dydk/jac
            e2(3) = dzdk/jac
!
! ---   CALCUL DU VECTEUR E3 NORMAL A E1 ET E2
!
            e3(1) = e1(2)*e2(3)-e1(3)*e2(2)
            e3(2) = e1(3)*e2(1)-e1(1)*e2(3)
            e3(3) = e1(1)*e2(2)-e1(2)*e2(1)
!
!
! ---   CALCUL DES TERMES EN E3*E3
!
            e3xx = e3(1)*e3(1)*jacpo2
            e3xy = e3(1)*e3(2)*jacpo2
            e3xz = e3(1)*e3(3)*jacpo2
            e3yy = e3(2)*e3(2)*jacpo2
            e3yz = e3(2)*e3(3)*jacpo2
            e3zz = e3(3)*e3(3)*jacpo2
!
! ---   COORDONNEES DU POINT D'INTEGRATION COURANT :
!       ------------------------------------------
            xgau = zero
            ygau = zero
            zgau = zero
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
                zgau = zgau+zr(ivf+ldec+ino-1)*zr(i+3)
            end do
!
! --- CALCUL DE VECT1(I)
!
            do ino = 1, nno
                i = igeom+3*(ino-1)-1
!
                zr(ivect1+6*(ino-1)+1-1) = zr(ivect1+6*(ino-1)+1-1)+zr(ivf+ldec+ino-1)*(xgau-xg &
                                                                                        )*jacpoi
                zr(ivect1+6*(ino-1)+2-1) = zr(ivect1+6*(ino-1)+2-1)+zr(ivf+ldec+ino-1)*(ygau-yg &
                                                                                        )*jacpoi
                zr(ivect1+6*(ino-1)+3-1) = zr(ivect1+6*(ino-1)+3-1)+zr(ivf+ldec+ino-1)*(zgau-zg &
                                                                                        )*jacpoi
                zr(ivect1+6*(ino-1)+4-1) = zr(ivect1+6*(ino-1)+4-1)+zr(ivf+ldec+ino-1)*jacpoi
                zr(ivect1+6*(ino-1)+5-1) = zr(ivect1+6*(ino-1)+5-1)+zr(ivf+ldec+ino-1)*jacpo2
!
!            PRODUIT VECTORIEL N.(THETA.N).
!
                zr(ivect2+6*(ino-1)+1-1) = zr( &
                                           ivect2+6*(ino-1)+1-1)+zr(ivf+ldec+ino-1)*(e3yy+e3zz)
                zr(ivect2+6*(ino-1)+2-1) = zr(ivect2+6*(ino-1)+2-1)-zr(ivf+ldec+ino-1)*e3xy
                zr(ivect2+6*(ino-1)+3-1) = zr(ivect2+6*(ino-1)+3-1)-zr(ivf+ldec+ino-1)*e3xz
                zr(ivect2+6*(ino-1)+4-1) = zr( &
                                           ivect2+6*(ino-1)+4-1)+zr(ivf+ldec+ino-1)*(e3zz+e3xx)
                zr(ivect2+6*(ino-1)+5-1) = zr(ivect2+6*(ino-1)+5-1)-zr(ivf+ldec+ino-1)*e3yz
                zr(ivect2+6*(ino-1)+6-1) = zr( &
                                           ivect2+6*(ino-1)+6-1)+zr(ivf+ldec+ino-1)*(e3yy+e3xx)
            end do
        end do
!
! ---  FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
! ---  ET FIN DE L'OPTION 'CARA_SECT_POUT4'
!     ---------------------------
! --- - OPTION : CARA_SECT_POUT5-
!     ---------------------------
!
    else if (iopt .eq. 5) then
        do ipg = 1, npg
            dxdk = zero
            dydk = zero
            dzdk = zero
!
! ---   DERIVEES DES FONCTION DE FORME SUR L'ELEMENT REEL :
!       -------------------------------------------------
            ldec = (ipg-1)*nno
            do i = 1, nno
                dxdk = dxdk+zr(igeom+3*(i-1)+1-1)*zr(idfdk+ldec+i-1)
                dydk = dydk+zr(igeom+3*(i-1)+2-1)*zr(idfdk+ldec+i-1)
                dzdk = dzdk+zr(igeom+3*(i-1)+3-1)*zr(idfdk+ldec+i-1)
            end do
!
! ---   CALCUL DU RAYON
!
            xn1(1) = zr(igeom+1-1)
            xn1(2) = zr(igeom+2-1)
            xn1(3) = zr(igeom+3-1)
            do i = 1, 3
                gn1(i) = xn1(i)-zr(iorig-1+i)
            end do
            call normev(gn1, rayon)
!
! ---   JACOBIEN :
!       --------
            jac = sqrt(dxdk*dxdk+dydk*dydk+dzdk*dzdk)
            jacpoi = jac*zr(ipoids+ipg-1)
            jacpoi = jacpoi/rayon/pi
!
! ---   COORDONNEES DU POINT D'INTEGRATION COURANT :
!       ------------------------------------------
            xpg = zero
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
            call provec(gpg, gp0, vsin)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            sinphi = ddot(b_n, e1, b_incx, vsin, b_incy)
            phi0 = atan2(sinphi, cosphi)
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
                    zr(ivect1+6*(ino-1)+ii-1) = &
                        zr(ivect1+6*(ino-1)+ii-1)+cosmfi*pgl(1, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect1+6*(ino-1)+3+ii-1) = &
                        zr(ivect1+6*(ino-1)+3+ii-1)+sinmfi*pgl(1, ii)*zr(ivf+ldec+ino-1)*jacpoi
!
! CALCUL DE VECT2(I) : TERMES EN VMI(COS(M.PHI)) ET VMO (SIN(M.PHI))
!
                    zr(ivect2+6*(ino-1)+ii-1) = &
                        zr(ivect2+6*(ino-1)+ii-1)+cosmfi*pgl(2, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect2+6*(ino-1)+3+ii-1) = &
                        zr(ivect2+6*(ino-1)+3+ii-1)+sinmfi*pgl(2, ii)*zr(ivf+ldec+ino-1)*jacpoi
!
! CALCUL DE VECT3(I) : TERMES EN WMI(COS(M.PHI)) ET WMO (SIN(M.PHI))
!
                    zr(ivect3+6*(ino-1)+ii-1) = &
                        zr(ivect3+6*(ino-1)+ii-1)+cosmfi*pgl(3, ii)*zr(ivf+ldec+ino-1)*jacpoi
                    zr(ivect3+6*(ino-1)+3+ii-1) = &
                        zr(ivect3+6*(ino-1)+3+ii-1)+sinmfi*pgl(3, ii)*zr(ivf+ldec+ino-1)*jacpoi
                end do
            end do
        end do
    end if
!
end subroutine
