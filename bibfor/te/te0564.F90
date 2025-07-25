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
subroutine te0564(option, nomte)
!.......................................................................
!
!     BUT:
!       1) POUR L'OPTION : 'CARA_SECT_POUT3 ' :
!          CALCUL DU CHAMP ELEMENTAIRE A 6 COMPOSANTES :
!          SOMME/S_ELEMENT(DS,X.DS,Y.DS,X*X.DS,Y*Y.DS,X*Y.DS)
!          SUR LES ELEMENTS DE BORD DES ELEMENTS 2D
!
!          CES 6 QUANTITES GEOMETRIQUES SONT NOTEES :
!          A1 = S,AX,AY,AXX,AYY,AXY
!
!       2) POUR L'OPTION : 'CARA_SECT_POUT4 ' :
!          CALCUL DES 2 VECTEURS DEFINIS AUX NOEUDS DES ELEMENTS
!          AYANT POURS VALEURS AU NOEUD I DE L'ELEMENT:
!          POUR LE PREMIER VECTEUR
!               SOMME/S_ELEMENT(NI.DS,0,0)
!          POUR LE SECOND VECTEUR
!               SOMME/S_ELEMENT(X*NI.DS,Y*NI.DS)
!          SUR LES ELEMENTS DE BORD DES ELEMENTS 2D
!
!          AVEC X = XM - XG = NJ*XJ - XG
!               Y = YM - YG = NJ*YJ - YG
!               Z = ZM - ZG = NJ*ZJ - ZG
!          OU (XG,YG,ZG) SONT LES COORDONNEES DU CENTRE GEOMETRIQUE
!                        DU LIGREL DES MAILLES DE BORD (SEG) TRAITEES
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: elrefe
    character(len=16) :: nomte, option
    real(kind=8) :: jac, jacpoi, zero, r
    real(kind=8) :: xg, yg
    real(kind=8) :: dxdk, dydk, axgau, aygau
    real(kind=8) :: xgau, ygau, axxgau, ayygau
    real(kind=8) :: axygau
    integer(kind=8) :: nno, nnos, jgano, ndim, ipg, npg, idfdk, iopt
    integer(kind=8) :: ldec, isect, i, iorig, ivect1
    integer(kind=8) :: ivect2, ino, ipoids, ivf, igeom
    aster_logical :: laxi
!
!
!
!
    call elref1(elrefe)
!
    zero = 0.0d0
    iopt = 0
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
! --- RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!     ----------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
    if (option .eq. 'CARA_SECT_POUT3') then
        call jevech('PCASECT', 'E', isect)
        iopt = 3
        do i = 1, 6
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
        do i = 1, 2*nno
            zr(ivect1+i-1) = zero
            zr(ivect2+i-1) = zero
        end do
!
    end if
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
!
! ---   DERIVEES DES FONCTION DE FORME SUR L'ELEMENT REEL :
!       -------------------------------------------------
            do ino = 1, nno
                i = igeom+2*(ino-1)-1
                dxdk = dxdk+zr(i+1)*zr(idfdk+ldec+ino-1)
                dydk = dydk+zr(i+2)*zr(idfdk+ldec+ino-1)
            end do
!
! ---   JACOBIEN :
!       --------
            jac = sqrt(dxdk*dxdk+dydk*dydk)
            if (jac .le. r8prem()) then
                call utmess('F', 'ELEMENTS4_34')
            end if
!
!           dans le cas AXIS la prise en compte du rayon ne doit pas intervenir
!           a cette etape
            jacpoi = jac*zr(ipoids+ipg-1)
!
! ---   CALCUL DE AX, AY = SOMME(X.DS, Y.DS) :
!       ----------------------------------------------
            axgau = zero
            aygau = zero
!
            do ino = 1, nno
                i = igeom+2*(ino-1)-1
                axgau = axgau+zr(ivf+ldec+ino-1)*zr(i+1)
                aygau = aygau+zr(ivf+ldec+ino-1)*zr(i+2)
            end do
!
! ---   CALCUL DE  AXX, AYY, AXY = SOMME(X*X.DS, Y*Y.DS, X*Y.DS) :
!       -------------------------------------------------------
            xgau = zero
            ygau = zero
!
            do ino = 1, nno
                i = igeom+2*(ino-1)-1
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
            end do
!
            axxgau = xgau*xgau
            ayygau = ygau*ygau
            axygau = xgau*ygau
!
!---  CALCUL DE A1 = S
            zr(isect+1-1) = zr(isect+1-1)+jacpoi
!---  AX
            zr(isect+2-1) = zr(isect+2-1)+axgau*jacpoi
!---  AY
            zr(isect+3-1) = zr(isect+3-1)+aygau*jacpoi
!---  AXX
            zr(isect+4-1) = zr(isect+4-1)+axxgau*jacpoi
!---  AYY
            zr(isect+5-1) = zr(isect+5-1)+ayygau*jacpoi
!---  AXY
            zr(isect+6-1) = zr(isect+6-1)+axygau*jacpoi
!
        end do
!
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
        do ipg = 1, npg
!
            ldec = (ipg-1)*nno
!
            dxdk = zero
            dydk = zero
!
! ---   DERIVEES DES FONCTION DE FORME SUR L'ELEMENT REEL :
!       -------------------------------------------------
            do ino = 1, nno
                i = igeom+2*(ino-1)-1
                dxdk = dxdk+zr(i+1)*zr(idfdk+ldec+ino-1)
                dydk = dydk+zr(i+2)*zr(idfdk+ldec+ino-1)
            end do
!
! ---   JACOBIEN :
!       --------
            jac = sqrt(dxdk*dxdk+dydk*dydk)
            if (jac .le. r8prem()) then
                call utmess('F', 'ELEMENTS4_34')
            end if
!
            if (laxi) then
                r = 0.d0
                do ino = 1, nno
                    r = r+zr(igeom+2*(ino-1))*zr(ivf+ldec+ino-1)
                end do
                jac = jac*r
            end if
            jacpoi = jac*zr(ipoids+ipg-1)
!
!---    CALCUL DE VECT1(I) = SOMME(NI.DS, 0)
!       ---------------------------------------
            do ino = 1, nno
                zr(ivect1+2*(ino-1)+1-1) = zr(ivect1+2*(ino-1)+1-1)+zr(ivf+ldec+ino-1)*jacpoi
            end do
!
!---    CALCUL DE VECT2(I) = SOMME(X*NI.DS, Y*NI.DS)
!       --------------------------------------------
            xgau = zero
            ygau = zero
!
            do ino = 1, nno
                i = igeom+2*(ino-1)-1
                xgau = xgau+zr(ivf+ldec+ino-1)*zr(i+1)
                ygau = ygau+zr(ivf+ldec+ino-1)*zr(i+2)
            end do
!
            do ino = 1, nno
                zr(ivect2+2*(ino-1)) = zr(ivect2+2*(ino-1))+zr(ivf+ldec+ino-1)*(xgau-xg &
                                                                                )*jacpoi
                zr(ivect2+2*(ino-1)+1) = zr(ivect2+2*(ino-1)+1)+zr(ivf+ldec+ino-1)*(ygau-yg &
                                                                                    )*jacpoi
            end do
!
!
!
        end do
!
! ---  FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION
! ---  ET FIN DE L'OPTION 'CARA_SECT_POUT4'
!
    end if
!
!
end subroutine
