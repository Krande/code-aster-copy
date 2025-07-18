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
subroutine clffch(alias, type, nno, xi, yi, &
                  zi, xin, yin, zin, tn, &
                  ajx, ajy, ajz, bjxx, bjyy, &
                  bjzz, bjxy, bjxz, bjyz, ider)
!
!......................................................................C
!......................................................................C
!                                                                      C
! BUT:   CALCUL DES FONCTIONS DE FORMES ET DE LEURS DERIVEES           C
!        AU POINT DE COORDONNEES XI,YI,ZI                              C
!                                                                      C
! ATTENTION: CETTE ROUTINE EST SPECIFIQUE DES ELEMENTS HOMOGENEISE     C
!                            INI100                                    C
!                                                                      C
! ENTREES                                                              C
!      NNO         : NOMBRE DE NOEUDS                                  C
!      ALIAS       : NOM D'ALIAS DE L'ELEMENT :HEXA20 OU HEXA8         C
!      TYPE        : TYPE DE LA FONCTION DE FORME : POUTRE OU FLUIDE   C
!      XI,YI,ZI    : POINT DE CALCUL DES F FORMES ET DERIVEES          C
!      XIN,YIN,ZIN : COORDONNEES INTRINSEQUES                          C
!      IDER        : INDICATEUR DE CALCUL DES DERIVEES                 C
!                    IDER = 0 CALCUL DE TN                             C
!                    IDER = 1 CALCUL DE TN, AJ*                        C
!                    IDER = 2 CALCUL DE TN AJ* BJ*                     C
!                                                                      C
! SORTIES                                                              C
!      TN  : FONCTIONS DE FORMES EN XI,YI,ZI                           C
!      AJ (AJX, AJY, AJZ) : DERIVEES DES F FORMES EN XI,YI,ZI          C
!      BJ (BJXX ... BJYZ) : DERIVEES SECONDES DES F FORMES EN XI,YI,ZI C
!......................................................................C
!......................................................................C
!
    implicit none
#include "asterfort/utmess.h"
    character(len=6) :: alias, type
    real(kind=8) :: tn(*), ajx(*), ajy(*), ajz(*), xin(20), yin(20), zin(20)
    real(kind=8) :: bjxx(*), bjyy(*), bjzz(*), bjxy(*), bjxz(*), bjyz(*), xi, yi
    real(kind=8) :: zi
    integer(kind=8) :: ider, nno
!----------------------------------------------------------------------
    real(kind=8) :: x0, y0, z0, fxy, fxydx, fxydy, fxydxy, fz, fzdz, fzd2z, f
    real(kind=8) :: fdx, fdy, fdz
    integer(kind=8) :: i
!----------------------------------------------------------------------
!
    if (type .eq. 'POUTRE') then
!
!     -------------------------------------------------------------
!     --- F. DE FORME ASSOCIEES AUX DDLS DE FLEXION DES POUTRES ---
!     -------------------------------------------------------------
!
        do i = 1, nno
!
            x0 = xi*xin(i)
            y0 = yi*yin(i)
            z0 = zi*zin(i)
!
            fxy = (1.d0+x0)*(1.d0+y0)*0.25d0
            fxydx = xin(i)*(1.d0+y0)*0.25d0
            fxydy = (1.d0+x0)*yin(i)*0.25d0
            fxydxy = xin(i)*yin(i)*0.25d0
!
!          FONCTIONS DE FORME A DERIVEE NULLE AU BORD
!
            fz = -0.25d0*z0*z0*z0+0.75d0*z0+0.5d0
            fzdz = -0.75d0*z0*z0*zin(i)+0.75d0*zin(i)
            fzd2z = -1.5d0*z0*zin(i)*zin(i)
!
            tn(i) = fxy*fz
!
            if (ider .gt. 0) then
                ajx(i) = fxydx*fz
                ajy(i) = fxydy*fz
                ajz(i) = fxy*fzdz
            end if
!
            if (ider .gt. 1) then
                bjxx(i) = 0.d0
                bjyy(i) = 0.d0
                bjzz(i) = fxy*fzd2z
                bjxy(i) = fxydxy*fz
                bjxz(i) = fxydx*fzdz
                bjyz(i) = fxydy*fzdz
            end if
!
!          FONCTIONS DE FORME A VALEUR NULLE AU BORD
!
            fz = 0.25d0*(-zin(i)-zi+zin(i)*zi*zi+zi*zi*zi)
            fzdz = 0.25d0*(-1.d0+2.d0*zin(i)*zi+3.d0*zi*zi)
            fzd2z = 0.25d0*(2.d0*zin(i)+6.d0*zi)
!
            tn(i+8) = fxy*fz
!
            if (ider .gt. 0) then
                ajx(i+8) = fxydx*fz
                ajy(i+8) = fxydy*fz
                ajz(i+8) = fxy*fzdz
            end if
!
            if (ider .gt. 1) then
                bjxx(i+8) = 0.d0
                bjyy(i+8) = 0.d0
                bjzz(i+8) = fxy*fzd2z
                bjxy(i+8) = fxydxy*fz
                bjxz(i+8) = fxydx*fzdz
                bjyz(i+8) = fxydy*fzdz
            end if
        end do
!
    else if (type .eq. 'FLUIDE') then
!     --------------------------------------------------------------
        if (alias .eq. 'HEXA20') then
!
!     --------------------------------------------------------------
!     --- FONCTIONS DE FORMES AUX DDL FLUIDE DE LA MAILLE HEXA20 ---
!     --------------------------------------------------------------
!
            do i = 1, nno
                x0 = xi*xin(i)
                y0 = yi*yin(i)
                z0 = zi*zin(i)
!
                if (i .le. 8) then
                    f = 0.125d0*(1.d0+x0)*(1.d0+y0)*(1.d0+z0)*(x0+y0+z0-2.d0)
                    fdx = 0.125d0*(1.d0+y0)*(1.d0+z0)*(2.d0*xi+xin(i)*(z0+y0-1.d0))
                    fdy = 0.125d0*(1.d0+x0)*(1.d0+z0)*(2.d0*yi+yin(i)*(x0+z0-1.d0))
                    fdz = 0.125d0*(1.d0+x0)*(1.d0+y0)*(2.d0*zi+zin(i)*(x0+y0-1.d0))
!
                else if ((i .eq. 9) .or. (i .eq. 11) .or. (i .eq. 17) .or. ( &
                         i .eq. 19)) then
                    f = 0.25d0*(1.d0-xi*xi)*(1.d0+y0)*(1.d0+z0)
                    fdx = -0.5d0*xi*(1.d0+y0)*(1.d0+z0)
                    fdy = 0.25d0*(1.d0-xi*xi)*yin(i)*(1.d0+z0)
                    fdz = 0.25d0*(1.d0-xi*xi)*(1.d0+y0)*zin(i)
!
                else if ((i .eq. 10) .or. (i .eq. 12) .or. (i .eq. 18) .or. ( &
                         i .eq. 20)) then
                    f = 0.25d0*(1.d0-yi*yi)*(1.d0+x0)*(1.d0+z0)
                    fdy = -0.5d0*yi*(1.d0+x0)*(1.d0+z0)
                    fdx = 0.25d0*(1.d0-yi*yi)*xin(i)*(1.d0+z0)
                    fdz = 0.25d0*(1.d0-yi*yi)*(1.d0+x0)*zin(i)
!
                else if ((i .eq. 13) .or. (i .eq. 14) .or. (i .eq. 15) .or. ( &
                         i .eq. 16)) then
                    f = 0.25d0*(1.d0-zi*zi)*(1.d0+x0)*(1.d0+y0)
                    fdz = -0.5d0*zi*(1.d0+x0)*(1.d0+y0)
                    fdx = 0.25d0*(1.d0-zi*zi)*xin(i)*(1.d0+y0)
                    fdy = 0.25d0*(1.d0-zi*zi)*(1.d0+x0)*yin(i)
                end if
!
                tn(i) = f
                if (ider .gt. 0) then
                    ajx(i) = fdx
                    ajy(i) = fdy
                    ajz(i) = fdz
                end if
            end do
!
        else if (alias .eq. 'HEXA8 ') then
!
!     ----------------------------------------------------
!     --- FONCTIONS DE FORMES ASSOCIEES A LA GEOMETRIE ---
!     ---    ET AUX DDL FLUIDE DE LA MAILLE HEXA8      ---
!     ----------------------------------------------------
!
            do i = 1, nno
                x0 = xi*xin(i)
                y0 = yi*yin(i)
                z0 = zi*zin(i)
!
                tn(i) = (1.d0+x0)*(1.d0+y0)*(1.d0+z0)*0.125d0
!
                if (ider .gt. 0) then
                    ajx(i) = xin(i)*(1.d0+y0)*(1.d0+z0)*0.125d0
                    ajy(i) = yin(i)*(1.d0+x0)*(1.d0+z0)*0.125d0
                    ajz(i) = zin(i)*(1.d0+x0)*(1.d0+y0)*0.125d0
                end if
                if (ider .gt. 1) then
                    bjxx(i) = 0.d0
                    bjyy(i) = 0.d0
                    bjzz(i) = 0.d0
                    bjxy(i) = xin(i)*yin(i)*(1.d0+z0)*0.125d0
                    bjxz(i) = xin(i)*zin(i)*(1.d0+y0)*0.125d0
                    bjyz(i) = yin(i)*zin(i)*(1.d0+x0)*0.125d0
                end if
            end do
        end if
!     ---------------------------------------------------------
    else
        call utmess('F', 'ELEMENTS_20')
    end if
!
end subroutine
