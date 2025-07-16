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

subroutine fun2(xi1, xi2, pp, xkk, qq, vt, n)
    implicit none
    integer(kind=8) :: n
    real(kind=8) :: xi1, xi2, pp, xkk, qq, vt
!
! --------------------------------------------------------------------------------------------------
!
!   calcule le moment d'inertie équivalent d'une poutre droite a
!   section variable sous l’hypothèse de variation linéaire des coordonnées
!
! --------------------------------------------------------------------------------------------------
!
!    TYPE  !   NOM  !              SIGNIFICATION
!    --------------------------------------------------------
! IN R*8   !  XI1   !  moment d'inertie section initiale
! IN R*8   !  XI2   !  moment d'inertie section finale
! IN R*8   !  PP    !  coeff. de cisaillement dépendant du matériau (p = e/(g*as*l*l) )
! IN  I    !   N    !  type de calcul  (n =1  3 ou 4) degré du polynôme du dénominateur
! OUT R*8  !  XKK   !  l’élément diagonal a e/l3 pres
! OUT R*8  !   Q    !
! OUT R*8  !   VT   !
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: alg, di, phi, a3, rho
    real(kind=8) :: xx, xx25, xx50, xx75, pxi1
    real(kind=8) :: q_lim, xkk_lim, vt_lim
    real(kind=8) :: uu, dd
!
! --------------------------------------------------------------------------------------------------
    !
    !   Réécriture des expressions pour traiter le cas xi1 = xi2 + epsi (epsi petit)
    !
    ! Pour :
    !   n=1 : passage par un développement limité, à cause du log()
    !           pour q, vt, xkk
    !   n=3 : passage par un développement limité, à cause du log()
    !           pour xkk
    !   n=4 : on peut réécrire et simplifier. Pas besoin de développement limité
    !
    !
    !   n = 1
    !       alg = log(xi2/xi1); di= xi2-xi1; phi= p*di
    !       q   = 1.0/alg - xi1/di
    !       xkk = di/((xi2+xi1)/(2.0*di)-1.0/alg+phi)
    !       vt  = di/alg
    !   n = 3
    !       q   = xi1**(1.0/3.0)/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
    !       vt  = 2.0*((xi1*xi2)**(2.0/3.0))/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
    !       !
    !       a   = (xi2/xi1)**(1.0/3.0); a3=(a-1.0)**3; phi=p*xi1*a3
    !       rho = log(a)-2.0*q*(a-1.0)+phi
    !       xkk = xi1*a3/rho
    !   n = 4
    !       a   = (xi2/xi1)**0.250; a3=a*a*a; a2=(a-1.0)**2;
    !       phi = p*xi1*a2; c = a3-3.0*a+2.0
    !       q   = c/(2.0*(a-1.0)*(a3-1.0))
    !       rho = a2/(3.0*a3)-c*q/(6.0*a3)+phi
    !       xkk = a2*xi1/rho
    !       vt  = 3.0*xi1*a3/(a*a + a + 1.0)
    !
    if (n .ge. 4) then
        xx = xi2/xi1
        pxi1 = pp*xi1
        xx25 = xx**0.25; xx50 = xx**0.50; xx75 = xx**0.75
        qq = 0.50*(xx25+2.0)/(xx50+xx25+1.0)
        xkk = 4.0*xi1*(xx75+xx50+xx25)/ &
              (4.0*pxi1*(xx75+xx50+xx25)+1.0)
        vt = 3.0*xi1*xx75/(xx50+xx25+1.0)
    else
        !
        pxi1 = pp*xi1
        q_lim = 0.50
        xkk_lim = 12.0*xi1/(12.0*pxi1+1.0)
        vt_lim = xi1
        !
        dd = (xi2-xi1)/(xi2+xi1)
        if (abs(dd) <= 1.0E-10) then
            ! On considère que xi2=xi1
            qq = q_lim
            xkk = xkk_lim
            vt = vt_lim
        else if (abs(dd) <= 1.0E-06) then
            ! Développement limité au 3ème ordre
            uu = dd/(12.0*pxi1+1.0)
            if (n .le. 1) then
                qq = q_lim*(1.0-dd/3.0-4.0*(dd**3)/45.0)
                xkk = xkk_lim*(1.0+uu-(48.0*pxi1-11.0)*(uu**2)/15.0+ &
                               (576.0*(pxi1**2)+11.0)*(uu**3)/15.0)
                vt = vt_lim*(1.0+dd+2.0*(dd**2)/3.0+2.0*(dd**3)/3.0)
            else if (n .eq. 3) then
                qq = xi1**(1.0/3.0)/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
                vt = 2.0*((xi1*xi2)**(2.0/3.0))/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
                !
                xkk = xkk_lim*(1.0+uu-(24.0*pxi1-3.0)*(uu**2)/5.0+ &
                               (288.0*(pxi1**2)+3.0)*(uu**3)/5.0)
            end if
        else
            if (n .le. 1) then
                alg = log(xi2/xi1)
                di = xi2-xi1
                phi = pp*di
                qq = 1.0/alg-xi1/di
                xkk = di/((xi2+xi1)/(2.0*di)-1.0/alg+phi)
                vt = di/alg
            else if (n .eq. 3) then
                qq = xi1**(1.0/3.0)/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
                vt = 2.0*((xi1*xi2)**(2.0/3.0))/(xi2**(1.0/3.0)+xi1**(1.0/3.0))
                !
                xx = (xi2/xi1)**(1.0/3.0)
                a3 = (xx-1.0)**3
                phi = pp*xi1*a3
                rho = log(xx)-2.0*qq*(xx-1.0)+phi
                xkk = xi1*a3/rho
            end if
        end if
    end if
end subroutine
