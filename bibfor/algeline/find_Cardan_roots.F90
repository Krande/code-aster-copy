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

subroutine find_Cardan_roots(coef, roots, nbroot)
!
    implicit none
#include "asterc/r8pi.h"
#include "asterfort/cubic_root.h"
!
    real(kind=8), intent(in) :: coef(4)
    real(kind=8), intent(out) :: roots(3)
    integer(kind=8), intent(out) :: nbroot
!
!
!    RACINES REELLES POLYNOME DEGRE 3 : A*X**3 + B*X**2 + C*X + D = 0
!    PAR LES FORMULES DE CARDAN
!
! IN  COEF : COEFFICIENTS EN COMMENCANT PAR LE PLUS HAUT DEGRE
! OUT ROOTS : RACINES
! OUT NBROOT : NOMBRE DE RACINES
!
    real(kind=8) :: a, b, c, d, p, q, delta
    real(kind=8) :: u1, v1, z1, z2, z3, prec
    real(kind=8) :: pi, acos1, argument, racine, ref_u, arg_1, arg_2
    real(kind=8) :: temp1, temp2, temp3, ref_p, ref_q, ref_delta
!
    prec = 1.d-13
    pi = r8pi()
!
    a = coef(1)
    b = coef(2)
    c = coef(3)
    d = coef(4)
!
    temp1 = -b*b/(3.d0*a*a)
    temp2 = c/a
    p = temp1+temp2
!
!     doit-on considerer que p=0?
    ref_p = max(abs(temp1), abs(temp2))
    if (abs(p) .lt. prec*ref_p) then
        p = 0.d0
    end if
!
    temp1 = 2.d0*b*b*b/(27.d0*a*a*a)
    temp2 = -9.d0*b*c/(27.d0*a*a)
    temp3 = d/a
    q = temp1+temp2+temp3
!
!     doit-on considerer que q=0?
    ref_q = max(abs(temp1), abs(temp2), abs(temp3))
    if (abs(q) .lt. prec*ref_q) then
        q = 0.d0
    end if
!
!     discriminant
    delta = -(4.d0*p*p*p+27.d0*q*q)
!     doit-on considerer que delta = 0?
!     reference pour delta
    ref_delta = max(abs(4.d0*p*p*p), 27.d0*q*q)
!
    if (delta .gt. prec*ref_delta) then
        ! trois reelles distinctes
        nbroot = 3
        argument = -0.5d0*q*sqrt(27.d0/(-p*p*p))
        if (abs(argument-1.d0) .lt. prec) then
            argument = 1.d0
        else if (abs(argument+1.d0) .lt. prec) then
            argument = -1.d0
        end if
        acos1 = 1.d0/3.d0*acos(argument)
        z1 = 2*sqrt(-p/3.d0)*cos(acos1)
        roots(1) = z1-b/(3.d0*a)
        z2 = 2*sqrt(-p/3.d0)*cos(acos1+2*pi/3.d0)
        roots(2) = z2-b/(3.d0*a)
        z3 = 2*sqrt(-p/3.d0)*cos(acos1-2*pi/3.d0)
        roots(3) = z3-b/(3.d0*a)
!
    else if (abs(delta) .lt. prec*ref_delta) then
        !
        ! on considere que delta = 0
        if (p .eq. 0.d0 .and. q .eq. 0.d0) then
            ! 0 est racine triple
            nbroot = 1
            roots(1) = 0.d0
        else
            ! deux racines reelles dont une double
            nbroot = 2
            ! la racine double
            if (p .ne. 0.d0) then
                !
                z1 = -3.d0*q/(2.d0*p)
            else
                ! securite, ms ne doit pas arriver normalement!!
                z1 = -cubic_root(-q/2.d0)
            end if
            roots(1) = z1-b/(3.d0*a)
            ! la racine simple
            ! solution plus simple si p different de 0
            if (p .ne. 0.d0) then
                z2 = 3.d0*q/p
            else
                ! securite, ms ne doit pas arriver normalement!!
                z2 = 2.d0*cubic_root(-q/2.d0)
            end if
            roots(2) = z2-b/(3.d0*a)
        end if
!
    else
        ! une seule racine reelle
        nbroot = 1
        racine = sqrt(-delta/27.d0)
        ref_u = max(abs(q), abs(racine))
        arg_1 = -q+racine
        if (abs(arg_1) .lt. prec*ref_u) then
            arg_1 = 0.d0
        end if
        arg_2 = -q-racine
        if (abs(arg_2) .lt. prec*ref_u) then
            arg_2 = 0.d0
        end if
        u1 = cubic_root(0.5d0*arg_1)
        v1 = cubic_root(0.5d0*arg_2)
        z1 = u1+v1
        roots(1) = z1-b/(3.d0*a)
    end if
!
end subroutine
