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

subroutine xderfk_wrap(kappa, mu, r, theta, ndim, dfkdpo, option, istano)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "asterfort/xderfk.h"
#include "asterc/r8depi.h"
#include "asterc/r8prem.h"
!
    integer(kind=8) :: ndim, istano
    real(kind=8) :: r, theta, dfkdpo(ndim, ndim, 2), kappa, mu
    character(len=*) :: option
!
!
!     BUT : DERIVEES DES FONCTIONS D'ENRICHISSEMENT
!           DANS LA BASE POLAIRE (R,THETA)
!
! IN  R      : PREMIERE COORDONNEE DANS LA BASE POLAIRE
! IN  THETA  : SECONDE COORDONNEE DANS LA BASE POLAIRE
! OUT DFKDPO : DERIVEES DES FONCTIONS D'ENRICHISSEMENT
!   -- FORMAT DE STOCKAGE DES DERIVEES --
!       DFKDPO(i, j, l)
!            i <-> Ki
!            j <-> Kij=Ki.ej
!            l <-> [dKij/dr dKij/dtheta]
!----------------------------------------------------------------
!
    character(len=8) :: pref
!
    pref = option
!    if (istano.eq.-2) pref='SMOOTH'
!
    call xderfk(kappa, mu, r, theta, ndim, dfkdpo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (pref .eq. 'DEFAULT') then
        goto 999
!
    elseif (pref .eq. 'BASIC') then
        dfkdpo(1:ndim, 1:ndim, 1:2) = 0.
        dfkdpo(1, 1, 1:2) = [1/sqrt(r)*cos(theta/2.d0), -0.5*sqrt(r)*sin(theta/2.d0)]
        dfkdpo(1, 2, 1:2) = [1/sqrt(r)*sin(theta/2.d0), 0.5*sqrt(r)*cos(theta/2.d0)]
        dfkdpo(2, 1, 1:2) = [1/sqrt(r)*sin(theta/2.d0), 0.5*sqrt(r)*cos(theta/2.d0)]
        dfkdpo(2, 2, 1:2) = [1/sqrt(r)*cos(theta/2.d0), -0.5*sqrt(r)*sin(theta/2.d0)]
        if (ndim .eq. 3) then
            dfkdpo(3, 3, 1:2) = [1/sqrt(r)*sin(theta/2.d0), 0.5*sqrt(r)*cos(theta/2.d0)]
        end if
!
    elseif (pref .eq. 'SMOOTH') then
        dfkdpo(1:ndim, 1:ndim, 1) = 0.
        if (r .gt. r8prem()) then
            dfkdpo(1:ndim, 1:ndim, 2) = dfkdpo(1:ndim, 1:ndim, 2)/sqrt(r)
        end if
!
    elseif (pref .eq. 'JUMP') then
        dfkdpo = 2.d0*mu*sqrt(r8depi())*dfkdpo/(kappa+1)
        !dfkdpo(1:2,1:2,1:2)=2.d0*mu*sqrt(r8depi())*dfkdpo(1:2,1:2,1:2)/(kappa+1)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
999 continue
!
end subroutine
