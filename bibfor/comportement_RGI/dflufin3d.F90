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

subroutine dflufin3d(sige6, bw, pw, bg, pg, dsw6, delta, rc, &
                     xflu, dfin, cmp1, dfmx2)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================

!   endommagement par fluage
    implicit none
#include "rgi_module.h"

    real(kind=8), intent(in) :: sige6(6), bw, pw, bg, pg, dsw6(6)
    real(kind=8), intent(in) :: delta, rc, xflu
    real(kind=8), intent(out) :: dfin, cmp1
    real(kind=8), intent(in) :: dfmx2

    real(kind=8) :: tauflu1, sigs, sig0, sigd6(6), sigeq, sigs3, taulim
    real(kind=8) :: taueq, tauflu0, xdenom, dfmx
    integer(kind=8) i

!     initialisation
    taueq = 0.d0
    dfmx = 0.d0
!***********************************************************************
!     la non linearite est estimee avec le deviateur de la contrainte
!     macroscopique
!***********************************************************************

!     calcul du deviateur des contraintes macroscopiques
    sigs = 0.d0
!     deduction des pressions intra poreuse
    sig0 = bw*pw+bg*pg
!     deduction de la surpression par fluage de dessiccation (dsw6=s/sfld*d(bwpw))
    do i = 1, 3
        sigs = sigs+sige6(i)+dsw6(i)-sig0
    end do
    sigeq = 0.d0
    sigs3 = sigs/3.d0
    do i = 1, 3
        sigd6(i) = (sige6(i)+dsw6(i)-sig0)-sigs3
        sigeq = sigeq+sigd6(i)**2.d0
    end do
    do i = 4, 6
        sigeq = sigeq+2.d0*((sige6(i)+dsw6(i))**2.d0)
    end do
    sigeq = dsqrt(sigeq/2.d0)
    taulim = rc*(1.d0/dsqrt(3.d0)-delta/3.d0)
    taueq = dmin1(sigeq+delta*sigs3, taulim)
    tauflu0 = (taueq/taulim)
    tauflu1 = dmax1(tauflu0, 0.d0)

!     calcul du coeff de fluage non lineaire
    if (xflu .gt. 1.d0) then
!         le 1.5 suivant provient du taux de charge caracteristique pour
!         caracteriser l amplification non lineaire
        dfmx = 1.5d0*(1.d0-1.d0/xflu)
        if (taueq .gt. 0.d0) then
            xdenom = taulim-dfmx*taueq
            if (xdenom .gt. (taulim/XMAX)) then
                cmp1 = taulim/xdenom
!               calcul de l'endommagement asymptotique avec le taux de charge actuel
                dfin = dmin1(tauflu1*dfmx2, 1.d0-EPSIL)
            else
                cmp1 = XMAX
                dfin = dfmx2
            end if
        else
            cmp1 = 1.d0
            dfin = 0.d0
        end if
    else
        cmp1 = 1.d0
        dfin = 0.d0
    end if
end subroutine
