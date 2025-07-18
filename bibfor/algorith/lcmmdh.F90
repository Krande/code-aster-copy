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
subroutine lcmmdh(coeft, ifa, nmat, nbcomm, alphap, &
                  nfs, nsg, hsr, nbsys, is, &
                  nuecou, hs, soms1, soms2, soms3)
    implicit none
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/lcmmdc.h"
    integer(kind=8) :: ifa, nmat, nbcomm(nmat, 3), is, nbsys, nfs, nsg
    real(kind=8) :: coeft(*), alphap(12), hs, hsr(nsg, nsg), soms1, soms2, soms3
! person_in_charge: jean-michel.proix at edf.fr
! ======================================================================
!  CALCUL DE LA FONCTION H(OMEGA) POUR LA LOI D'ECOULEMENT  DD-CFC
!       IN  COEFT   :  PARAMETRES MATERIAU
!           IFA     :  NUMERO DE FAMILLE
!           NBCOMM  :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!           NMAT    :  NOMBRE DE MATERIAUX
!           ALPHAP  :  ALPHA =RHO*B**2 (TOTAL) A T+DT
!     OUT:
!           HS      :  FONCTION D'EVOLUTION DENSITE DISLOCATION
!           SOMS1    :  SOMME(j=1,12)(SQRT(a_sj omega_j))
!           SOMS2    :  SOMME(j=forest(s))(SQRT(a_sj) omega_j)
!           SOMS3    :  SOMME(j=copla(s))(SQRT(a_sj omega_j))
!     ----------------------------------------------------------------
    real(kind=8) :: a, b, y, termea, termeb, termey, denom, ceff, rmin, beta
    real(kind=8) :: numer
    real(kind=8) :: alphas, dcdals, unsurd, gc0, k
    integer(kind=8) :: iei, iu, iv, ifl, is3, iv3, nuecou
!     ----------------------------------------------------------------
!
!
    rmin = r8miem()
    ifl = nbcomm(ifa, 1)
    iei = nbcomm(ifa, 3)
!
!     LOI D'ECOULEMENT DD-CFC
!
    if ((nuecou .eq. 5) .or. (nuecou .eq. 8)) then
        a = coeft(ifl+3)
        b = coeft(ifl+4)
        y = coeft(ifl+6)
!
        beta = coeft(iei+2)
!         NUMHSR=NINT(COEFT(IEI+5))
!
!        EVOLUTION DE LA DENSITE DE DISLO
        termea = 0.d0
        denom = 0.d0
        numer = 0.d0
        do iu = 1, 12
!            PARTIE POSITIVE DE ALPHA
            if (alphap(iu) .gt. 0.d0) then
                denom = denom+sqrt(hsr(is, iu)*alphap(iu))
            end if
        end do
!        SOMME SUR FOREST(S)
        if (denom .gt. rmin) then
!           TERME AU NUMERATEUR SUR FOREST(S)
            termea = 0.d0
            do iv = 1, 12
                is3 = (is-1)/3
                iv3 = (iv-1)/3
                if (is3 .ne. iv3) then
!                 PARTIE POSITIVE DE ALPHA
                    if (alphap(iv) .gt. 0.d0) then
                        numer = numer+sqrt(hsr(is, iv))*alphap(iv)
                    end if
                end if
            end do
            termea = a*numer/denom
        end if
!
!        SOMME SUR COPLA(S)
        termeb = 0.d0
        if (nbsys .eq. 12) then
            do iv = 1, 12
                is3 = (is-1)/3
                iv3 = (iv-1)/3
!           PARTIE POSITIVE DE ALPHA
                if (is3 .eq. iv3) then
                    if (alphap(iv) .gt. 0.d0) then
                        termeb = termeb+sqrt(hsr(is, iv)*alphap(iv))
                    end if
                end if
            end do
        else if (nbsys .eq. 1) then
            alphas = alphap(is)
!           PARTIE POSITIVE DE ALPHA
            if (alphas .gt. 0.d0) then
                termeb = termeb+sqrt(hsr(is, is)*alphas)
            end if
        else
            ASSERT(.false.)
        end if
!
        call lcmmdc(coeft, ifa, nmat, nbcomm, alphap, &
                    is, ceff, dcdals)
!
!        TERME -Y*RHO_S
        if (alphap(is) .gt. 0.d0) then
            termey = -y*alphap(is)/beta
        else
            termey = 0.d0
        end if
        hs = (termea+termeb*b*ceff+termey)
        soms1 = denom
        soms2 = numer
        soms3 = termeb
    end if
!
!     LOI D'ECOULEMENT ECP-CFC
!
    if (nuecou .eq. 6) then
        beta = coeft(ifl+3)
        unsurd = coeft(ifl+4)
        gc0 = coeft(ifl+6)
        k = coeft(ifl+7)
!
!        NUMHSR=NINT(COEFT(IEI+2))
!
!       EVOLUTION DE LA DENSITE DE DISLO
!
        denom = 0.d0
        do iu = 1, 12
            if ((iu .ne. is) .and. (alphap(iu) .gt. 0.d0)) then
                denom = denom+alphap(iu)
            end if
        end do
        denom = sqrt(denom)
!
        hs = beta*unsurd+denom/k
!
        if (alphap(is) .gt. 0.d0) then
            hs = hs-gc0*alphap(is)/beta
        end if
!
        soms1 = 0.d0
        soms2 = 0.d0
        soms3 = 0.d0
    end if
!
end subroutine
