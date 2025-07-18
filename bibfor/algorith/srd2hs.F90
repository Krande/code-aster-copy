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

subroutine srd2hs(nmat, materf, devsig, sii, rcos3t, d2hds2)

!

!!!
!!! MODELE LKR : CALCUL DE LA DERIVEE 2NDE DE H PAR RAPPORT A DEVIATEUR SIGMA
!!!

! ===================================================================================
! IN  : NMAT           : DIMENSION TABLE DES PARAMETRES MATERIAU
!     : MATERF(NMAT,2) : PARAMETRES MATERIAU A T+DT
!     : DEVSIG(6)      : DEVIATEUR DES CONTRAINTES
!     : SII            : 2EME INVARIANT DU DEVIATEUR
!     : RCOS3T         : COS(3THETA) = SQRT(54)*DET(DEVSIG)/SII**3
!     : DHDS(6)        : DERIVEE DE H PAR RAPPORT A DEVSIG
! OUT : D2HDS2(6,6)    :  DERIVEE 2NDE H PAR RAPPORT A SIGMA (NDT X NDT)
! ===================================================================================

    implicit none

#include "asterc/r8pi.h"
#include "asterfort/cjst.h"
#include "asterfort/lcprte.h"
#include "asterfort/srd2de.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nmat
    real(kind=8) :: devsig(6), rcos3t, sii, d2hds2(6, 6), materf(nmat, 2)

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: ndi, ndt, i
    real(kind=8) :: mident(6, 6), gamma, beta, pdcds(6, 6), drdcos, sxs(6, 6)
    real(kind=8) :: dcosds(6)
    real(kind=8) :: ddetds(6), fact4, dcds1(6), dcds2(6)
    real(kind=8) :: d2dets(6, 6), r54, dets, fact3
    real(kind=8) :: pi, denom, fact1, fact2, d2rdc2
    real(kind=8) :: fact5(6, 6), fact6(6, 6), fact7(6, 6)
    real(kind=8) :: d2cds2(6, 6)
    common/tdim/ndt, ndi

    ddetds(:) = 0.d0
    dcds1(:) = 0.d0
    dcds2(:) = 0.d0
    pdcds(:, :) = 0.d0
    fact5(:, :) = 0.d0
    sxs(:, :) = 0.d0
    d2dets(:, :) = 0.d0
    mident(:, :) = 0.d0
    d2hds2(:, :) = 0.d0

    !!!
    !!! Recup. des para. mater
    !!!

    beta = materf(4, 2)
    gamma = materf(5, 2)
    pi = r8pi()
    r54 = sqrt(54.d0)

    !!!
    !!! Construction de d**2r/dcos(3t)**2
    !!!

    denom = 9.d0*(1.d0-(gamma*rcos3t)**2.d0)**(3.d0/2.d0)
    fact1 = -gamma**2.d0*sqrt(1.d0-(gamma*rcos3t)**2.d0)* &
            cos(beta*pi/6.d0-1.d0/3.d0*acos(gamma*rcos3t))
    fact2 = -3.d0*(gamma**3.d0)*rcos3t* &
            sin(beta*pi/6.d0-1.d0/3.d0*acos(gamma*rcos3t))
    d2rdc2 = (fact1+fact2)/denom

    !!!
    !!! Construction de d cos3t / ds x d cos3t / ds
    !!!

    dets = sii**3.d0*rcos3t/r54
    fact3 = r54/sii**3.d0
    fact4 = -3.d0*r54*dets/sii**5.d0
    call cjst(devsig, ddetds)
    dcds1(1:ndt) = fact3*ddetds(1:ndt)
    dcds2(1:ndt) = fact4*devsig(1:ndt)
    dcosds(1:ndt) = dcds1(1:ndt)+dcds2(1:ndt)
    call lcprte(dcosds, dcosds, pdcds)

    !!!
    !!! Construction de d(r)/d(cos3t)
    !!!

    drdcos = -gamma*sin(beta*pi/6.d0-1.d0/3.d0*acos(gamma*rcos3t))/ &
             3.d0/sqrt(1.d0-(gamma*rcos3t)**2.d0)

    !!!
    !!! Construction de d**2 cos3t / ds**2
    !!!

    call lcprte(devsig, ddetds, fact5)
    fact5(1:ndt, 1:ndt) = (-1.d0*2.d0)*fact5(1:ndt, 1:ndt)
    call lcprte(devsig, devsig, sxs)
    fact6(1:ndt, 1:ndt) = (5.d0*dets/sii**2.d0)*sxs(1:ndt, 1:ndt)
    call srd2de(devsig, d2dets)
    fact7(1:ndt, 1:ndt) = (sii**2.d0/3.d0)*d2dets(1:ndt, 1:ndt)

    do i = 1, ndt
        mident(i, i) = 1.d0
    end do

    d2cds2 = (3.d0*r54/sii**5.d0)*(fact5+fact6+fact7-dets*mident)

    !!!
    !!! Assemblage final
    !!!

    d2hds2 = (d2rdc2*pdcds)+(drdcos*d2cds2)

end subroutine
