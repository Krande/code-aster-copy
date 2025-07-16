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

subroutine srds2h(nbmat, mater, s, dhds, ds2hds, retcom)

!

!!!
!!! MODELE LKR : CALCUL DES DERIVEES D(SII*H(THETA))/DSIG
!!!

! ===================================================================================
! IN  : NBMAT  : NOMBRE DE PARAMETRES DU MODELE
!     : MATER  : PARAMETRES DU MODELE
!     : INVAR  :  INVARIANT DES CONTRAINTES
!     : S      :  DEVIATEUR DES CONTRAINTES
!     : DHDS   : Dh(THETA)/DS
! OUT : DS2HDS : D(SII*H(THETA))/DSIG
! ===================================================================================

    implicit none

#include "asterc/r8miem.h"
#include "asterfort/cos3t.h"
#include "asterfort/srhtet.h"
#include "asterfort/r8inir.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nbmat, retcom
    real(kind=8) :: mater(nbmat, 2), s(6), dhds(6), ds2hds(6)

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: ndt, ndi, i, k
    real(kind=8) :: pref, r0c, rtheta
    real(kind=8) :: kron(6), iden6(6, 6)
    real(kind=8) :: a(6), b(6, 6), bt(6, 6)
    real(kind=8) :: sii, rcos3t, ptit
    common/tdim/ndt, ndi

    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/

    data iden6/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
               &0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
               &0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,&
               &0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,&
               &0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0,&
               &0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/

    !!!
    !!! Recuperation des parametres du modele
    !!!

    pref = mater(1, 2)

    !!!
    !!! Calcul du deviateur et verification qu'il n'est pas nul
    !!!

    retcom = 0
    ds2hds = 0.d0
    ptit = r8miem()

    sii = norm2(s(1:ndt))

    if (sii .lt. ptit) then
        retcom = 1
        goto 1000
    end if

    !!!
    !!! Recuperation de r(theta) et r0c
    !!!

    rcos3t = cos3t(s, pref, 1.0d-8)
    call srhtet(nbmat, mater, rcos3t, r0c, rtheta)

    !!!
    !!! Calcul du premier terme
    !!!

    call r8inir(6, 0.d0, a, 1)
    do i = 1, ndt
        a(i) = dhds(i)*sii+rtheta*s(i)/sii
    end do

    !!!
    !!! Calcul du second terme
    !!!

    call r8inir(6*6, 0.d0, b, 1)

    do i = 1, ndt
        do k = 1, ndt
            b(i, k) = iden6(i, k)-kron(i)*kron(k)/3.d0
        end do
    end do

    !!!
    !!! Resultat final
    !!!

    call r8inir(6, 0.d0, ds2hds, 1)
    bt(1:ndt, 1:ndt) = transpose(b(1:ndt, 1:ndt))
    ds2hds(1:ndt) = matmul(bt(1:ndt, 1:ndt), a(1:ndt))

1000 continue

end subroutine
