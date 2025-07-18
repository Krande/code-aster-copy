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
subroutine lkds2h(nbmat, mater, invar, s, dhds, &
                  ds2hds, retcom)
!
    implicit none
#include "asterc/r8miem.h"
#include "asterfort/cos3t.h"
#include "asterfort/lkhtet.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nbmat, retcom
    real(kind=8) :: invar, mater(nbmat, 2), s(6), dhds(6), ds2hds(6)
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE--------------------------
! =================================================================
! --- BUT : CALCUL DES DERICEES d(sII*h(THETA))/dsig --------------
! =================================================================
! IN  : NBMAT  : NOMBRE DE PARAMETRES DU MODELE -------------------
! --- : MATER  : PARAMETRES DU MODELE -----------------------------
! --- : INVAR :  INVARIANT DES CONTRAINTES ------------------------
! --- : S     :  DEVIATEUR DES CONTRAINTES ------------------------
!     : DHDS   : dh(THETA)/ds--------------------------------------
! OUT : DS2HDS: d(sII*h(THETA))/dsig-------------------------------
! =================================================================
    integer(kind=8) :: ndt, ndi, i, k
    real(kind=8) :: h0ext, pref, h0e, h0c, htheta
    real(kind=8) :: kron(6), iden6(6, 6)
    real(kind=8) :: a(6), b(6, 6), bt(6, 6)
    real(kind=8) :: sii, rcos3t
    real(kind=8) :: zero, un, trois, lgleps, ptit
    real(kind=8) :: fact1
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(zero=0.0d0)
    parameter(un=1.0d0)
    parameter(trois=3.0d0)
    parameter(lgleps=1.0d-8)
! -----------------------------------------------------------------
    common/tdim/ndt, ndi
! -----------------------------------------------------------------
    data kron/un, un, un, zero, zero, zero/
    data iden6/un, zero, zero, zero, zero, zero,&
     &                 zero, un, zero, zero, zero, zero,&
     &                 zero, zero, un, zero, zero, zero,&
     &                 zero, zero, zero, un, zero, zero,&
     &                 zero, zero, zero, zero, un, zero,&
     &                 zero, zero, zero, zero, zero, un/
! =================================================================
! --- RECUPERATION DE PARAMETRES DU MODELE ------------------------
! =================================================================
    h0ext = mater(4, 2)
    pref = mater(1, 2)
! =================================================================
! --- CALCUL DU DEVIATEUR ET DU PREMIER INVARIANT DES CONTRAINTES -
! =================================================================
    retcom = 0
    ds2hds = 0.d0
    ptit = r8miem()
    sii = norm2(s(1:ndt))
    if (sii .lt. ptit) then
        retcom = 1
        goto 1000
    end if
!
! =================================================================
! --- CALCUL DE h(THETA), H0E ET H0C, -----------------------------
! =================================================================
    rcos3t = cos3t(s, pref, lgleps)
    call lkhtet(nbmat, mater, rcos3t, h0e, h0c, &
                htheta)
!
    fact1 = (h0c-h0ext)/(h0c-h0e)
! =================================================================
! --- CALCUL DU PREMIER TERME
! =================================================================
!
    call r8inir(6, 0.d0, a, 1)
    do i = 1, ndt
        a(i) = fact1*dhds(i)*sii+htheta*s(i)/sii
    end do
!
! =================================================================
! --- CALCUL DU SECOND  TERME
! =================================================================
    call r8inir(6*6, 0.d0, b, 1)
    do i = 1, ndt
        do k = 1, ndt
            b(i, k) = iden6(i, k)-kron(i)*kron(k)/trois
        end do
    end do
!
! =================================================================
! --- RESULTAT FINAL
! =================================================================
    call r8inir(6, 0.d0, ds2hds, 1)
!
    bt(1:ndt, 1:ndt) = transpose(b(1:ndt, 1:ndt))
    ds2hds(1:ndt) = matmul(bt(1:ndt, 1:ndt), a(1:ndt))
!
! =================================================================
1000 continue
end subroutine
