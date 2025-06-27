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
subroutine majsig(materf, se, seq, i1e, alpha, &
                  dp, plas, sig)
    implicit none
    real(kind=8) :: materf(5, 2), dp, se(*), seq, i1e, alpha, sig(6)
! =====================================================================
! --- MISE A JOUR DES CONTRAINTES -------------------------------------
! =====================================================================
    integer(kind=8) :: ii, ndt, ndi
    real(kind=8) :: young, nu, troisk, deuxmu, trois, deux, un
    real(kind=8) :: i1, dev(6), plas
    parameter(trois=3.0d0)
    parameter(deux=2.0d0)
    parameter(un=1.0d0)
! =====================================================================
    common/tdim/ndt, ndi
! =====================================================================
! --- AFFECTATION DES VARIABLES ---------------------------------------
! =====================================================================
    young = materf(1, 1)
    nu = materf(2, 1)
    troisk = young/(un-deux*nu)
    deuxmu = young/(un+nu)
! =====================================================================
! --- MISE A JOUR DU DEVIATEUR ----------------------------------------
! =====================================================================
    if (plas .eq. 2.0d0) then
        do ii = 1, ndt
            dev(ii) = 0.0d0
        end do
    else
        do ii = 1, ndt
            dev(ii) = se(ii)*(un-trois*deuxmu/deux*dp/seq)
        end do
    end if
!
! =====================================================================
! --- MISE A JOUR DU PREMIER INVARIANT --------------------------------
! =====================================================================
    i1 = i1e-trois*troisk*alpha*dp
! =====================================================================
! --- MISE A JOUR DU VECTEUR DE CONTRAINTES ---------------------------
! =====================================================================
    do ii = 1, ndt
        sig(ii) = dev(ii)
    end do
    do ii = 1, ndi
        sig(ii) = sig(ii)+i1/trois
    end do
! =====================================================================
end subroutine
