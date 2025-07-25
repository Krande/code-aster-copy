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

subroutine derpro(a, da, b, db, dab)
!
    implicit none
    real(kind=8) :: a, da, b, db, dab
! --- BUT : CALCUL DE DAB = A*DB+B*DA ----------------------------------
! ======================================================================
! IN  : A      : VALEUR DE LA FONCTION A -------------------------------
! --- : DA     : VALEUR DE LA DERIVEE DE A -----------------------------
! --- : B      : VALEUR DE LA FONCTION B -------------------------------
! --- : DB     : VALEUR DE LA DERIVEE DE B -----------------------------
! OUT : DAB    : VALEUR DE LA DERIVEE DE A*B ---------------------------
! ======================================================================
    dab = a*db+b*da
! ======================================================================
end subroutine
