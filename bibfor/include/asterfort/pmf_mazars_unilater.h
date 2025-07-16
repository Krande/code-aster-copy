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
!
interface
    subroutine pmf_mazars_unilater(for_pmf, nf,    nbvalc, &
                                   compor,  crit,  defam, defap, varim, &
                                   varimp,  contm, defm,  ddefp, modf,  &
                                   sigf,    varip, codret)
        use pmfcom_type
        type(pmfcom_user), intent(in) :: for_pmf
        !
        integer(kind=8) :: nf
        integer(kind=8) :: nbvalc
        character(len=24) :: compor(*)
        real(kind=8) :: crit(*)
        real(kind=8) :: defam(*)
        real(kind=8) :: defap(*)
        real(kind=8) :: varim(nbvalc*nf)
        real(kind=8) :: varimp(nbvalc*nf)
        real(kind=8) :: contm(nf)
        real(kind=8) :: defm(nf)
        real(kind=8) :: ddefp(nf)
        real(kind=8) :: modf(nf)
        real(kind=8) :: sigf(nf)
        real(kind=8) :: varip(nbvalc*nf)
        integer(kind=8) :: codret
    end subroutine pmf_mazars_unilater
end interface
