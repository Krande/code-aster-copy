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
subroutine cmarge(N, M, Ned, Med, Ned_g, Med_g, margin,Nrd,Mrd,load_rate,c0_c,c0_crd)
       real(kind=8), pointer :: N(:)
       real(kind=8), pointer :: M(:)
       real(kind=8) :: Ned
       real(kind=8) :: Med
       real(kind=8) :: Ned_g
       real(kind=8) :: Med_g
       real(kind=8) :: margin
       real(kind=8) :: Nrd
       real(kind=8) :: Mrd
       real(kind=8) :: load_rate
       real(kind=8) :: c0_c
       real(kind=8) :: c0_crd
end subroutine cmarge
end interface
