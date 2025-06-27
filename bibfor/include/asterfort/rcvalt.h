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
            subroutine rcvalt(fami,kpg,ksp,poum,jmat,nomat,mfact,nbpar, &
     &nompar,valpar,nbres,valres,icodre,iarret)
              integer(kind=8), intent(in) :: nbpar
              character(len=*), intent(in) :: fami
              integer(kind=8), intent(in) :: kpg
              integer(kind=8), intent(in) :: ksp
              character(len=1), intent(in) :: poum
              integer(kind=8), intent(in) :: jmat
              character(len=*), intent(in) :: nomat
              character(len=*), intent(in) :: mfact
              character(len=*), intent(in) :: nompar(nbpar)
              real(kind=8), intent(in) :: valpar(nbpar)
              integer(kind=8), intent(in) :: nbres
              real(kind=8), intent(out) :: valres(*)
              integer(kind=8), intent(out) :: icodre(*)
              integer(kind=8), intent(in) :: iarret
            end subroutine rcvalt
          end interface 
