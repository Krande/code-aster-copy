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
            subroutine focrr2(nomfon,resu,base,nomcha,maille,noeud,cmp, &
     &npoint,nusp,ivari,nomvari,ier)
              character(len=19), intent(in) :: nomfon
              character(len=19), intent(in) :: resu
              character(len=1), intent(in) :: base
              character(len=16), intent(in) :: nomcha
              character(len=8), intent(in) :: maille
              character(len=8), intent(in) :: noeud
              character(len=8), intent(in) :: cmp
              integer(kind=8) :: npoint
              integer(kind=8) :: nusp
              integer(kind=8) :: ivari
              character(len=16), intent(in) :: nomvari
              integer(kind=8) :: ier
            end subroutine focrr2
          end interface 
