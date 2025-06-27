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
            subroutine noligr(ligrz,igrel,numel,nunoeu,code,inema, &
                       &nbno,jlgns,rapide,jliel0,jlielc,jnema0,jnemac,l_lag1)
              character(len=*), intent(in) :: ligrz
              integer(kind=8), intent(in) :: igrel
              integer(kind=8), intent(in) :: numel
              integer(kind=8), intent(in) :: nunoeu
              integer(kind=8), intent(in) :: code
              integer(kind=8), intent(inout) :: inema
              integer(kind=8), intent(inout) :: nbno
              integer(kind=8), intent(in) :: jlgns
              character(len=3) ,optional, intent(in) :: rapide
              integer(kind=8) ,optional, intent(in) :: jliel0
              integer(kind=8) ,optional, intent(in) :: jlielc
              integer(kind=8) ,optional, intent(in) :: jnema0
              integer(kind=8) ,optional, intent(in) :: jnemac
              aster_logical, intent(in), optional :: l_lag1
            end subroutine noligr
          end interface 
