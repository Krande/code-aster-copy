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
            subroutine trasst(modgen,numsst,isst1,lisint,nbeq1,nbmod,   &
     &nbint,solveu)
              character(len=8) :: modgen
              integer(kind=8) :: numsst
              integer(kind=8) :: isst1
              character(len=24) :: lisint
              integer(kind=8) :: nbeq1
              integer(kind=8) :: nbmod
              integer(kind=8) :: nbint
              character(len=19) :: solveu
            end subroutine trasst
          end interface 
