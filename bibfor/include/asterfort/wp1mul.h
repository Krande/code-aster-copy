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
            subroutine wp1mul(lmasse,lamor,lraide,ptorig,tolf,nitf,     &
     &nbfreq,mxresf,nprec,resufi,resufr,solveu)
              integer(kind=8) :: mxresf
              integer(kind=8) :: lmasse
              integer(kind=8) :: lamor
              integer(kind=8) :: lraide
              complex(kind=8) :: ptorig(3,*)
              real(kind=8) :: tolf
              integer(kind=8) :: nitf
              integer(kind=8) :: nbfreq
              integer(kind=8) :: nprec
              integer(kind=8) :: resufi(mxresf,*)
              real(kind=8) :: resufr(mxresf,*)
              character(len=19) :: solveu
            end subroutine wp1mul
          end interface 
