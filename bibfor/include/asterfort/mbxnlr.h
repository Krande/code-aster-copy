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
    subroutine mbxnlr(option,fami,nddl,nno,ncomp,kpg,ipoids,igeom,&
                  imate,ideplm,ideplp,ivectu,icontp,&
                  imatuu,dff,alpha,beta,&
                  vecteu,matric)
    character(len=16) :: option
    character(len=4) :: fami
    integer(kind=8) :: nddl, nno, ncomp
    integer(kind=8) :: kpg
    integer(kind=8) :: ipoids
    integer(kind=8) :: igeom, imate
    integer(kind=8) :: ideplm, ideplp
    integer(kind=8) :: ivectu, icontp, imatuu
    real(kind=8) :: dff(2, nno), alpha, beta
    aster_logical :: vecteu, matric
    end subroutine mbxnlr
end interface
