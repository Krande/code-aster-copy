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
#include "asterf_types.h"
!
interface 
    subroutine xhmddl(ndim, nfh, ddls, nddl, nno, nnos,&
                      stano, matsym, option, nomte, mat,&
                      vect, ddlm, nfiss, jfisno, lcontx, contac)
        integer(kind=8) :: ndim
        integer(kind=8) :: nfh
        integer(kind=8) :: ddls
        integer(kind=8) :: nddl
        integer(kind=8) :: nno
        integer(kind=8) :: nnos
        integer(kind=8) :: stano(*)
        aster_logical :: matsym
        character(len=16) :: option
        character(len=16) :: nomte
        real(kind=8) :: mat(*)
        real(kind=8) :: vect(*)
        integer(kind=8) :: ddlm
        integer(kind=8) :: nfiss
        integer(kind=8) :: jfisno
        aster_logical :: lcontx
        integer(kind=8) :: contac
    end subroutine xhmddl
end interface 
