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
#include "asterf_types.h"
!
interface
    subroutine xmulhm(contac, ddls, ddlc, ddlm, jaint, ifiss,&
                      jheano, vstnc, lact, lcalel, lelim,&
                      nfh, nfiss, ninter,&
                      nlact, nno, nnol, nnom, nnos,&
                      pla, pos, typma, jstano)
        integer(kind=8) :: contac
        integer(kind=8) :: ddls
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: jaint
        integer(kind=8) :: ifiss
        integer(kind=8) :: jheano
        integer(kind=8) :: vstnc(*)
        integer(kind=8) :: lact(16)
        aster_logical :: lcalel
        aster_logical :: lelim
        integer(kind=8) :: nfh
        integer(kind=8) :: nfiss
        integer(kind=8) :: ninter
        integer(kind=8) :: nlact(2)
        integer(kind=8) :: nno
        integer(kind=8) :: nnol
        integer(kind=8) :: nnom
        integer(kind=8) :: nnos
        integer(kind=8) :: pla(27)
        integer(kind=8) :: pos(16)
        character(len=8) :: typma
        integer(kind=8), optional, intent(in) :: jstano
    end subroutine xmulhm
end interface
