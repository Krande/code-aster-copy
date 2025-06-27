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
    subroutine xlmail(fiss, nmaen1, nmaen2, nmaen3, nmafon,&
                      jmaen1, jmaen2, jmaen3, jmafon, nfon,&
                      jfon, jnofaf, nbfond, jbas, jtail, jfonmu,&
                      ndim, goinop)
        character(len=8) :: fiss
        integer(kind=8) :: nmaen1
        integer(kind=8) :: nmaen2
        integer(kind=8) :: nmaen3
        integer(kind=8) :: nmafon
        integer(kind=8) :: jmaen1
        integer(kind=8) :: jmaen2
        integer(kind=8) :: jmaen3
        integer(kind=8) :: jmafon
        integer(kind=8) :: nfon
        integer(kind=8) :: jfon
        integer(kind=8) :: jnofaf
        integer(kind=8) :: nbfond
        integer(kind=8) :: jbas
        integer(kind=8) :: jtail
        integer(kind=8) :: jfonmu
        integer(kind=8) :: ndim
        aster_logical :: goinop
    end subroutine xlmail
end interface
