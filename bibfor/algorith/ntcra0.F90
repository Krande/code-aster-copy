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

subroutine ntcra0(sddisc)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmcrpx.h"
#include "asterfort/wkvect.h"
    character(len=19) :: sddisc
!
! ----------------------------------------------------------------------
!
! ROUTINE THER_* (STRUCTURES DE DONNES)
!
! CREATION SD ARCHIVAGE SPECIAL SI CALCUL NON TRANSITOIRE
!
! ----------------------------------------------------------------------
!
!
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
!
! ----------------------------------------------------------------------
!
    integer(kind=8), parameter :: iocc = 0
    character(len=16), parameter :: factorKeyword = ' ', keywStep = ' '
    character(len=1), parameter :: base = 'V'
    character(len=19) :: sdarch
    character(len=24) :: sdarchAinfJv, sdarchAexcJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    character(len=16), pointer :: sdarchAexc(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()

! - Name of datastructures to store
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    sdarchAexcJv = sdarch(1:19)//'.AEXC'

! - Create datastructures
    call wkvect(sdarchAexcJv, 'V V K16', 1, vk16=sdarchAexc)
    call wkvect(sdarchAinfJv, 'V V I', 4, vi=sdarchAinf)

! - Get parameters from ARCHIVAGE
    call nmcrpx(factorKeyword, keywStep, iocc, sdarch, base)
!
    call jedema()
!
end subroutine
