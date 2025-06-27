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

subroutine cnocns_wrap(cnoz, basez, cnsz)
! A_UTIL
    implicit none
#include "asterf_types.h"
#include "asterfort/cnocns.h"
!
    character(len=*) :: cnoz, basez, cnsz
!
! --------------------------------------------------------------------------------------------------
!
! BUT : TRANSFORMER UN CHAM_NO (CNOZ) EN CHAM_NO_S (CNSZ)
!
! --------------------------------------------------------------------------------------------------
!
!     ARGUMENTS:
! CNOZ    IN/JXIN  K19 : SD CHAM_NO A TRANSFORMER
! BASEZ   IN       K1  : BASE DE CREATION POUR CNSZ : G/V/L
! CNSZ    IN/JXOUT K19 : SD CHAM_NO_S A CREER
!
! --------------------------------------------------------------------------------------------------
!
    call cnocns(cnoz, basez, cnsz, ASTER_TRUE)
end subroutine
