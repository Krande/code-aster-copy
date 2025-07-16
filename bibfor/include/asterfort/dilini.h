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
#include "asterf_types.h"
!
interface
    subroutine dilini(option, ivf, ivf2, idfde,&
                      idfde2, jgano, ndim, ipoids, ipoid2,&
                      npi, dimdef, nddls, nddlm, lgpg,&
                      dimcon, typmod, dimuel, nnom, nnos, ds_dil)
        use dil_type
        character(len=16) :: option
        integer(kind=8) :: ivf
        integer(kind=8) :: ivf2
        integer(kind=8) :: idfde
        integer(kind=8) :: idfde2
        integer(kind=8) :: jgano
        integer(kind=8) :: ndim
        integer(kind=8) :: ipoids
        integer(kind=8) :: ipoid2
        integer(kind=8) :: npi
        integer(kind=8) :: dimdef
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: lgpg
        integer(kind=8) :: dimcon
        character(len=8) :: typmod(2)
        integer(kind=8) :: dimuel
        integer(kind=8) :: nnom
        integer(kind=8) :: nnos
        type(dil_modelisation) :: ds_dil
    end subroutine dilini
end interface
