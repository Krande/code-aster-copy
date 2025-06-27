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
#include "LoadTypes_type.h"
!
interface
    subroutine char_crea_cart(phenom    , loadType     , load    , mesh, valeType,&
                              nbMap     , map          , nbCmp   ,&
                              createMap_, physQuantity_, cmpName_)
        character(len=*), intent(in) :: phenom
        character(len=16), intent(in) :: loadType
        character(len=8), intent(in) :: load, mesh
        character(len=4), intent(in) :: valeType
        integer(kind=8), intent(out) :: nbMap
        character(len=19), intent(out) :: map(LOAD_MAP_NBMAX)
        integer(kind=8), intent(out) :: nbCmp(LOAD_MAP_NBMAX)
        aster_logical, optional, intent(in) :: createMap_
        character(len=8), optional, intent(out) :: physQuantity_(LOAD_MAP_NBMAX)
        character(len=8), optional, intent(out) :: cmpName_(LOAD_MAP_NBMAX, LOAD_MAP_NBCMPMAX)
    end subroutine char_crea_cart
end interface
