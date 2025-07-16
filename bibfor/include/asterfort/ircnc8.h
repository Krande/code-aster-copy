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
    subroutine ircnc8(fileUnit , realFormat  , cplxFormat  ,&
                      nodeListNb, nodeListNume, nodeListName,&
                      lMeshCoor, meshDime    , meshCoor    ,&
                      cmpCataNb, cmpCataName ,&
                      cmpListNb, cmpListIndx ,&
                      nec      , nueq        ,&
                      prno     , codeInte    ,&
                      lmax     , lmin        ,&
                      lsup     , borsup      ,&
                      linf     , borinf      ,&
                      vale)
        integer(kind=8), intent(in) :: fileUnit
        character(len=8), intent(in) :: realFormat, cplxFormat
        integer(kind=8), intent(in) :: nodeListNb
        integer(kind=8), pointer :: nodeListNume(:)
        character(len=8), pointer :: nodeListName(:)
        aster_logical, intent(in) :: lMeshCoor
        integer(kind=8), intent(in) :: meshDime
        real(kind=8), pointer :: meshCoor(:)
        integer(kind=8), intent(in) :: cmpCataNb
        character(len=8), pointer :: cmpCataName(:)
        integer(kind=8), intent(in) :: cmpListNb
        integer(kind=8), pointer :: cmpListIndx(:)
        integer(kind=8), intent(in) :: nec
        integer(kind=8), pointer :: nueq(:), prno(:), codeInte(:)
        aster_logical, intent(in) :: lsup, linf, lmax, lmin
        real(kind=8),  intent(in) :: borsup, borinf
        complex(kind=8), pointer  :: vale(:)
    end subroutine ircnc8
end interface
