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
    subroutine irecri(fileUnit   , dsName        , lResu     ,&
                      titleKeywf , titleKeywfIocc,&
                      storeNb    , storeListIndx ,&
                      fieldListNb, fieldListType , &
                      paraNb     , paraName      , paraFormat,&
                      cmpUserNb  , cmpUserName   ,&
                      cellUserNb , cellUserNume  ,&
                      nodeUserNb , nodeUserNume  ,&
                      lMeshCoor  , lmax          , lmin,&
                      lsup       , borsup        ,&
                      linf       , borinf        ,&
                      realFormat , cplxFormat)
        integer(kind=8), intent(in) :: fileUnit
        character(len=*), intent(in) :: dsName, titleKeywf
        integer(kind=8), intent(in) :: titleKeywfIocc
        aster_logical, intent(in) :: lResu
        integer(kind=8), intent(in) :: storeNb
        integer(kind=8) , pointer :: storeListIndx(:)
        integer(kind=8), intent(in) :: fieldListNb
        character(len=*), pointer :: fieldListType(:)
        integer(kind=8), intent(in) :: paraNb
        character(len=*), pointer :: paraName(:)
        character(len=1), intent(in) :: paraFormat
        integer(kind=8), intent(in) :: cmpUserNb
        character(len=8), pointer :: cmpUserName(:)
        integer(kind=8), intent(in) :: nodeUserNb
        integer(kind=8) , pointer :: nodeUserNume(:)
        integer(kind=8), intent(in) :: cellUserNb
        integer(kind=8) , pointer :: cellUserNume(:)
        aster_logical, intent(in) :: lMeshCoor
        aster_logical, intent(in) :: lsup, linf, lmax, lmin
        real(kind=8), intent(in) :: borsup, borinf
        character(len=*), intent(in) :: realFormat, cplxFormat
    end subroutine irecri
end interface
