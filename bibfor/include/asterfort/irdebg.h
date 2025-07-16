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
    subroutine irdebg(dsName,&
                      fileFormat, fileUnit, fileVersion,&
                      lResu, lMesh, lField,&
                      lMeshCoor,  paraFormat,realFormat,cplxFormat,&
                      lsup, linf, lmax, lmin,&
                      borsup, borinf,&
                      storeListNb, storeListIndx,&
                      fieldListNb, fieldListType, fieldMedListType,&
                      paraListNb, paraListName,&
                      cmpListNb, cmpListName,&
                      nodeListNb, nodeListNume,&
                      cellListNb, cellListNume)
        character(len=8), intent(in) :: dsName
        character(len=8), intent(in) :: fileFormat
        integer(kind=8), intent(in) ::  fileUnit, fileVersion
        aster_logical, intent(in) :: lResu, lMesh, lField
        aster_logical, intent(in) :: lMeshCoor
        character(len=1), intent(in) :: paraFormat
        character(len=16), intent(in) :: realFormat
        character(len=4), intent(in) :: cplxFormat
        aster_logical, intent(in) :: lsup, linf, lmax, lmin
        real(kind=8), intent(in) :: borsup, borinf
        integer(kind=8), intent(in) :: storeListNb
        integer(kind=8) , pointer :: storeListIndx(:)
        integer(kind=8), intent(in) :: fieldListNb
        character(len=16), pointer :: fieldListType(:)
        character(len=80), pointer :: fieldMedListType(:)
        integer(kind=8), intent(in) :: paraListNb
        character(len=16), pointer :: paraListName(:)
        integer(kind=8), intent(in) :: cmpListNb
        character(len=8), pointer :: cmpListName(:)
        integer(kind=8), intent(in) :: nodeListNb
        integer(kind=8) , pointer :: nodeListNume(:)
        integer(kind=8), intent(in) :: cellListNb
        integer(kind=8) , pointer :: cellListNume(:)
    end subroutine irdebg
end interface
