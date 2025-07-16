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
    subroutine iremed(fileUnit   , dsNameZ      , lResu           ,&
                      fieldListNb, fieldListType, fieldMedListType,&
                      storeListNb, storeListIndx,&
                      paraListNb , paraListName ,&
                      cmpListNb  , cmpListName  ,&
                      cellUserNb , cellUserNume ,&
                      nodeUserNb , nodeUserNume ,&
                      cplxFormat , lVariName    , caraElem,&
                      lfichUniq, lNomCas)
        integer(kind=8), intent(in) :: fileUnit
        character(len=19), intent(in) :: dsNameZ
        aster_logical, intent(in) :: lResu
        integer(kind=8), intent(in) :: fieldListNb
        character(len=16), pointer :: fieldListType(:)
        character(len=80), pointer :: fieldMedListType(:)
        integer(kind=8), intent(in) :: storeListNb
        integer(kind=8), pointer :: storeListIndx(:)
        integer(kind=8), intent(in) :: paraListNb
        character(len=16), pointer :: paraListName(:)
        integer(kind=8), intent(in) :: cmpListNb
        character(len=8), pointer :: cmpListName(:)
        integer(kind=8), intent(in) :: cellUserNb
        integer(kind=8), pointer :: cellUserNume(:)
        integer(kind=8), intent(in) :: nodeUserNb
        integer(kind=8), pointer :: nodeUserNume(:)
        character(len=*), intent(in) ::  cplxFormat
        aster_logical, intent(in) :: lVariName
        character(len=8), intent(in) :: caraElem
        aster_logical, intent(in) :: lfichUniq, lNomCas
    end subroutine iremed
end interface
