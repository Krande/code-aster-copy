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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine dbrParaInfoPod(operation, paraPod)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/romFieldPrintInfo.h"
#include "asterfort/romSnapInfo.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: operation
    type(ROM_DS_ParaDBR_POD), intent(in) :: paraPod
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Print informations about parameters - For POD methods
!
! --------------------------------------------------------------------------------------------------
!
! In  paraPod          : datastructure for parameters (POD)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: cmpFilterName
    character(len=16) :: variFilterName
    character(len=24) :: fieldName
    real(kind=8) :: toleSVD, toleIncr
    integer(kind=8) :: nbModeMaxi, nbCmpToFilter, iCmpToFilter, nbVariToFilter
    type(ROM_DS_Field) :: field
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_48')
    end if
!
! - Get parameters
!
    toleSVD = paraPod%toleSVD
    toleIncr = paraPod%toleIncr
    field = paraPod%field
    fieldName = paraPod%fieldName
    nbModeMaxi = paraPod%nbModeMaxi
    nbCmpToFilter = paraPod%nbCmpToFilter
    nbVariToFilter = paraPod%nbVariToFilter
!
! - Print - General for POD
!
    if (niv .ge. 2) then
        if (nbModeMaxi .ne. 0) then
            call utmess('I', 'ROM18_49', si=nbModeMaxi)
        end if
        call utmess('I', 'ROM18_51', sr=toleSVD)
        if (operation .eq. 'POD_INCR') then
            call utmess('I', 'ROM18_52', sr=toleIncr)
        end if
        call utmess('I', 'ROM18_50', sk=fieldName)
        call utmess('I', 'ROM18_64')
        call romFieldPrintInfo(field)
        if (nbCmpToFilter .eq. 0) then
            call utmess('I', 'ROM18_65')
        else
            call utmess('I', 'ROM18_66')
            do iCmpToFilter = 1, nbCmpToFilter
                cmpFilterName = paraPod%cmpToFilter(iCmpToFilter)
                if (nbVariToFilter .ne. 0) then
                    variFilterName = paraPod%variToFilter(iCmpToFilter)
                    call utmess('I', 'ROM18_68', si=iCmpToFilter, &
                                nk=2, valk=[cmpFilterName, variFilterName])
                else
                    call utmess('I', 'ROM18_67', si=iCmpToFilter, sk=cmpFilterName)
                end if
            end do
        end if
    end if
!
! - Print about snapshots selection
!
    if (niv .ge. 2) then
        call romSnapInfo(paraPod%snap)
    end if
!
end subroutine
