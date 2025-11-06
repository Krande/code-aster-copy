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
subroutine vdxtemp(kInf, kMoy, kSup, &
                   kpgsn, ksi3, &
                   hasTempRefe, tempRefe, &
                   hasTemp, tempKpg)
!
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: kInf, kMoy, kSup
    integer(kind=8), intent(in) :: kpgsn
    real(kind=8), intent(in) :: ksi3
    aster_logical, intent(in) :: hasTempRefe
    real(kind=8), intent(in) :: tempRefe
    aster_logical, intent(out) :: hasTemp
    real(kind=8), intent(out) :: tempKpg
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Get temperature and dilatation coefficient at current integration point
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret1, iret2, iret3
    real(kind=8) :: tempInf, tempMoy, tempSup
    real(kind=8) :: p1xi3, p2xi3, p3xi3
!
! --------------------------------------------------------------------------------------------------
!
    hasTemp = ASTER_FALSE
    tempKpg = r8nnem()

! - Get temperature
    call rcvarc(' ', 'TEMP', '+', 'MASS', kpgsn, &
                kInf, tempInf, iret1)
    call rcvarc(' ', 'TEMP', '+', 'MASS', kpgsn, &
                kMoy, tempMoy, iret2)
    call rcvarc(' ', 'TEMP', '+', 'MASS', kpgsn, &
                kSup, tempSup, iret3)
    p1xi3 = 1-ksi3*ksi3
    p2xi3 = -ksi3*(1-ksi3)/2.d0
    p3xi3 = ksi3*(1+ksi3)/2.d0
    if (hasTempRefe) then
        if ((iret1+iret2+iret3) .eq. 0) then
            tempKpg = tempMoy*p1xi3+tempInf*p2xi3+tempSup*p3xi3
            tempKpg = tempKpg-tempRefe
            hasTemp = ASTER_TRUE
        else
            call utmess('F', 'SHELL1_1')
        end if
    end if
!
end subroutine
