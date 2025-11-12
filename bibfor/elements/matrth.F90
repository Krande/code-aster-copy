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
! aslint: disable=W0413
!
subroutine matrth(famiZ, &
                  elasID, elasKeywordZ, jvMaterCode, &
                  hasTemp, tempMoy, alpha, &
                  young_, nu_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: famiZ
    integer(kind=8), intent(in) :: elasID
    character(len=*), intent(in) :: elasKeywordZ
    integer(kind=8), intent(in) :: jvMaterCode
    aster_logical, intent(in) :: hasTemp
    real(kind=8), intent(in) :: tempMoy
    real(kind=8), intent(out) :: alpha
    real(kind=8), optional, intent(out) :: young_, nu_
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Get elasticity parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    integer(kind=8), parameter :: nbPropIsot = 3
    character(len=16), parameter :: propNameIsot(nbPropIsot) = &
                                    (/'E    ', 'NU   ', 'ALPHA'/)
    integer(kind=8), parameter :: nbPropOrth = 2
    character(len=16), parameter :: propNameOrth(nbPropOrth) = &
                                    (/'ALPHA_L', 'ALPHA_T'/)
    real(kind=8) :: propVale(nbPropIsot)
    integer(kind=8) :: propCode(nbPropIsot)
    real(kind=8) :: young, nu
    character(len=16) :: elasKeyword
    aster_logical :: hasAlphaTher
!
! --------------------------------------------------------------------------------------------------
!
    hasAlphaTher = ASTER_TRUE
    elasKeyword = elasKeywordZ
    young = 0.d0
    alpha = 0.d0
    nu = 0.d0
    if (elasID .eq. ELAS_ISOT) then
        call rcvala(jvMaterCode, ' ', elasKeyword, &
                    nbPara, paraName, [tempMoy], &
                    nbPropIsot, propNameIsot, &
                    propVale, propCode, 1)
        if (propCode(3) .ne. 0) then
            hasAlphaTher = ASTER_FALSE
        else
            young = propVale(1)
            nu = propVale(2)
            alpha = propVale(3)
        end if

    else if (elasID .eq. ELAS_ORTH) then
        call rcvalb(famiZ, 1, 1, '+', &
                    jvMaterCode, ' ', elasKeyword, &
                    0, paraName, [tempMoy], &
                    nbPropOrth, propNameOrth, &
                    propVale, propCode, 1)
        if (propCode(1) .ne. 0) then
            hasAlphaTher = ASTER_FALSE
        else
            if ((propVale(1) .eq. 0.d0) .and. (propVale(2) .eq. 0.d0)) then
                hasAlphaTher = ASTER_FALSE
            else
                call utmess('F', 'SHELL1_2')
            end if
        end if

    else
        call utmess('F', 'SHELL1_3', sk=elasKeyword)
    end if
    if (hasTemp .and. .not. hasAlphaTher) then
        call utmess('F', 'SHELL1_4')
    end if
    if (present(young_)) then
        young_ = young
    end if
    if (present(nu_)) then
        nu_ = nu
    end if
!
end subroutine
