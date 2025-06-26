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
subroutine getMainPara(phenom, &
                       model, materField, mateco, caraElem, listLoad)
!
    use listLoad_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/nmdoch.h"
#include "asterfort/ntdoch.h"
!
    character(len=4), intent(in) :: phenom
    character(len=24), intent(out) :: mateco
    character(len=8), intent(out) :: model, materField, caraElem
    character(len=24), intent(out) :: listLoad
!
! --------------------------------------------------------------------------------------------------
!
! Get main parameters from user
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    integer(kind=8) :: iret
    type(ListLoad_Prep) :: listLoadPrep
    aster_logical :: lTher
! --------------------------------------------------------------------------------------------------
!
    caraElem = " "
    materField = " "
    mateco = " "
    model = " "

! - Flag for thermal-dependent material parameters
    if (phenom .eq. "MECA") then
        lTher = ASTER_FALSE
    else if (phenom .eq. "THER") then
        lTher = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if

    call getvid(' ', 'MODELE', scal=model)
    call getvid(' ', 'CHAM_MATER', scal=materField)
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=iret)
    if (iret .eq. 0) then
        caraElem = ' '
    end if
    if (materField .eq. " ") then
        mateco = ' '
    else
        call rcmfmc(materField, mateco, l_ther_=lTher)
    end if

! - Get loads/BC
    listLoad = '&&MAINPARA.LISCHA'
    listLoadPrep%model = model(1:8)
    listLoadPrep%lHasPilo = ASTER_FALSE
! - Pas normal ! Que fait-on en harmonique ?
    listLoadPrep%funcIsCplx = ASTER_FALSE
    if (phenom .eq. "MECA") then
        call nmdoch(listLoadPrep, listLoad, jvBase)
    else if (phenom .eq. "THER") then
        call ntdoch(listLoadPrep, listLoad, jvBase)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
