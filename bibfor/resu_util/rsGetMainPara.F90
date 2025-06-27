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
subroutine rsGetMainPara(phenom, resultZ, numeStore, &
                         listLoadZ, model, materField, mateco, caraElem, &
                         noLoads)
!
    use listLoad_type
    use listLoad_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/nmdoch.h"
#include "asterfort/ntdoch.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rslesd.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=4), intent(in) :: phenom
    character(len=*), intent(in) :: resultZ
    integer(kind=8), intent(in) :: numeStore
    character(len=*), intent(in) :: listLoadZ
    character(len=24), intent(out) :: mateco
    character(len=8), intent(out) :: model, materField, caraElem
    aster_logical, intent(out) :: noLoads
!
! --------------------------------------------------------------------------------------------------
!
! Results datastructure - Utility
!
! Get main parameters in results datastructure
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    character(len=8) :: result
    integer(kind=8) :: jvPara, nbLoadUser, nbLoadResu
    character(len=24) :: listLoadResu, listLoad
    type(ListLoad_Prep) :: listLoadPrep
    aster_logical :: lTher, lConsistent
!
! --------------------------------------------------------------------------------------------------
!
    result = resultZ
    caraElem = " "
    materField = " "
    mateco = " "
    model = " "
    listLoad = listLoadZ
    noLoads = ASTER_FALSE

! - Flag for thermal-dependent material parameters
    if (phenom .eq. "MECA") then
        lTher = ASTER_FALSE
    else if (phenom .eq. "THER") then
        lTher = ASTER_TRUE
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get model, cara_elem, etc. from result
    call rslesd(result, numeStore, model, materField, caraElem)
    if (materField .eq. " ") then
        mateco = ' '
    else
        call rcmfmc(materField, mateco, l_ther_=lTher)
    end if

! - Get loads from result datastructure
    call rsadpa(result, 'L', 1, 'EXCIT', numeStore, 0, sjv=jvPara)
    listLoadResu = zk24(jvPara)

! - Get loads/BC from user
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

! - Check consistency
    call getNbLoadsFromList(listLoad, nbLoadUser)
    nbLoadResu = 0

    if (listLoadResu .ne. " ") then
        call getNbLoadsFromList(listLoadResu, nbLoadResu)
        if (nbLoadUser .ne. 0) then
            call checkConsistency(listLoad, listLoadResu, lConsistent)
            if (.not. lConsistent) then
                call utmess('A', 'CHARGES9_53')
            end if
        end if
    end if
    if (nbLoadUser .eq. 0 .and. nbLoadResu .eq. 0) then
        noLoads = ASTER_TRUE
    end if
    if (nbLoadUser .eq. 0 .and. (listLoadResu .ne. " ")) then
        call copisd('LISTE_CHARGES', jvBase, listLoadResu, listLoad)
    end if
!
end subroutine
