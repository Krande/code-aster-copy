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
subroutine ntarch(numins, model, materField, caraElem, para, &
                  sddisc, ds_inout, force, ds_algorom_, l_dry_)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/didern.h"
#include "asterfort/diinst.h"
#include "asterfort/dinuar.h"
#include "asterfort/nmarce.h"
#include "asterfort/ntarc0.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpg.h"
#include "asterfort/romAlgoNLTableSave.h"
!
    integer(kind=8), intent(in) :: numins
    character(len=8), intent(in) :: model, materField, caraElem
    real(kind=8), intent(in) :: para(*)
    character(len=19), intent(in) :: sddisc
    type(NL_DS_InOut), intent(in) :: ds_inout
    aster_logical, intent(inout) :: force
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
    aster_logical, optional, intent(in) :: l_dry_
!
! --------------------------------------------------------------------------------------------------
!
! THER_*  - Algorithm
! SECH_*  - Algorithm
!
! Storing results
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  ds_algorom       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: instan
    integer(kind=8) :: iret
    integer(kind=8) :: nume_store
    character(len=19) :: k19bid
    character(len=8) :: result
    character(len=24) :: listLoadResu
    aster_logical :: l_dry
!
! --------------------------------------------------------------------------------------------------
!
    result = ds_inout%result
    listLoadResu = ds_inout%listLoadResu

    l_dry = ASTER_FALSE
    if (present(l_dry_)) l_dry = l_dry_
!
! - Print timer
!
    call uttcpg('IMPR', 'INCR')
!
! - Last step => storing
!
    if (didern(sddisc, numins)) then
        force = .true.
    end if
!
! - Get index for storing
!
    call dinuar(result, sddisc, numins, force, &
                nume_store)
!
! - Current time
!
    instan = diinst(sddisc, numins)
!
! - Storing
!
    if (nume_store .ge. 0) then
!
! ----- Print head
!
        call utmess('I', 'ARCHIVAGE_5')
!
! ----- Increased result datastructure if necessary
!
        if (l_dry) then
            call rsexch(' ', result, 'SECH', nume_store, k19bid, &
                        iret)
        else
            call rsexch(' ', result, 'TEMP', nume_store, k19bid, &
                        iret)
        end if
        if (iret .eq. 110) then
            call rsagsd(result, 0)
        end if
!
! ----- Storing parameters
!
        call ntarc0(result, model, materField, caraElem, listLoadResu, &
                    para, nume_store, instan)
!
! ----- Storing fields
!
        call nmarce(ds_inout, result, sddisc, instan, nume_store, &
                    force)
!
! ----- Storing reduced parameters table (ROM)
!
        if (present(ds_algorom_)) then
            if (ds_algorom_%l_rom) then
                if (nume_store .gt. 0) then
                    call romAlgoNLTableSave(nume_store, instan, ds_algorom_)
                end if
            end if
        end if
    end if
!
end subroutine
