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

subroutine ntnoli(model, mate, cara_elem, l_stat, l_evol, &
                  para, sddisc, ds_inout)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ntarch.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsrusd.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: mate
    character(len=24), intent(in) :: cara_elem
    aster_logical, intent(in) :: l_stat
    aster_logical, intent(in) :: l_evol
    real(kind=8), intent(in) :: para(*)
    character(len=19), intent(in) :: sddisc
    type(NL_DS_InOut), intent(inout) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! THER_LINEAIRE - Init
!
! Prepare storing
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  mate             : name of material characteristics (field)
! In  cara_elem        : name of elementary characteristics (field)
! In  l_stat           : .true. is stationnary
! In  l_evol           : .true. if transient
! In  para             : parameters for time
!                            (1) THETA
!                            (2) DELTAT
! In  sddisc           : datastructure for time discretization
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: sdarch
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    integer(kind=8) :: nume_store, nume_inst
    aster_logical :: force, lreuse
    character(len=8) :: result
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<THERMIQUE> PREPARATION DE LA SD EVOL_THER'
    end if
!
! - INSTANT INITIAL
!
    nume_inst = 0
    force = .true.
!
! - Get parameters in input/ouput management datastructure
!
    result = ds_inout%result
    lreuse = ds_inout%l_reuse
!
! --- ACCES SD ARCHIVAGE
!
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
!
! - Current storing index
!
    call jeveuo(sdarchAinfJv, 'L', vi=sdarchAinf)
    nume_store = sdarchAinf(1)
!
! --- CREATION DE LA SD EVOL_THER OU NETTOYAGE DES ANCIENS NUMEROS
!
    if (lreuse) then
        ASSERT(nume_store .ne. 0)
        call rsrusd(result, nume_store)
    else
        ASSERT(nume_store .eq. 0)
        call rscrsd('G', result, 'EVOL_THER', 100)
    end if
!
! - Stroing initial state
!
    if ((.not. lreuse) .and. (.not. l_stat) .and. l_evol) then
        call utmess('I', 'ARCHIVAGE_4')
        call ntarch(nume_inst, model, mate, cara_elem, para, &
                    sddisc, ds_inout, force)
    end if
!
end subroutine
