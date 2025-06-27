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
subroutine nmetl2(model, i_field, ds_inout)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/nmetcv.h"
#include "asterfort/nmetnc.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
!
    integer(kind=8), intent(in) :: i_field
    character(len=8), intent(in) :: model
    type(NL_DS_InOut), intent(inout) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Read field for ETAT_INIT - Field by field
!
! --------------------------------------------------------------------------------------------------
!
! In  i_field          : field index
! In  model            : model
! IO  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    aster_logical :: l_field_read
    character(len=24) :: field_read, field_read_cv, field_algo
    character(len=24) :: fieldType
    character(len=4) :: init_type, disc_type
    character(len=24) :: algo_name, fieldInDisc, init_name
!
! --------------------------------------------------------------------------------------------------
!
    field_read_cv = '&&NMETL2.CHAMP.CONV'

! - Field to read ?
    if (ds_inout%l_field_acti(i_field) .and. ds_inout%field(i_field)%l_read_init) then
! ----- Name of field (type) in results datastructure
        fieldType = ds_inout%field(i_field)%type

! ----- Name of field for initial state
        init_name = ds_inout%field(i_field)%init_name

! ----- Spatial discretization of field
        disc_type = ds_inout%field(i_field)%disc_type

! ----- Name of field in algorithm
        algo_name = ds_inout%field(i_field)%algo_name
        call nmetnc(algo_name, field_algo)

! ----- Actual state of field
        init_type = ds_inout%field(i_field)%init_type
!
! ----- Informations about field read in ETAT_INIT
!
        field_read = ds_inout%field(i_field)%field_read
        l_field_read = ds_inout%l_field_read(i_field)
!
! ----- Read initial field
!
        if (l_field_read) then
! --------- Get discretization of input field
            call dismoi('TYPE_CHAMP', field_read, 'CHAMP', repk=fieldInDisc, arret='C', ier=iret)

! --------- Try to convert field (discretization) if necessary and copy it
            call nmetcv(model, init_name, &
                        fieldType, field_read, fieldInDisc, field_read_cv, disc_type)
            if (disc_type .eq. 'NOEU') then
                call vtcopy(field_read_cv, field_algo, iret)
                if (iret .ne. 0) then
                    call utmess('A', 'MECANONLINE5_80', sk=fieldType)
                end if
            else if ((disc_type .eq. 'ELGA') .or. (disc_type .eq. 'ELEM') .or. &
                     (disc_type .eq. 'ELNO')) then
                call copisd('CHAMP_GD', 'V', field_read_cv, field_algo)
            else
                write (6, *) 'DISCRETISATION NON TRAITEE: ', fieldInDisc
                ASSERT(.false.)
            end if

! --------- New state of field
            ds_inout%field(i_field)%init_type = 'READ'
        end if
!
! ----- Copy initial field
!
        if (.not. l_field_read) then
            if (init_name .ne. ' ' .and. ds_inout%field(i_field)%init_type .eq. ' ') then
                call copisd('CHAMP', 'V', init_name, field_algo)
                ds_inout%field(i_field)%init_type = 'ZERO'
            end if
        end if
    end if
!
    call detrsd('CHAMP', field_read_cv)
!
end subroutine
