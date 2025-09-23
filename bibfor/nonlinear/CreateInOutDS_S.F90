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

subroutine CreateInOutDS_S(ds_inout, l_temp_nonl)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(inout) :: ds_inout
    aster_logical, intent(in) :: l_temp_nonl
!
! --------------------------------------------------------------------------------------------------
!
! SECH_NON_LINE - Input/output management
!
! Create input/output datastructure
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_inout         : datastructure for input/output management
! In  l_temp_nonl      : flag for THER_NON_LINE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: nb_field_defi = 2
    integer(kind=8) :: i_field
! - Name of field (type) in results datastructure (add one -> don't forget to modify rscrsd.F90)
    character(len=16), parameter :: field_type(nb_field_defi) = &
                                    (/'SECH            ', 'COMPORTHER      '/)
! - Type of GRANDEUR for field
    character(len=8), parameter :: gran_name(nb_field_defi) = &
                                   (/'TEMP_R  ', 'COMPOR  '/)
! - Keyword for initial state (ETAT_INIT)
    character(len=8), parameter :: init_keyw(nb_field_defi) = &
                                   (/'CHAM_NO ', '        '/)
! - Spatial discretization of field
    character(len=4), parameter :: disc_type(nb_field_defi) = &
                                   (/'NOEU', 'ELGA'/)
! - TRUE if field can been read for initial state (ETAT_INIT)
    aster_logical, parameter :: l_read_init(nb_field_defi) = &
        (/.true._1, .false._1/)
! - TRUE if field can been store (ARCHIVAGE)
    aster_logical, parameter :: l_store(nb_field_defi) = &
        (/.true._1, .true._1/)
    ! - TRUE if field can been followed (OBSERVATION/SUIVI_DDL)
    aster_logical, parameter :: l_obsv(nb_field_defi) = &
        (/.true._1, .false._1/)
! - Keyword for OBSERVATION
    character(len=16), parameter :: obsv_keyw(nb_field_defi) = &
                                    (/'SECH            ', '                '/)
! - Variable (JEVEUX name) for field (#H# for hat variable)
    character(len=24), parameter :: algo_name(nb_field_defi) = &
                                    (/'#H#VALINC#SECH  ', 'XXXXXXXXXXXXXXXX'/)
! - Variable (JEVEUX name) for init field
    character(len=24), parameter :: init_name(nb_field_defi) = &
                                    (/'&&NTETCR.SECH0  ', 'XXXXXXXXXXXXXXXX'/)
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<DRYING> . Create input/output management datastructure'
    end if
!
! - Check
!
    ds_inout%l_temp_nonl = l_temp_nonl
    ds_inout%nb_field = nb_field_defi
    ASSERT(ds_inout%nb_field .le. ds_inout%nb_field_maxi)
!
! - Set list of fields
!
    do i_field = 1, nb_field_defi
        ds_inout%field(i_field)%type = field_type(i_field)
        ds_inout%field(i_field)%field_read = ' '
        ds_inout%field(i_field)%gran_name = gran_name(i_field)
        ds_inout%field(i_field)%obsv_keyw = obsv_keyw(i_field)
        ds_inout%field(i_field)%init_keyw = init_keyw(i_field)
        ds_inout%field(i_field)%disc_type = disc_type(i_field)
        ds_inout%field(i_field)%l_read_init = l_read_init(i_field)
        ds_inout%field(i_field)%l_store = l_store(i_field)
        ds_inout%field(i_field)%l_obsv = l_obsv(i_field)
        ds_inout%field(i_field)%algo_name = algo_name(i_field)
        ds_inout%field(i_field)%init_name = init_name(i_field)
        ds_inout%field(i_field)%init_type = ' '
        ds_inout%l_field_read(i_field) = .false._1
        ds_inout%l_field_acti(i_field) = .false._1
    end do
!
end subroutine
