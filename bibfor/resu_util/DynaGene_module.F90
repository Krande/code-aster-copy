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

module DynaGene_module

    implicit none

    private
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelibe.h"

    integer(kind=8), parameter :: nb_field = 7
  integer(kind=8), parameter :: ordr = 1, disc = 2, ptem = 3, depl = 4, vite = 5, acce = 6, vint = 7
    character(len=8), parameter :: field_names(nb_field) = (/'   .ORDR', '   .DISC', '   .PTEM', \
    '   .DEPL', '   .VITE', '   .ACCE', '.NL.VINT'/)

    type DynaGene
        integer(kind=8) :: ordr = ordr
        integer(kind=8) :: disc = disc
        integer(kind=8) :: ptem = ptem
        integer(kind=8) :: depl = depl
        integer(kind=8) :: vite = vite
        integer(kind=8) :: acce = acce
        integer(kind=8) :: vint = vint
        integer(kind=8) :: length = 0
        integer(kind=8) :: n_bloc = -1
        character(len=8), private :: sdname
        integer(kind=8), private, dimension(:), pointer :: v_shift => null()
        real(kind=8), private, dimension(:), pointer :: v_bloc => null()
        integer(kind=8), private, dimension(:), pointer :: v_blo2 => null()
        integer(kind=8), private, dimension(:), pointer :: v_ordr => null()
        real(kind=8), private, dimension(:), pointer :: v_disc => null()
        real(kind=8), private, dimension(:), pointer :: v_ptem => null()
        real(kind=8), private, dimension(:), pointer :: v_depl => null()
        real(kind=8), private, dimension(:), pointer :: v_vite => null()
        real(kind=8), private, dimension(:), pointer :: v_acce => null()
        real(kind=8), private, dimension(:), pointer :: v_vint => null()
        integer(kind=8), private :: current_bloc(nb_field) = -1
    contains
        procedure, public, pass :: init => init_dyna_gene
        procedure, public, pass :: free => free_dyna_gene
        procedure, public, pass :: get_values_by_index
        procedure, public, pass :: get_values_by_disc
        procedure, public, pass :: get_values_by_ordr
        procedure, public, pass :: get_values
        procedure, public, pass :: get_current_bloc
        procedure, public, pass :: has_field
        procedure, private, pass :: get_field_name
        procedure, private, pass :: free_field
        procedure, private, pass :: load_field
        procedure, private, pass :: update_field
        procedure, private, pass :: get_shift_length_pointer
    end type DynaGene

    public :: DynaGene

contains

! --------------------------------------------------------------------------------------------------

    subroutine get_field_name(this, field_type, i_bloc, field_name)
!
        implicit none
!
        class(DynaGene), intent(in) :: this
        integer(kind=8), intent(in) :: field_type, i_bloc
        character(len=24), intent(out) :: field_name

        character(len=7) :: intk7
        character(len=8) :: intk8

        if (i_bloc .ge. 0) then
            if (i_bloc .eq. 0) then
                intk8 = '        '
            else
                call codent(i_bloc, 'D0', intk7)
                intk8 = '.'//intk7
            end if
            field_name = this%sdname//intk8//field_names(field_type)
        else
            field_name = ' '
        end if

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine free_field(this, field_type)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type

        character(len=24) :: field_name

        call this%get_field_name(field_type, this%current_bloc(field_type), field_name)
        if (field_name .ne. ' ') then
            call jelibe(field_name)
            this%current_bloc(field_type) = -1
        end if

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine load_field(this, field_type, i_bloc)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type, i_bloc

        character(len=24) :: field_name

        this%current_bloc(field_type) = i_bloc
        call this%get_field_name(field_type, this%current_bloc(field_type), field_name)
        if (field_name .ne. ' ') then
            select case (field_type)
            case (ordr)
                call jeveuo(field_name, 'L', vi=this%v_ordr)
            case (disc)
                call jeveuo(field_name, 'L', vr=this%v_disc)
            case (ptem)
                call jeveuo(field_name, 'L', vr=this%v_ptem)
            case (depl)
                call jeveuo(field_name, 'L', vr=this%v_depl)
            case (vite)
                call jeveuo(field_name, 'L', vr=this%v_vite)
            case (acce)
                call jeveuo(field_name, 'L', vr=this%v_acce)
            case (vint)
                call jeveuo(field_name, 'L', vr=this%v_vint)
            case default
                ASSERT(.false.)
            end select
        end if

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine update_field(this, field_type, i_bloc)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type, i_bloc

        if (this%current_bloc(field_type) .ne. i_bloc) then
            call this%free_field(field_type)
            call this%load_field(field_type, i_bloc)
        end if

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine get_shift_length_pointer(this, field_type, shift, length, vr, vi)
!
        implicit none
!
        class(DynaGene), intent(in) :: this
        integer(kind=8), intent(in) :: field_type
        integer(kind=8), intent(out), optional :: shift, length
        integer(kind=8), intent(out), optional, dimension(:), pointer :: vi
        real(kind=8), intent(out), optional, dimension(:), pointer :: vr

        integer(kind=8):: i_bloc

        i_bloc = this%current_bloc(field_type)
        if (present(shift)) then
            if (i_bloc .eq. 0) then
                shift = 0
            else
                shift = this%v_shift(i_bloc)
            end if
        end if
        if (present(length)) then
            if (i_bloc .eq. 0) then
                length = this%length
            else
                length = this%v_shift(i_bloc+1)+1-this%v_shift(i_bloc)
            end if
        end if

        ASSERT(field_type .ne. ordr .and. present(vr) .or. field_type .eq. ordr .and. present(vi))

        select case (field_type)
        case (ordr)
            vi => this%v_ordr
        case (disc)
            vr => this%v_disc
        case (ptem)
            vr => this%v_ptem
        case (depl)
            vr => this%v_depl
        case (vite)
            vr => this%v_vite
        case (acce)
            vr => this%v_acce
        case (vint)
            vr => this%v_vint
        case default
            ASSERT(.false.)
        end select

    end subroutine

! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------

    subroutine init_dyna_gene(this, sdname)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        character(len=8), intent(in) :: sdname

        integer(kind=8) :: iret, length, i_bloc
        character(len=24) :: field_name

        this%sdname = sdname
        call jeexin(this%sdname//'           .BLOC', iret)
        if (iret .ne. 0) then
            call jeveuo(this%sdname//'           .BLOC', 'L', vr=this%v_bloc)
            call jelira(this%sdname//'           .BLOC', 'LONMAX', this%n_bloc)
            call jeveuo(this%sdname//'           .BLO2', 'L', vi=this%v_blo2)

            AS_ALLOCATE(this%n_bloc+1, vi=this%v_shift)
            this%v_shift(1) = 0

            do i_bloc = 1, this%n_bloc
                call this%get_field_name(ordr, i_bloc, field_name)
                call jelira(field_name, 'LONMAX', length)
                this%v_shift(i_bloc+1) = this%v_shift(i_bloc)+length-1
            end do

            this%length = this%v_shift(this%n_bloc+1)+1
        else
            call jeexin(this%sdname//'           .DISC', iret)
            if (iret .ne. 0) then
                this%n_bloc = 0
                call jelira(this%sdname//'           .DISC', 'LONMAX', this%length)
            end if
        end if

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine free_dyna_gene(this)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this

        integer(kind=8) :: field_type

        if (this%n_bloc .gt. 0) then
            call jelibe(this%sdname//'           .BLOC')
            call jelibe(this%sdname//'           .BLO2')
            AS_DEALLOCATE(vi=this%v_shift)
        end if
        this%n_bloc = -1

        do field_type = 1, nb_field
            call this%free_field(field_type)
        end do

    end subroutine

! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------

    subroutine get_values_by_index(this, field_type, idx, shift, length, vr, vi)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type, idx
        integer(kind=8), intent(out), optional :: shift, length
        real(kind=8), intent(out), optional, dimension(:), pointer :: vr
        integer(kind=8), intent(out), optional, dimension(:), pointer :: vi
!
        integer(kind=8) :: i_bloc

        ASSERT(idx .ge. 1 .and. idx .le. this%length)
        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        if (this%n_bloc .gt. 0) then
            do i_bloc = 1, this%n_bloc
                if (idx .le. this%v_shift(i_bloc+1)+1) exit
            end do
        else
            i_bloc = 0
        end if

        call this%update_field(field_type, i_bloc)
        call this%get_shift_length_pointer(field_type, shift, length, vr, vi)

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine get_values_by_disc(this, field_type, disc, shift, length, vr, vi)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type
        real(kind=8), intent(in) :: disc
        integer(kind=8), intent(out), optional :: shift, length
        real(kind=8), intent(out), optional, dimension(:), pointer :: vr
        integer(kind=8), intent(out), optional, dimension(:), pointer :: vi
!
        integer(kind=8) :: i_bloc

        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        if (this%n_bloc .gt. 0) then
            do i_bloc = 1, this%n_bloc-1
                if (disc .lt. this%v_bloc(i_bloc)) exit
            end do
        else
            i_bloc = 0
        end if

        call this%update_field(field_type, i_bloc)
        call this%get_shift_length_pointer(field_type, shift, length, vr, vi)
!
    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine get_values_by_ordr(this, field_type, ordr, shift, length, vr, vi)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type
        integer(kind=8), intent(in) :: ordr
        integer(kind=8), intent(out), optional :: shift, length
        real(kind=8), intent(out), optional, dimension(:), pointer :: vr
        integer(kind=8), intent(out), optional, dimension(:), pointer :: vi
!
        integer(kind=8) :: i_bloc

        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        if (this%n_bloc .gt. 0) then
            do i_bloc = 1, this%n_bloc
                if (ordr .le. this%v_blo2(i_bloc)) exit
            end do
            ASSERT(i_bloc .ne. this%n_bloc+1)
        else
            i_bloc = 0
        end if

        call this%update_field(field_type, i_bloc)
        call this%get_shift_length_pointer(field_type, shift, length, vr, vi)
!
    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine get_values(this, field_type, i_bloc, shift, length, vr, vi)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type, i_bloc
        integer(kind=8), intent(out), optional :: shift, length
        real(kind=8), intent(out), optional, dimension(:), pointer :: vr
        integer(kind=8), intent(out), optional, dimension(:), pointer :: vi
!

        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        if (this%n_bloc .gt. 0) then
            ASSERT(i_bloc .ge. 1 .and. i_bloc .le. this%n_bloc)
        else
            ASSERT(i_bloc .eq. 0)
        end if

        call this%update_field(field_type, i_bloc)
        call this%get_shift_length_pointer(field_type, shift, length, vr, vi)
!
    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine get_current_bloc(this, field_type, i_bloc)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type
        integer(kind=8), intent(out) :: i_bloc
!

        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        i_bloc = this%current_bloc(field_type)

    end subroutine

! --------------------------------------------------------------------------------------------------

    subroutine has_field(this, field_type, iret)
!
        implicit none
!
        class(DynaGene), intent(inout) :: this
        integer(kind=8), intent(in) :: field_type
        integer(kind=8), intent(out) :: iret

        character(len=24) :: field_name

        ASSERT(field_type .ge. 1 .and. field_type .le. nb_field)

        call this%get_field_name(field_type, this%n_bloc, field_name)
        call jeexin(field_name, iret)

    end subroutine

end module DynaGene_module
