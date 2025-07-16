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

module HHO_compor_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/rcangm.h"
#include "asterfort/tecach.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics - Behaviour
!
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: HHO_Compor_State
    public :: isLargeStrain
    private :: initialize_compor, typmod1
!
    type HHO_Compor_State
!
        aster_logical      :: l_debug = ASTER_FALSE
! -----
        character(len=4)    :: fami = ' '
        character(len=8)    :: typmod(2) = [" ", " "]
        character(len=16)   :: option = " "
        character(len=16)   :: mult_comp = " "
!
        aster_logical       :: l_largestrain = ASTER_FALSE
        aster_logical       :: c_plan = ASTER_FALSE
        aster_logical       :: axis = ASTER_FALSE
!
        integer(kind=8)             :: nbsigm = 0
        integer(kind=8)             :: lgpg = 0
        integer(kind=8)             :: codret = 0
        integer(kind=8)             :: imater = 0
!
        real(kind=8)        :: angl_naut(3)
! --- pointer
        character(len=16), pointer :: compor(:) => null()
        real(kind=8), pointer :: carcri(:) => null()
        real(kind=8), pointer :: vari_prev(:) => null()
        real(kind=8), pointer :: sig_prev(:) => null()
        real(kind=8), pointer :: vari_curr(:) => null()
        real(kind=8), pointer :: sig_curr(:) => null()
! ----- member function
    contains
        procedure, pass :: initialize => initialize_compor
!
    end type HHO_Compor_State
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_compor(this, fami, option, ndim, bary)
!
        implicit none
!
        class(HHO_Compor_State), intent(inout)  :: this
        character(len=4), intent(in) :: fami
        character(len=16), intent(in) :: option
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: bary(3)
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a HHO_Compor_State type
!   In this     : a HHo Compor
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret, jmate, jtab(7)
        character(len=16), pointer :: v_mult(:) => null()
!
        this%fami = fami
        this%option = option
        this%typmod(1) = typmod1(ndim)
        this%typmod(2) = "HHO"
!
!
        this%axis = this%typmod(1) .eq. 'AXIS'
        this%c_plan = this%typmod(1) .eq. 'C_PLAN'
        this%nbsigm = nbsigm()
!
        call jevech('PMATERC', 'L', jmate)
        this%imater = zi(jmate-1+1)
!
        if (this%option .ne. "RIGI_MECA" .and. this%option .ne. "FORC_NODA" &
            .and. this%option .ne. "REFE_FORC_NODA") then
            call jevech('PCOMPOR', 'L', vk16=this%compor)
            call jevech('PCARCRI', 'L', vr=this%carcri)
            call jevech('PCONTMR', 'L', vr=this%sig_prev)
            call jevech('PVARIMR', 'L', vr=this%vari_prev)
            call jevech('PMULCOM', 'L', vk16=v_mult)
!
            call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
            ASSERT(iret .eq. 0)
            this%lgpg = max(jtab(6), 1)*jtab(7)
            this%mult_comp = v_mult(1)
            this%l_largestrain = isLargeStrain(this%compor(DEFO))
            call rcangm(ndim, bary, this%angl_naut)
        else
            this%l_largestrain = ASTER_FALSE
            if (this%option == "FORC_NODA") then
                call jevech('PSIEFR', 'L', vr=this%sig_prev)
                call jevech('PCOMPOR', 'L', vk16=this%compor)
                this%l_largestrain = isLargeStrain(this%compor(DEFO))
            end if
            call rcangm(ndim, bary, this%angl_naut)
        end if

        if (L_SIGM(option)) then
            call jevech('PCONTPR', 'E', vr=this%sig_curr)
        end if
!
        if (L_VARI(option)) then
            call jevech('PVARIPR', 'E', vr=this%vari_curr)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function isLargeStrain(defo_comp) result(bool)
!
        implicit none
!
        character(len=16), intent(in) :: defo_comp
        aster_logical :: bool
!
! --------------------------------------------------------------------------------------------------
!
!   To know if the model is in large deformations
!   In defo_comp : type of deformation
! --------------------------------------------------------------------------------------------------
!
        bool = ASTER_FALSE
!
        if (defo_comp(1:5) .eq. 'PETIT') then
            bool = ASTER_FALSE
        elseif (defo_comp .eq. 'GDEF_LOG') then
            bool = ASTER_TRUE
        elseif (defo_comp .eq. 'GREEN_LAGRANGE') then
            bool = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function typmod1(ndim) result(typmod)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        character(len=8) :: typmod
!
! --------------------------------------------------------------------------------------------------
!
!   Select typmod 1
! --------------------------------------------------------------------------------------------------
!
        select case (ndim)
        case (3)
            typmod = '3D'
        case (2)
            if (lteatt('AXIS', 'OUI')) then
                ASSERT(ASTER_FALSE)
                typmod = 'AXIS'
            else if (lteatt('C_PLAN', 'OUI')) then
                ASSERT(ASTER_FALSE)
                typmod = 'C_PLAN'
            else if (lteatt('D_PLAN', 'OUI')) then
                typmod = 'D_PLAN'
            else
                ASSERT(ASTER_FALSE)
            end if
        case default
            ASSERT(ASTER_FALSE)
        end select
    end function
!
end module
