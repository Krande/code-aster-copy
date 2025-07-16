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

module czm_elas_module

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics
    type MATERIAL
        real(kind=8) :: kn_t, kn_c, kt
        aster_logical:: an, at, cn
    end type MATERIAL

    ! CZM_ELAS class
    type CONSTITUTIVE_LAW
        integer(kind=8)       :: exception = 0
        integer(kind=8)       :: ndim
        real(kind=8)  :: r
        real(kind=8), dimension(:), allocatable:: phi
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(ndim, fami, kpg, ksp, imate, t, su) &
        result(self)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate, ndim
        real(kind=8), intent(in)     :: t(:), su(:)
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW)      :: self
! --------------------------------------------------------------------------------------------------
! ndim      displacement dimension
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! t         cohesive forces (local co-ordinates)
! su        displacement jump (local co-ordinates)
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter   :: nbel = 6, nblg = 1
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iok(nbel+nblg)
        real(kind=8)        :: valel(nbel), vallg(nblg)
        character(len=16)   :: nomel(nbel), nomlg(nblg)
! --------------------------------------------------------------------------------------------------
        data nomel/'RIGI_NOR', 'RIGI_TAN', 'ADHE_NOR', 'ADHE_TAN', 'RIGI_NOR_COMP', 'RIGI_NOR_TRAC'/
        data nomlg/'PENA_LAGR_ABSO'/
! --------------------------------------------------------------------------------------------------

        ! Dimension controls
        ASSERT(size(t) .eq. ndim)
        ASSERT(size(su) .eq. ndim)

        ! Parametres generaux
        self%ndim = ndim

        ! Augmentation coefficient
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_ELAS', 0, ' ', [0.d0], nblg, nomlg, &
                    vallg, iok, 2)
        self%r = vallg(1)
        ASSERT(self%r .gt. 0)

        ! Material parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_ELAS', 0, ' ', [0.d0], nbel, nomel, &
                    valel, iok, 0)

        self%mat%kn_c = merge(0.d0, valel(1), iok(1) .ne. 0)
        self%mat%kn_c = merge(self%mat%kn_c, valel(5), iok(5) .ne. 0)

        self%mat%kn_t = merge(0.d0, valel(1), iok(1) .ne. 0)
        self%mat%kn_t = merge(self%mat%kn_t, valel(6), iok(6) .ne. 0)

        self%mat%kt = merge(0.d0, valel(2), iok(2) .ne. 0)

        select case (nint(valel(3)))
        case (0)
            self%mat%cn = ASTER_FALSE
            self%mat%an = ASTER_FALSE
        case (1)
            self%mat%cn = ASTER_TRUE
            self%mat%an = ASTER_FALSE
        case (2)
            self%mat%cn = ASTER_FALSE
            self%mat%an = ASTER_TRUE
        case default
            ASSERT(ASTER_FALSE)
        end select

        select case (nint(valel(4)))
        case (0)
            self%mat%at = ASTER_FALSE
        case (1)
            self%mat%at = ASTER_TRUE
        case default
            ASSERT(ASTER_FALSE)
        end select

        ! Constitutive input phi = tau + r*su
        allocate (self%phi(ndim))
        self%phi = t+self%r*su

    end function Init

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! =====================================================================

    subroutine Integrate(self, delta, dphi_delta, vi)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(out)         :: delta(:), dphi_delta(:, :), vi(:)
! --------------------------------------------------------------------------------------------------
! dela          Gauss point displacement jump
! dphi_delta    derivative d(delta)/d(phi)
! vi            internal variable (post-treatment only)
! --------------------------------------------------------------------------------------------------
        integer(kind=8)         :: i
        real(kind=8)    :: rigi_norm
! --------------------------------------------------------------------------------------------------

        delta = 0
        dphi_delta = 0
        vi = 0

        ! Normal displacement
        if (.not. self%mat%an .and. .not. self%mat%cn) then
            rigi_norm = merge(self%mat%kn_t, self%mat%kn_c, self%phi(1) .gt. 0)
            delta(1) = self%phi(1)/(rigi_norm+self%r)
            dphi_delta(1, 1) = 1/(rigi_norm+self%r)

        else if (.not. self%mat%an .and. self%mat%cn) then
            delta(1) = max(0.d0, self%phi(1))/(self%mat%kn_t+self%r)
            dphi_delta(1, 1) = merge(1.d0, 0.d0, self%phi(1) .ge. 0)/(self%mat%kn_t+self%r)

        end if

        ! Tangential displacement
        if (.not. self%mat%at) then
            do i = 2, self%ndim
                delta(i) = self%phi(i)/(self%mat%kt+self%r)
                dphi_delta(i, i) = 1/(self%mat%kt+self%r)
            end do
        end if

        ! Internal variables storage
        vi(1:self%ndim) = delta

    end subroutine Integrate

end module czm_elas_module
