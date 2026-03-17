! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

module czm_frot_mix_module

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, Integrate

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics
    type MATERIAL
        real(kind=8) :: kn, kt, cohe, frot
        aster_logical :: ad
    end type MATERIAL

    ! CZM_FROT_MIX class
    type CONSTITUTIVE_LAW
        integer(kind=8) :: exception = 0
        aster_logical :: elas, rigi, pred
        integer(kind=8) :: ndim
        real(kind=8) :: r
        real(kind=8), dimension(:), allocatable:: phi, deltap
        integer(kind=8) :: statep
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(ndim, option, fami, kpg, ksp, imate, t, su, vim) result(self)

        implicit none

        integer(kind=8), intent(in) :: kpg, ksp, imate, ndim
        real(kind=8), intent(in) :: t(:), su(:), vim(:)
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW) :: self
! --------------------------------------------------------------------------------------------------
! ndim      displacement jump dimension
! option    computation option
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! t         cohesive forces (local co-ordinates)
! su        jump jump (local co-ordinates)
! vim       internal variables at the beginning of the time step
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter :: nbel = 2, nbpl = 3, nblg = 1
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: iokel(nbel), iokpl(nbpl), ioklg(nblg)
        real(kind=8) :: valel(nbel), valpl(nbpl), vallg(nblg)
        character(len=16) :: nomel(nbel), nompl(nbpl), nomlg(nblg)
! --------------------------------------------------------------------------------------------------
        data nomel/'RIGI_NOR', 'RIGI_TAN'/
        data nompl/'ADHE', 'COHESION', 'COEF_FROT'/
        data nomlg/'PENA_LAGR_ABSO'/
! --------------------------------------------------------------------------------------------------

        ! Dimension controls
        ASSERT(size(t) .eq. ndim)
        ASSERT(size(su) .eq. ndim)

        ! General parameters
        self%ndim = ndim

        ! Options
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'

        ! Material elastic parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nbel, nomel, &
                    valel, iokel, 0)

        ! Material plastic parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nbpl, nompl, &
                    valpl, iokpl, 2)

        ! Augmentation coefficient
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'CZM_FROT_MIX', 0, ' ', [0.d0], nblg, nomlg, &
                    vallg, ioklg, 2)

        self%mat%kn = valel(1)
        self%mat%kt = valel(2)

        select case (nint(valpl(1)))
        case (0)
            self%mat%ad = ASTER_FALSE
        case (1)
            self%mat%ad = ASTER_TRUE
        end select
        self%mat%cohe = valpl(2)
        self%mat%frot = valpl(3)

        self%r = vallg(1)

        ! Check that RIGI_NOR and RIGI_TAN are found if ADHE='ELAS'
        if ((.not. self%mat%ad) .and. &
            (iokel(1) .eq. 1 .or. iokel(2) .eq. 1)) call utmess('F', 'MECANONLINE3_4')

        ! Constitutive input phi = tau + r*su
        allocate (self%phi(ndim))
        self%phi = t+self%r*su

        ! Previous plastic jump
        allocate (self%deltap(ndim-1))
        self%deltap = vim(1:ndim-1)

        ! Previous plastic state
        self%statep = vim(6)

    end function Init

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! =====================================================================

    subroutine Integrate(self, delta, dphi_delta, vi)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout) :: self
        real(kind=8), intent(out) :: delta(:), dphi_delta(:, :), vi(:)
! --------------------------------------------------------------------------------------------------
! delta         Gauss point jump jump
! dphi_delta    derivative d(delta)/d(phi)
! vi            internal variables
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: i, j
        real(kind=8) :: tel(self%ndim), nel(self%ndim), id(self%ndim, self%ndim)
        real(kind=8) :: fel, dka, telq, teln
        integer(kind=8) :: state
!--------------------------------------------------------------------------------------------------

        ! Initialisation
        delta = 0.d0
        dphi_delta = 0.d0
        dka = 0.d0
        tel = 0.d0
        id = 0.d0
        forall (i=1:self%ndim) id(i, i) = 1.d0

! ======================================================================
!  ELASTIC PREDICTION
! ======================================================================

        if (.not. self%mat%ad) then
            tel(1) = self%phi(1)/(1.d0+self%r/self%mat%kn)
            if (self%mat%frot .gt. 0.d0) tel(1) = min(tel(1), self%mat%cohe/self%mat%frot)
            do i = 2, self%ndim
                tel(i) = (self%phi(i)-self%r*self%deltap(i-1))/(1.d0+self%r/self%mat%kt)
            end do
        else
            tel(1) = self%phi(1)
            if (self%mat%frot .gt. 0.d0) tel(1) = min(tel(1), self%mat%cohe/self%mat%frot)
            do i = 2, self%ndim
                tel(i) = self%phi(i)-self%r*self%deltap(i-1)
            end do
        end if
        teln = tel(1)
        telq = sqrt(dot_product(tel(2:self%ndim), tel(2:self%ndim)))
        fel = telq+self%mat%frot*teln-self%mat%cohe

! ======================================================================
!  COMPUTATION OF THE PLASTIC JUMP
! ======================================================================

        ! Elastic regime
        if (fel .lt. 0.d0) then

            ! Elastic state
            state = 0

            ! Compute the plastic jump (the true internal variables)
            vi(1:self%ndim-1) = self%deltap(1:self%ndim-1)

            ! Plastic regime
        else

            ! Plastic state
            state = 1

            ! Compute the plastic jump (the true internal variables)
            if (.not. self%mat%ad) then
                dka = fel/(self%r/(1.d0+self%r/self%mat%kt))
            else
                dka = fel/self%r
            end if
            do i = 2, self%ndim
                nel(i) = tel(i)/telq
                vi(i-1) = self%deltap(i-1)+dka*nel(i)
            end do

        end if

! ======================================================================
!  COMPUTATION OF DELTA
! ======================================================================

        if (.not. self%mat%ad) then
            if (self%mat%frot*teln .lt. self%mat%cohe) then
                delta(1) = self%phi(1)/(self%mat%kn+self%r)
            else
                delta(1) = (self%phi(1)-self%mat%cohe/self%mat%frot)/self%r
            end if
            do i = 2, self%ndim
                delta(i) = (self%phi(i)+self%mat%kt*vi(i-1))/(self%mat%kt+self%r)
            end do
        else
            if (self%mat%frot*teln .lt. self%mat%cohe) then
                delta(1) = 0.d0
            else
                delta(1) = (self%phi(1)-self%mat%cohe/self%mat%frot)/self%r
            end if
            do i = 2, self%ndim
                delta(i) = vi(i-1)
            end do
        end if

! ======================================================================
!  COMPUTATION OF THE POSTPROCESSING INTERNAL VARIABLES
! ======================================================================

        vi(3:self%ndim+2) = delta(1:self%ndim)
        vi(6) = state
        if (self%mat%frot*teln .lt. self%mat%cohe) then
            vi(7) = 0
        else
            vi(7) = 1
        end if

! ======================================================================
!  COMPUTATION OF THE TANGENT MATRIX
! ======================================================================

        if (.not. self%rigi) goto 999

        ! Tangent matrix selection for the prediction
        if (self%pred) then
            state = self%statep
        end if

        if (state .eq. 0 .or. self%elas) then

            if (.not. self%mat%ad) then
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    dphi_delta(1, 1) = 1.d0/(self%mat%kn+self%r)
                else
                    dphi_delta(1, 1) = 1.d0/self%r
                end if
                do i = 2, self%ndim
                    dphi_delta(i, i) = 1.d0/(self%mat%kt+self%r)
                end do
            else
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    dphi_delta(1, 1) = 0.d0
                else
                    dphi_delta(1, 1) = 1.d0/self%r
                end if
                do i = 2, self%ndim
                    dphi_delta(i, i) = 0.d0
                end do
            end if

        else

            ! d(delta(1))/d(phi(1))
            if (.not. self%mat%ad) then
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    dphi_delta(1, 1) = 1.d0/(self%mat%kn+self%r)
                else
                    dphi_delta(1, 1) = 1.d0/self%r
                end if
            else
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    dphi_delta(1, 1) = 0.d0
                else
                    dphi_delta(1, 1) = 1.d0/self%r
                end if
            end if

            ! d(delta(i))/d(phi(1)), for i=2 to ndim
            if (.not. self%mat%ad) then
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    do i = 2, self%ndim
                        dphi_delta(i, 1) = self%mat%frot*nel(i)/self%r &
                                           /(1.d0+self%r/self%mat%kn)
                    end do
                end if
            else
                if (self%mat%frot*teln .lt. self%mat%cohe) then
                    do i = 2, self%ndim
                        dphi_delta(i, 1) = self%mat%frot*nel(i)/self%r
                    end do
                end if
            end if

            ! d(delta(i))/d(phi(j)), for i,j=2 to ndim
            if (.not. self%mat%ad) then
                do i = 2, self%ndim
                    do j = 2, self%ndim
                        dphi_delta(i, j) = id(i, j)/(self%mat%kt+self%r) &
                                           +(nel(i)*nel(j)+fel/telq*(id(i, j)-nel(i)*nel(j))) &
                                           /self%r/(1.d0+self%r/self%mat%kt)
                    end do
                end do
            else
                do i = 2, self%ndim
                    do j = 2, self%ndim
                        dphi_delta(i, j) = (nel(i)*nel(j)+fel/telq*(id(i, j)-nel(i)*nel(j)))/self%r
                    end do
                end do
            end if

        end if

999     continue

    end subroutine Integrate

end module czm_frot_mix_module
