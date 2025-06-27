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
! ==================================================================================================
!
! Module for the management of integration of behaviour (MGIS)
!
! ==================================================================================================
!
module BehaviourMGIS_module
! ==================================================================================================
    use BehaviourMGIS_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: getMGISDime
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/BehaviourMGIS_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getMGISDime
!
! Get dimensions for MGIS
!
! In  lGreenLagr        : Green-Lagrange model
!  input: F tensor
!  output: stress PK2
!          matrix dPK2/dC
! In  lCZM              : CZM model
! In  lGradVari         : Grad_Vari model
! In  ndim              : space dimension (2 or 3)
! In  neps              : size of strain tensor for finite element (Voigt)
! In  nsig              : size of stress tensor for finite element (Voigt)
! In  nvi               : size of external state variables vector for finite element
! In  ndsde             : size of tangent matrix for finite element
! Out nstran            : size of strain tensor for MGIS (Voigt)
! Out nforc             : size of stress tensor for MGIS (Voigt)
! Out nstatv            : size of external state variables vector for MGIS
! Out nmatr             : size of tangent matrix for MGIS (always square)
!
! --------------------------------------------------------------------------------------------------
    subroutine getMGISDime(lGreenLagr, lCZM, lGradVari, ndim, &
                           neps, nsig, nvi, ndsde, &
                           nstran, nforc, nstatv, nmatr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: lGreenLagr, lCZM, lGradVari
        integer(kind=8), intent(in) :: ndim, neps, nsig, nvi, ndsde
        integer(kind=8), intent(out) :: nstran, nforc, nstatv, nmatr
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
        if (lGreenLagr) then
            ASSERT(neps .eq. 9)
            ASSERT(nsig .eq. 6)
            ASSERT(ndsde .eq. 36)
        elseif (lCZM) then
            ASSERT(neps .le. 6)
            ASSERT(nsig .le. 6)
            ASSERT(ndsde .le. 36)
        elseif (lGradVari) then
            ASSERT(neps .eq. 2*ndim+2+ndim)
            ASSERT(nsig .eq. 2*ndim+2+ndim)
            ASSERT(ndsde .eq. neps*nsig)
        else
            ASSERT(nsig .ge. 2*ndim)
            ASSERT(neps .ge. 2*ndim)
! --------- Because of DKT/DeBorst => not an equality
            ASSERT(neps*nsig .le. ndsde)
        end if

! ----- Deformation gradient
        if (lGreenLagr) then
            if (ndim .eq. 2) then
                nstran = 5
            elseif (ndim .eq. 3) then
                nstran = 9
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (lCZM) then
            nstran = ndim
        elseif (lGradVari) then
            nstran = 2*ndim+1
        else
            nstran = 2*ndim
        end if

! ----- Thermodynamic force
        if (lGreenLagr) then
            nforc = 2*ndim
        elseif (lCZM) then
            nforc = ndim
        elseif (lGradVari) then
            nstran = 2*ndim+1
        else
            nforc = 2*ndim
        end if

! ----- Internal state variables
        nstatv = nvi

! ----- Matrix
        nmatr = 2*ndim
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module BehaviourMGIS_module
