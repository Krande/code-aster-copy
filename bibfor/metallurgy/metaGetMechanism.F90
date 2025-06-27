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
subroutine metaGetMechanism(metaRela, metaGlob, &
                            l_plas, l_visc, &
                            l_hard_isotline, l_hard_isotnlin, &
                            l_hard_kine, l_hard_line, l_anneal, &
                            l_plas_tran)
!
    implicit none
!
#include "asterf_types.h"
!
    character(len=16), intent(in) :: metaRela, metaGlob
    aster_logical, optional, intent(out) :: l_plas
    aster_logical, optional, intent(out) :: l_visc
    aster_logical, optional, intent(out) :: l_hard_isotline, l_hard_isotnlin
    aster_logical, optional, intent(out) :: l_hard_kine
    aster_logical, optional, intent(out) :: l_hard_line
    aster_logical, optional, intent(out) :: l_anneal
    aster_logical, optional, intent(out) :: l_plas_tran
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Characteristics of comportment law
!
! --------------------------------------------------------------------------------------------------
!
! In  metaRela     : behaviour for each phase
! In  metaGlob     : global behaviour
! Out l_plas       : ASTER_TRUE if plasticity
! Out l_visc       : ASTER_TRUE if visco-plasticity
! Out l_hard_isot  : ASTER_TRUE if isotropic hardening
! Out l_hard_kine  : ASTER_TRUE if kinematic hardening
! Out l_hard_line  : ASTER_TRUE if linear hardening
! Out l_anneal     : ASTER_TRUE if restoration hardening
! Out l_plas_tran  : ASTER_TRUE if transformation plasticity
!
! --------------------------------------------------------------------------------------------------
!
    if (present(l_plas)) then
        l_plas = ASTER_FALSE
        if (metaRela(6:6) .eq. 'P') then
            l_plas = ASTER_TRUE
        end if
    end if
!
    if (present(l_visc)) then
        l_visc = ASTER_FALSE
        if (metaRela(6:6) .eq. 'V') then
            l_visc = ASTER_TRUE
        end if
    end if
!
    if (present(l_anneal)) then
        l_anneal = ASTER_FALSE
        if (metaGlob(12:16) .eq. '_RE  ' .or. metaGlob(12:16) .eq. '_PTRE') then
            l_anneal = ASTER_TRUE
        end if
    end if
!
    if (present(l_plas_tran)) then
        l_plas_tran = ASTER_FALSE
        if (metaGlob(12:16) .eq. '_PT  ' .or. metaGlob(12:16) .eq. '_PTRE') then
            l_plas_tran = ASTER_TRUE
        end if
    end if
!
    if (present(l_hard_isotline)) then
        l_hard_isotline = ASTER_FALSE
        if (metaRela(8:16) .eq. 'ISOT_LINE') then
            l_hard_isotline = ASTER_TRUE
        end if
    end if
!
    if (present(l_hard_isotnlin)) then
        l_hard_isotnlin = ASTER_FALSE
        if (metaRela(8:16) .eq. 'ISOT_TRAC') then
            l_hard_isotnlin = ASTER_TRUE
        end if
    end if
!
    if (present(l_hard_kine)) then
        l_hard_kine = ASTER_FALSE
        if (metaRela(8:16) .eq. 'CINE_LINE') then
            l_hard_kine = ASTER_TRUE
        end if
    end if
!
    if (present(l_hard_line)) then
        l_hard_line = ASTER_FALSE
        if (metaRela(8:16) .eq. 'ISOT_LINE' .or. metaRela(8:16) .eq. 'CINE_LINE') then
            l_hard_line = ASTER_TRUE
        end if
    end if
!
end subroutine
