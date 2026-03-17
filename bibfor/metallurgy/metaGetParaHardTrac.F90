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
!
subroutine metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                               l_temp, temp, &
                               epseq, h0, rp_, nbValeMaxi_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rctype.h"
#include "asterfort/rctrac.h"
#include "asterfort/rcfonc.h"
#include "asterfort/utmess.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: jvMaterCode
    integer(kind=8), intent(in) :: metaType
    integer(kind=8), intent(in) :: nbPhase
    aster_logical, intent(in) :: l_temp
    real(kind=8), intent(in) :: temp
    real(kind=8), intent(in) :: epseq(nbPhase)
    real(kind=8), intent(out) :: h0(nbPhase)
    real(kind=8), optional, intent(out) :: rp_(nbPhase)
    integer(kind=8), optional, intent(out) :: nbValeMaxi_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get hardening slope (non-linear)
!
! --------------------------------------------------------------------------------------------------
!
! In  jvMaterCode  : coded material address
! In  metaType     : type of metallurgy
! In  nbPhase      : total number of phases (cold and hot)
! In  l_temp       : .true. if temperature command variable is affected
! In  temp         : temperature
! In  epseq        : cumulated plastic strain
! Out h0           : current hardening slope
! Out rp           : current isotropic hardening value
! Out nbValeMaxi   : maximum number of points of traction curves for all phases
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: tractionName = 'META_TRACTION'
    integer(kind=8), parameter :: nbProp = 5
    character(len=8) :: propName(nbProp)
    integer(kind=8) :: iPhase
    integer(kind=8) :: jvProl, jvVale
    integer(kind=8) :: nbVale, nbValeMaxi
    character(len=8) :: paraName
    real(kind=8) :: paraVale, young
!
! --------------------------------------------------------------------------------------------------
!

! - Select name of properties for traction curve
    if (metaType .eq. META_STEEL) then
        propName(1) = 'SIGM_F1'
        propName(2) = 'SIGM_F2'
        propName(3) = 'SIGM_F3'
        propName(4) = 'SIGM_F4'
        propName(5) = 'SIGM_C'

    elseif (metaType .eq. META_ZIRC) then
        propName(1) = 'SIGM_F1'
        propName(2) = 'SIGM_F2'
        propName(3) = 'SIGM_C'

    else
        ASSERT(ASTER_FALSE)
    end if

!
    nbValeMaxi = -1
    do iPhase = 1, nbPhase
! ----- Get informations (value and type) about parameters for traction curve
        call rctype(jvMaterCode, 1, 'TEMP', [temp], &
                    paraVale, paraName, &
                    keyw_factz=tractionName, keywz=propName(iPhase))
        if ((paraName .eq. 'TEMP') .and. (.not. l_temp)) then
            call utmess('F', 'COMPOR5_5', sk=paraName)
        end if

! ----- Get access of traction curve
        call rctrac(jvMaterCode, 2, propName(iPhase), temp, &
                    jvProl, jvVale, nbVale, young)

        if (present(rp_)) then
            call rcfonc('V', 2, jvProl, jvVale, nbVale, &
                        p=epseq(iPhase), rp=rp_(iPhase), rprim=h0(iPhase))
        else
            call rcfonc('V', 2, jvProl, jvVale, nbVale, &
                        p=epseq(iPhase), rprim=h0(iPhase))
        end if
        if (nbVale .ge. nbValeMaxi) then
            nbValeMaxi = nbVale
        end if
    end do
    if (present(nbValeMaxi_)) then
        nbValeMaxi_ = nbValeMaxi
    end if
!
end subroutine
