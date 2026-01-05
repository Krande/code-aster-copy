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
subroutine verifh(famiZ, kpg, ksp, poumZ, jvMaterCode, epsyHydr)
!
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: famiZ
    integer(kind=8), intent(in) :: kpg, ksp
    character(len=*), intent(in) :: poumZ
    integer(kind=8), intent(in) :: jvMaterCode
    real(kind=8), intent(out) :: epsyHydr
!
! --------------------------------------------------------------------------------------------------
!
! Computation of autogenous shrinkage
! (inspired by verift.F90)
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  jvMaterCode      : adress for material parameters
! Out epsyHydr         : strain from autogenous shrinkage (retrait endogène)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: multiMater = " "
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCodePrev(nbProp), propCodeCurr(nbProp)
    character(len=8), parameter :: propName(nbProp) = "B_ENDOGE"
    real(kind=8) :: propVale(nbProp)
    character(len=8), parameter :: exteName = "HYDR"
    integer(kind=8) :: iretHydr, iretHydrPrev, iretHydrCurr
    real(kind=8) :: hydrPrev, hydrCurr
    real(kind=8) :: bendoPrev, bendoCurr
    character(len=1) :: poum
    character(len=16) :: valk(2)
    aster_logical :: lHasProp
    character(len=16) :: elasKeyword
    integer(kind=8) :: elasID
!
! --------------------------------------------------------------------------------------------------
!
    epsyHydr = 0.d0
    poum = poumZ

! - Detect external state variable
    iretHydr = 1
    call rcvarc(' ', exteName, '+', famiZ, kpg, &
                ksp, hydrCurr, iretHydr)

    if (iretHydr .eq. 0) then
! ----- Get values
        iretHydrPrev = 1
        hydrPrev = r8nnem()
        if (poum .eq. 'T' .or. poum .eq. '-') then
            call rcvarc(' ', exteName, '-', famiZ, kpg, &
                        ksp, hydrPrev, iretHydrPrev)
        end if
        iretHydrCurr = 1
        hydrCurr = r8nnem()
        if (poum .eq. 'T' .or. poum .eq. '+') then
            call rcvarc(' ', exteName, '+', famiZ, kpg, &
                        ksp, hydrCurr, iretHydrCurr)
        end if

! ----- Get parameters
        call get_elas_id(jvMaterCode, elasID, elasKeyword)
        propCodePrev = 1
        bendoPrev = r8nnem()
        if (poum .eq. 'T' .or. poum .eq. '-') then
            if (iretHydrPrev .eq. 0) then
                call rcvalb(famiZ, kpg, ksp, '-', &
                            jvMaterCode, multiMater, elasKeyword, &
                            0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCodePrev, 1)
                bendoPrev = propVale(1)
            end if
        end if
        propCodeCurr = 1
        bendoCurr = r8nnem()
        if (poum .eq. 'T' .or. poum .eq. '+') then
            if (iretHydrCurr .eq. 0) then
                call rcvalb(famiZ, kpg, ksp, '+', &
                            jvMaterCode, multiMater, elasKeyword, &
                            0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCodeCurr, 1)
                bendoCurr = propVale(1)
            end if
        end if

! ----- Test
        lHasProp = ASTER_FALSE
        if (poum .eq. 'T') then
            lHasProp = (propCodePrev(1)+propCodeCurr(1)) .eq. 0 .and. &
                       (iretHydrCurr+iretHydrPrev) .eq. 0
        elseif (poum .eq. '-') then
            lHasProp = propCodePrev(1) .eq. 0 .and. &
                       iretHydrPrev .eq. 0
        elseif (poum .eq. '+') then
            lHasProp = propCodeCurr(1) .eq. 0 .and. &
                       iretHydrCurr .eq. 0
        else
            ASSERT(ASTER_FALSE)
        end if
        if (.not. lHasProp) then
            valk(1) = exteName
            valk(2) = propName(1)
            call utmess('F', 'COMPOR5_32', nk=2, valk=valk)
        end if

! ----- Compute strains
        if (poum .eq. 'T') then
            if (iretHydrPrev+iretHydrCurr .eq. 0) then
                epsyHydr = (-bendoCurr*hydrCurr)-(-bendoPrev*hydrPrev)
            end if
        else if (poum .eq. '-') then
            if (iretHydrPrev .eq. 0) then
                epsyHydr = -bendoPrev*hydrPrev
            end if
        else if (poum .eq. '+') then
            if (iretHydrCurr .eq. 0) then
                epsyHydr = -bendoCurr*hydrCurr
            end if
        end if
    end if
!
end subroutine
