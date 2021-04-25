! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine cmePrep(optionz, modelz, timeCurr, timeIncr, chtime)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mecact.h"
!
character(len=*), intent(in) :: optionz, modelz
real(kind=8), intent(in) :: timeCurr, timeIncr
character(len=24), intent(out) :: chtime
!
! --------------------------------------------------------------------------------------------------
!
! CALC_MATR_ELEM
!
! Preparation
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  model            : name of the model
! In  timeCurr         : current time
! In  timeIncr         : time step
! Out chtime           : time parameters (field)
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbCmp = 6
    character(len=8), parameter :: cmpName(nbCmp) = (/'INST    ','DELTAT  ','THETA   ',&
                                                      'KHI     ','R       ','RHO     '/)
    real(kind=8) :: cmpVale(nbCmp)
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
    option       = optionz
    chtime       = '&&CHTIME'
    cmpVale(1:6) = [timeCurr, timeIncr, 1.d0, 0.d0, 0.d0, 0.d0]
!
    if ((option.eq.'RIGI_THER') .or. (option.eq.'MASS_THER')) then
        call mecact('V', chtime, 'MODELE', modelz, 'INST_R',&
                    ncmp=nbCmp, lnomcmp=cmpName, vr=cmpVale)
    endif
!
end subroutine
