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
#include "asterf_types.h"
!
interface
    subroutine dltlec(result, &
                      model, numeDOF, materField, mate, caraElem, &
                      jvMatr, masse, rigid, amort, lamort, &
                      nbLoad, nbVectAsse, listLoad, &
                      loadNameJv, loadInfoJv, loadFuncJv, jvVectAsse, jvVectFunc, &
                      nbWave, jvLoadWave, solveu, iinteg, t0, &
                      nume, numrep, ds_inout)
        use NonLin_Datastructure_type
        character(len=8), intent(out) :: result
        character(len=24), intent(out) :: model, numeDOF, materField, mate, caraElem
        integer(kind=8), intent(out) :: jvMatr(3)
        character(len=8), intent(out) :: masse, rigid, amort
        aster_logical, intent(out) :: lamort
        character(len=19), intent(out) :: listLoad, solveu
        character(len=24), intent(out) :: loadNameJv, loadInfoJv, loadFuncJv
        integer(kind=8), intent(out) :: nbVectAsse, nbLoad, nbWave
        integer(kind=8), intent(out) :: jvLoadWave, jvVectAsse, jvVectFunc
        integer(kind=8), intent(out) :: nume, numrep, iinteg
        real(kind=8), intent(out) :: t0
        type(NL_DS_InOut), intent(out) :: ds_inout
    end subroutine dltlec
end interface
