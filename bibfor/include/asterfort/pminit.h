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
#include "asterf_types.h"
!
interface
    subroutine pminit(jvMaterCode, nbVari, &
                      tablName, tablNbParaMaxi, tablNbPara, tablType, &
                      tablParaName, tablParaType, tablVale, &
                      anglNaut, pgl, lRota, &
                      epsiPrev, sigmPrev, vim, vip, &
                      loadEpsiType, loadType, loadFunc, coefImpo, &
                      coefMatrAdim, typeMatrPred, lMatrElas, matrElas, lPrintMatr, option, &
                      variName, nbVariTabl, &
                      sddisc, ds_conv, ds_algopara, sderro)
        use NonLin_Datastructure_type
        integer(kind=8), intent(in) :: jvMaterCode, nbVari
        character(len=8), intent(out) :: tablName
        integer(kind=8), intent(in) :: tablNbParaMaxi
        integer(kind=8), intent(out) :: tablNbPara, tablType
        character(len=16), intent(out) :: tablParaName(tablNbParaMaxi), tablParaType(tablNbParaMaxi)
        real(kind=8), intent(out) :: tablVale(tablNbParaMaxi)
        real(kind=8), intent(out) :: anglNaut(3), pgl(3, 3)
        aster_logical, intent(out) :: lRota
        real(kind=8), intent(out) :: epsiPrev(9), sigmPrev(6)
        real(kind=8), intent(out) :: vim(nbVari), vip(nbVari)
        integer(kind=8), intent(out) :: loadEpsiType, loadType(9)
        character(len=8), intent(out) :: loadFunc(9)
        real(kind=8), intent(out) :: coefImpo(6, 12), coefMatrAdim
        integer(kind=8), intent(out) :: typeMatrPred
        aster_logical, intent(out) :: lMatrElas
        real(kind=8), intent(out) :: matrElas(6, 6)
        aster_logical, intent(out) :: lPrintMatr
        character(len=16), intent(out) :: option
        character(len=8), intent(out) :: variName(nbVari)
        integer(kind=8), intent(out) :: nbVariTabl
        character(len=19), intent(out) :: sddisc
        type(NL_DS_Conv), intent(out) :: ds_conv
        type(NL_DS_AlgoPara), intent(out) :: ds_algopara
        character(len=24), intent(out) :: sderro
    end subroutine pminit
end interface
