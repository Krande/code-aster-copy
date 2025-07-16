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
subroutine calipj(load, model)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'LIAISON_PROJ'
!
!
!     Création des relations entre les ddls des noeuds d'un maillage esclave et les ddls des
!     noeuds des mailles d'un maillage maître.
!
!     les relations sont fabriquées à partir de la sd_corresp_2_mailla sortant de proj_champ
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use LoadKinematic_module
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
!
    character(len=8), intent(in) :: load, model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: valeType = 'REEL'
    character(len=16), parameter :: factorKeyword = 'LIAISON_PROJ'
    character(len=19) :: listLineRela
    integer(kind=8) :: nbOcc
!
! --------------------------------------------------------------------------------------------------
!
    call getfac(factorKeyword, nbOcc)
    if (nbOcc .ne. 0) then
        call kineLoadLinkProj(model, valeType, listLineRela)
        call aflrch(listLineRela, load, 'LIN', detr_lisrez=ASTER_TRUE)
    end if
!
end subroutine calipj
