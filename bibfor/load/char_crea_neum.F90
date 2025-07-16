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
subroutine char_crea_neum(load, model, mesh, geomDime, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/cachre.h"
#include "asterfort/char_eval_fonc.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: geomDime
    character(len=8), intent(in) :: model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Create Neumann loads
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh      : name of mesh
! In  load      : name of load
! In  geomDime  : space dimension
! In  model     : name of model
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: max_load_type
    parameter(max_load_type=7)
    integer(kind=8) :: nbocc(max_load_type)
    character(len=4) valeType2
    character(len=5) :: param(max_load_type)
    character(len=16) :: keywordfact(max_load_type)
!
    integer(kind=8) :: i
    character(len=5) :: curr_para
    data keywordfact/'FORCE_CONTOUR', 'FORCE_INTERNE', 'FORCE_ARETE',&
     &                 'FORCE_FACE', 'FORCE_POUTRE', 'FORCE_COQUE',&
     &                 'FORCE_COQUE_FO'/
    data param/'F1D2D', ' ', 'F1D3D',&
     &                 'F2D3D', 'F1D1D', ' ', ' '/
!
! --------------------------------------------------------------------------------------------------
!
!
!
! - Number of factor keywords
!
    do i = 1, max_load_type
        nbocc(i) = 0
        call getfac(keywordfact(i), nbocc(i))
    end do
!
! - Some checks: FORCE_FACE and FORCE_POUTRE prohibited en 2D
!
    if (geomDime .eq. 2) then
        if (nbocc(4) .ne. 0) then
            call utmess('F', 'CHARGES2_5', sk=keywordfact(4))
        end if
        if (nbocc(5) .ne. 0) then
            call utmess('F', 'CHARGES2_5', sk=keywordfact(5))
        end if
    end if
!
! - Load affectation
!
    valeType2 = valeType
!
    do i = 1, max_load_type
        if (nbocc(i) .ne. 0) then
            curr_para = param(i)
! --------- FORCE_INTERNE#2D
            if (keywordfact(i) .eq. 'FORCE_INTERNE' .and. geomDime .eq. 2) curr_para = 'F2D2D'
! --------- FORCE_INTERNE#3D
            if (keywordfact(i) .eq. 'FORCE_INTERNE' .and. geomDime .eq. 3) curr_para = 'F3D3D'
! --------- FORCE_COQUE#2D
            if (keywordfact(i) (1:11) .eq. 'FORCE_COQUE' .and. geomDime .eq. 2) curr_para = 'FCO2D'
! --------- FORCE_COQUE#3D
            if (keywordfact(i) (1:11) .eq. 'FORCE_COQUE' .and. geomDime .eq. 3) curr_para = 'FCO3D'
! --------- FORCE_COQUE_F
            if (keywordfact(i) .eq. 'FORCE_COQUE_FO') valeType2 = 'FONC'
!
            call cachre(load, model, mesh, geomDime, valeType2, &
                        curr_para, keywordfact(i))
!
            if (keywordfact(i) .eq. 'FORCE_COQUE_FO') then
                call char_eval_fonc(load, mesh, geomDime, curr_para)
            end if
        end if
    end do
!
end subroutine
