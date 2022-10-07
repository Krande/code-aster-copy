! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0356(nomopt, nomte)
!
use contact_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/niElem.h"
#include "asterfort/laParam.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
#include "contact_module.h"
#include "jeveux.h"
!
    character(len=16) :: nomte, nomopt
!
! --------------------------------------------------------------------------------------------------
!
!   CHAR_MECA_CONT and RIGI_CONT for Nitsche's method
!
! --------------------------------------------------------------------------------------------------
!
    type(ContactParameters) :: parameters
    type(ContactGeom) :: geom
    real(kind=8) :: vect_cont(MAX_NITS_DOFS), vect_fric(MAX_NITS_DOFS)
    real(kind=8) :: matr_cont(MAX_NITS_DOFS, MAX_NITS_DOFS)
    real(kind=8) :: matr_fric(MAX_NITS_DOFS, MAX_NITS_DOFS)
!
! - Informations about finite element
!
    call niElem(nomte, geom, parameters)
!
! - Get Parameters
!
    call laParam(parameters)
!
! - Verif
!
    ASSERT(parameters%algo_cont == CONT_ALGO_NITS)
    ASSERT(parameters%type_fric == FRIC_TYPE_NONE)
!
! - Computation
!
    if(nomopt == "CHAR_MECA_CONT") then
!
! --- Compute contact residual
!
!        call laVect(parameters, geom, vect_cont, vect_fric)
!
! --- Write vector
!
        call writeVector('PVECTCR', geom%nb_dofs, vect_cont)
!
        if(parameters%l_fric) then
            call writeVector('PVECTFR', geom%nb_dofs, vect_fric)
        end if
!
    elseif(nomopt == "RIGI_CONT") then
!
! --- Compute contact matrix
!
!        call laMatr(parameters, geom, matr_cont, matr_fric)
!
! - Write matrix
!
        call writeMatrix('PMATUUR', geom%nb_dofs, geom%nb_dofs, ASTER_TRUE, matr_cont)
!
        if(parameters%l_fric) then
            call writeMatrix('PMATUNS', geom%nb_dofs, geom%nb_dofs, ASTER_FALSE, matr_fric)
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
