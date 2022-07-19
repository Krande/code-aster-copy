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

subroutine te0355(nomopt, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/laelem.h"
#include "asterfort/laMatr.h"
#include "asterfort/laParam.h"
#include "asterfort/laQuantities.h"
#include "asterfort/laVect.h"
#include "asterfort/utmess.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16) :: nomte, nomopt
!
! --------------------------------------------------------------------------------------------------
!
!   CHAR_MECA_CONT and RIGI_CONT for augmented Lagrangian method (Mortar)
!
! --------------------------------------------------------------------------------------------------
!
    integer :: elem_dime
    aster_logical :: l_axis
    integer :: nb_lagr_c, indi_lagc(9), nb_dofs
    character(len=8) :: elem_slav_code, elem_mast_code
    integer :: nb_node_slav, nb_node_mast
    real(kind=8) :: mast_coor_init(3, 9), slav_coor_init(3, 9)
    real(kind=8) :: mast_coor_curr(3, 9), slav_coor_curr(3, 9)
    real(kind=8) :: mast_depl_curr(3, 9), slav_depl_curr(3, 9)
    real(kind=8) :: proj_tole, gamma_c_nodes(4), lagc_curr(4)
    real(kind=8) :: vect(MAX_CONT_DOFS), matr(MAX_CONT_DOFS, MAX_CONT_DOFS)
!
! - Informations about finite element
!
    call laelem(nomte         , elem_dime     ,&
                l_axis        , &
                nb_dofs       , nb_lagr_c     , indi_lagc   ,&
                elem_slav_code, nb_node_slav,&
                elem_mast_code, nb_node_mast)
!
! - Get Parameters
!
    call laParam(proj_tole, gamma_c_nodes)
!
! - Get Quantities
!
    call laQuantities(elem_dime, nb_node_slav, nb_node_mast, &
                    indi_lagc, &
                    slav_coor_init, mast_coor_init, &
                    slav_coor_curr, mast_coor_curr, &
                    slav_depl_curr, mast_depl_curr, lagc_curr)
!
! - Computation
!
    if(nomopt == "CHAR_MECA_CONT") then
!
! --- Compute contact residual
!
        call laVect(elem_dime   , l_axis        , nb_dofs, &
                    nb_lagr_c   , indi_lagc     , lagc_curr, &
                    gamma_c_nodes, &
                    nb_node_slav, elem_slav_code, slav_coor_init, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    proj_tole, vect)
!
! --- Write vector
!
        call writeVector('PVECTCR', nb_dofs, vect)
!
    elseif(nomopt == "RIGI_CONT") then
!
! --- Compute contact matrix
!
        call laMatr(elem_dime   , l_axis        , nb_dofs, &
                    nb_lagr_c   , indi_lagc     , lagc_curr, &
                    gamma_c_nodes, &
                    nb_node_slav, elem_slav_code, slav_coor_init, slav_coor_curr, &
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    proj_tole, matr)
!
! - Write matrix
!
        call writeMatrix('PMATUUR', nb_dofs, nb_dofs, ASTER_TRUE, matr)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
