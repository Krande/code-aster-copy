/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "astercxx.h"

#include <iostream>
#include <vector>

VectorReal B_p1_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c );
VectorReal B_p2_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c );
VectorReal B_p4_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c );
VectorReal B_p5_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c );

extern "C" {
void BP1_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A );
void BP2_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A );
void BP4_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A );
void BP5_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A );
}
