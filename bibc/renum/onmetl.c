/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org             */
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

#include "aster.h"

#include "aster_fort_utils.h"

#ifdef ASTER_HAVE_METIS
#include "GKlib.h"
#include "metis.h"

/* Prototypes of internal functions */
// int ReadGraphL( gk_graph_t *, int *, int *, int *, int *, int *, int * );
// int ComputeFillInL( gk_graph_t *, idx_t *, int *, int *, int *, int *, double *, int *, int * );
// int smbfctl( int, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *, int *, idx_t *, idx_t *, int *,
//              int *, int *, int *, int *, double * );

#endif

void DEFPPPPPPPPPPPPP( ONMETL, onmetl, nbnd, nadj, xadjd, adjncy, invpnd, permnd, supnd, parent,
                       nbsn, nbops, fctnzs, lgind, niv )

    ASTERINTEGER4 *nbnd,
    *nadj, *xadjd, *adjncy;
ASTERINTEGER4 *invpnd, *permnd, *supnd, *parent, *nbsn;
ASTERDOUBLE *nbops;
ASTERINTEGER4 *fctnzs, *lgind;
ASTERINTEGER *niv;
{
#ifdef ASTER_HAVE_METIS

#endif
}
