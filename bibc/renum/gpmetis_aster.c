/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org             */
/* Copyright 1994-2011, Regents of the University of Minnesota          */
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
#include "metis.h"

#endif

void DEFPPPPPP( GPMETIS_ASTER, gpmetis_aster, ASTERINTEGER *nbnd, ASTERINTEGER *nadj,
                ASTERINTEGER4 *xadjd, ASTERINTEGER4 *adjncy, ASTERINTEGER *nbpart,
                ASTERINTEGER *partout ) {
#ifdef ASTER_HAVE_METIS

    idx_t i;
    idx_t options[METIS_NOPTIONS];
    idx_t *part;
    idx_t objval;
    int status = 0;
    idx_t ncon = 1;
    idx_t mnbpart = (idx_t)*nbpart;
    idx_t mnbnd = (idx_t)*nbnd;
    long k, l;

    METIS_SetDefaultOptions( options );
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
    options[METIS_OPTION_NO2HOP] = METIS_OPTION_NO2HOP;
    options[METIS_OPTION_MINCONN] = 0;
    options[METIS_OPTION_CONTIG] = 0;
    options[METIS_OPTION_SEED] = -1;
    options[METIS_OPTION_NITER] = 10;
    options[METIS_OPTION_NCUTS] = 1;
    options[METIS_OPTION_UFACTOR] = -1;
    options[METIS_OPTION_DBGLVL] = 0;
    idx_t *xadj, *adjnc;
    xadj = malloc( ( ( *nbnd ) + 1 ) * sizeof( idx_t ) );
    adjnc = malloc( ( *nadj ) * sizeof( idx_t ) );
    part = malloc( ( *nbnd ) * sizeof( idx_t ) );

    for ( k = 0; k <= *nbnd; k++ ) {
        xadj[k] = (idx_t)xadjd[k] - 1; /* -1 on est en C */
    }
    for ( k = 0; k < *nadj; k++ ) {
        adjnc[k] = (idx_t)adjncy[k] - 1; /* -1 on est en C */
    }

    status = METIS_PartGraphKway( &mnbnd, &ncon, xadj, adjnc, NULL, NULL, NULL, &mnbpart, NULL,
                                  NULL, options, &objval, part );
    if ( status != METIS_OK ) {
        printf( "\n***Metis returned with an error.\n" );
    } else {
        for ( i = 0; i < *nbnd; i++ ) {
            partout[i] = part[i];
        }
    }
    free( xadj );
    free( part );
    free( adjnc );
#endif
}
