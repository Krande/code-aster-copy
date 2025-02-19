/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org             */
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
/*-----------------------------------------------------------------------------/
/ Ouverture d'un fichier HDF, renvoie éventuellement une erreur
/  Paramètres :
/   - in  nomfic : nom du fichier (char *)
/  Résultats :
/     identificateur du fichier, -1 sinon (hid_t = int)
/-----------------------------------------------------------------------------*/
#ifdef ASTER_HAVE_HDF5
#include <hdf5.h>
#else
typedef int hid_t;
#endif

hid_t DEFS( HDFOPF, hdfopf, char *nomfic, STRING_SIZE ln ) {
    hid_t iret = -1;
#ifdef ASTER_HAVE_HDF5
    hid_t idfic;
    int k;
    char *nomf;
    void *malloc( size_t size );

    nomf = (char *)malloc( ( ln + 1 ) * sizeof( char ) );
    for ( k = 0; k < ln; k++ ) {
        nomf[k] = nomfic[k];
    }
    k = ln - 1;
    while ( nomf[k] == ' ' ) {
        k--;
    }
    nomf[k + 1] = '\0';
    // désactive les impressions HDF dans stdout
    H5Eset_auto1( NULL, NULL );
    if ( ( idfic = H5Fopen( nomf, H5F_ACC_RDONLY, H5P_DEFAULT ) ) >= 0 )
        iret = idfic;
    free( nomf );
#else
    CALL_UTMESS( "F", "FERMETUR_3" );
#endif
    return iret;
}
