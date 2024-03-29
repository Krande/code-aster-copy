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
/ Ecriture sur un fichier HDF d'un segment de valeur associé à un objet JEVEUX
/  Paramètres :
/   - in  idfic  : identificateur du fichier (hid_t)
/   - in  nomg   : nom du groupe (char *)
/   - in  nomdts : nom du dataset (char *)
/   - in  type   : type des valeurs stockées (char *)
/   - in  sv     : valeurs associées
/  Résultats :
/     identificateur du fichier, -1 sinon (hid_t = int)
/-----------------------------------------------------------------------------*/
#ifdef ASTER_HAVE_HDF5
#include <hdf5.h>
#else
typedef int hid_t;
#endif
#include <stdlib.h>

ASTERINTEGER DEFPSSSPSP( HDFWSV, hdfwsv, hid_t *idf, char *nomg, STRING_SIZE lg, char *nomdts,
                         STRING_SIZE ln, char *type, STRING_SIZE lt, ASTERINTEGER *ltype, char *sv,
                         STRING_SIZE toto, ASTERINTEGER *lsv ) {
#ifdef ASTER_HAVE_HDF5
    hid_t idfic, datatype, dataspace, dataset, type_id;
    herr_t iret;
    hsize_t dimsf[1];
    int lg2, lmot;
    char *nomd, *vtype, *mot = NULL, *pmot;
    int k;
    void *malloc( size_t size );

    idfic = (hid_t)*idf;
    nomd = (char *)malloc( ( lg + ln + 2 ) * sizeof( char ) );
    for ( k = 0; k < lg; k++ ) {
        nomd[k] = nomg[k];
    }
    k = lg - 1;
    while ( k >= 0 ) {
        if ( nomd[k] == ' ' ) {
            k--;
        } else
            break;
    }
    nomd[k + 1] = '/';
    lg2 = k + 1 + 1;
    for ( k = 0; k < ln; k++ ) {
        nomd[lg2 + k] = nomdts[k];
    }
    k = lg2 + ln - 1;
    while ( k >= 0 ) {
        if ( nomd[k] == ' ' ) {
            k--;
        } else
            break;
    }
    nomd[k + 1] = '\0';

    vtype = (char *)malloc( ( lt + 1 ) * sizeof( char ) );
    for ( k = 0; k < lt; k++ ) {
        vtype[k] = type[k];
    }
    vtype[lt] = '\0';
    /*
     *   Type à déterminer en fonction de l'argument type
     */
    dimsf[0] = (hsize_t)*lsv;
    if ( strcmp( vtype, "R" ) == 0 ) {
        type_id = H5T_NATIVE_DOUBLE;
    } else if ( strcmp( vtype, "C" ) == 0 ) {
        type_id = H5T_NATIVE_DOUBLE;
        dimsf[0] = (hsize_t)*lsv;
    } else if ( strcmp( vtype, "I" ) == 0 ) {
#ifdef ASTER_HAVE_LONG_LONG
        type_id = H5T_NATIVE_LLONG;
#else
        type_id = H5T_NATIVE_LONG;
#endif
    } else if ( strcmp( vtype, "S" ) == 0 ) {
        type_id = H5T_NATIVE_INT;
    } else if ( strcmp( vtype, "L" ) == 0 ) {
        type_id = H5T_NATIVE_HBOOL;
    } else if ( strcmp( vtype, "K" ) == 0 ) {
        type_id = H5T_FORTRAN_S1;
        pmot = (char *)sv;
        lmot = (int)( *lsv * ( *ltype ) );
        mot = (char *)malloc( lmot * sizeof( char ) );
        for ( k = 0; k < *lsv; k++ ) {
            mot[k] = *pmot;
            pmot = pmot + ( *ltype );
        }
    } else {
        return -1;
    }
    if ( type_id == H5T_FORTRAN_S1 ) {
        if ( ( datatype = H5Tcopy( type_id ) ) < 0 )
            return -1;
        if ( type_id == H5T_FORTRAN_S1 ) {
            if ( ( iret = H5Tset_size( datatype, *ltype ) ) < 0 )
                return -1;
            if ( ( iret = H5Tset_strpad( datatype, H5T_STR_SPACEPAD ) ) < 0 )
                return -1;
        }
    } else {
        datatype = type_id;
    }

    if ( ( dataspace = H5Screate_simple( 1, dimsf, NULL ) ) < 0 )
        return -1;
    if ( ( dataset = H5Dcreate( idfic, nomd, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT ) ) < 0 )
        return -1;
    if ( ( iret = H5Dwrite( dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, sv ) ) < 0 )
        return -1;
    if ( ( iret = H5Dclose( dataset ) ) < 0 )
        return -1;
    if ( ( iret = H5Sclose( dataspace ) ) < 0 )
        return -1;
    if ( type_id == H5T_FORTRAN_S1 ) {
        if ( ( iret = H5Tclose( datatype ) ) < 0 )
            return -1;
    }

    free( nomd );
    free( vtype );
    if ( type_id == H5T_FORTRAN_S1 ) {
        free( mot );
    }
#else
    CALL_UTMESS( "F", "FERMETUR_3" );
#endif
    return 0;
}
