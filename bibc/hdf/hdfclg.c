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
/ Fermeture d'un groupe HDF, renvoie une erreur si le groupe ne peut être fermé
/  Paramètres :
/   - in idgrp : identificateur du groupe (hid_t)
/  Résultats :
/     0 = fermeture OK, -1 sinon (long)
/-----------------------------------------------------------------------------*/
#ifdef ASTER_HAVE_HDF5
#include <hdf5.h>
#else
typedef int hid_t;
#endif

ASTERINTEGER DEFP( HDFCLG, hdfclg, hid_t *idg ) {
#ifdef ASTER_HAVE_HDF5
    hid_t idgrp;
    herr_t icode;
    idgrp = (hid_t)*idg;
    if ( ( icode = H5Gclose( idgrp ) ) < 0 )
        return -1;
#else
    CALL_UTMESS( "F", "FERMETUR_3" );
#endif
    return 0;
}
