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

#ifndef ASTER_FORT_MESH_H_
#define ASTER_FORT_MESH_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALLO_ADDGROUPNODE( a, b ) CALLOP( ADDGROUPNODE, addgroupnode, a, b )
void DEFSP( ADDGROUPNODE, addgroupnode, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_ADDGROUPELEM( a, b ) CALLOP( ADDGROUPELEM, addgroupelem, a, b )
void DEFSP( ADDGROUPELEM, addgroupelem, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_ADDGRPMA( a, b, c, d, e ) CALLOOPPP( ADDGRPMA, addgrpma, a, b, c, d, e )
void DEFSSPPP( ADDGRPMA, addgrpma, const char *, STRING_SIZE, const char *, STRING_SIZE,
               const ASTERINTEGER *, ASTERINTEGER *, ASTERLOGICAL * );

#define CALLO_ADDGRPNO( a, b, c, d, e ) CALLOOPPP( ADDGRPNO, addgrpno, a, b, c, d, e )
void DEFSSPPP( ADDGRPNO, addgrpno, const char *, STRING_SIZE, const char *, STRING_SIZE,
               const ASTERINTEGER *, ASTERINTEGER *, ASTERLOGICAL * );

#define CALLO_CARGEO( a ) CALLO( CARGEO, cargeo, a )
void DEFS( CARGEO, cargeo, const char *, STRING_SIZE );

#define CALLO_INFOMA( a, b ) CALLOP( INFOMA, infoma, a, b )
void DEFSP( INFOMA, infoma, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_CHCKMA( a, b ) CALLOP( CHCKMA, chckma, a, b )
void DEFSP( CHCKMA, chckma, const char *, STRING_SIZE, ASTERDOUBLE * );

#define CALLO_LRMJOI_WRAP( a, b ) CALLOO( LRMJOI_WRAP, lrmjoi_wrap, a, b )
void DEFSS( LRMJOI_WRAP, lrmjoi_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_LRM_CLEAN_JOINT( a, b ) CALLOP( LRM_CLEAN_JOINT, lrm_clean_joint, a, b )
void DEFSP( LRM_CLEAN_JOINT, lrm_clean_joint, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALL_MDNOMA( a, b, c, d ) CALLSPSP( MDNOMA, mdnoma, a, b, c, d )
extern void DEFSPSP( MDNOMA, mdnoma, char *, STRING_SIZE, ASTERINTEGER *, char *, STRING_SIZE,
                     ASTERINTEGER * );

#define CALL_MDNOCH( a, b, c, d, e, f, g ) CALLSPPSSSP( MDNOCH, mdnoch, a, b, c, d, e, f, g )
extern void DEFSPPSSSP( MDNOCH, mdnoch, char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *, char *,
                        STRING_SIZE, char *, STRING_SIZE, char *, STRING_SIZE, ASTERINTEGER * );

#define CALL_GET_MED_TYPES( a, b ) CALLSS( GET_MED_TYPES, get_med_types, a, b )
#define CALLO_GET_MED_TYPES( a, b ) CALLOO( GET_MED_TYPES, get_med_types, a, b )
extern void DEFSS( GET_MED_TYPES, get_med_types, const char *, STRING_SIZE, const char *,
                   STRING_SIZE );

#define CALL_GET_MED_CONNECTIVITY( a, b ) CALLSS( GET_MED_CONNECTIVITY, get_med_connectivity, a, b )
#define CALLO_GET_MED_CONNECTIVITY( a, b )                                                         \
    CALLOO( GET_MED_CONNECTIVITY, get_med_connectivity, a, b )
extern void DEFSS( GET_MED_CONNECTIVITY, get_med_connectivity, const char *, STRING_SIZE,
                   const char *, STRING_SIZE );

#define CALLO_FILL_SPARSE_NAMED( a, b, c, d )                                                      \
    CALLOOOO( FILL_SPARSE_NAMED, fill_sparse_named, a, b, c, d )
extern void DEFSSSS( FILL_SPARSE_NAMED, fill_sparse_named, const char *, STRING_SIZE, const char *,
                     STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_CNCINV( a, b, c, d, e ) CALLSPPSS( CNCINV, cncinv, a, b, c, d, e )
#define CALLO_CNCINV( a, b, c, d, e ) CALLOPPOO( CNCINV, cncinv, a, b, c, d, e )
extern void DEFSPPSS( CNVINV, cncinv, const char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *,
                      const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_CHVENO( a, b, c ) CALLSSS( CHVENO, chveno, a, b, c )
#define CALLO_CHVENO( a, b, c ) CALLOOO( CHVENO, chveno, a, b, c )
extern void DEFSSS( CHVENO, chveno, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE );

#define CALL_CHECKNORMALS( a, b, c ) CALLSSS( CHECKNORMALS, checknormals, a, b, c )
#define CALLO_CHECKNORMALS( a, b, c ) CALLOOO( CHECKNORMALS, checknormals, a, b, c )
extern void DEFSSS( CHECKNORMALS, checknormals, const char *, STRING_SIZE, const char *,
                    STRING_SIZE, const char *, STRING_SIZE );

#define CALL_CNVOIS( a, b, c, d, e, f, g ) CALLOPOPPPO( CNVOIS, cnvois, a, b, c, d, e, f, g )
extern void DEFSPSPPPS( CNVOIS, cnvois, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                        STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, const char *,
                        STRING_SIZE );

#define CALL_CMBQBQ( a, b, c, d ) CALLOOPP( CMBQBQ, cmbqbq, a, b, c, d )
extern void DEFSSPP( CMBQBQ, cmbqbq, const char *, STRING_SIZE, const char *, STRING_SIZE,
                     ASTERINTEGER *, ASTERINTEGER * );

#define CALL_AJELLT( a, b, c, d, e, f, g, h, i )                                                   \
    CALLSSPSSSSPS( AJELLT, ajellt, a, b, c, d, e, f, g, h, i )
#define CALLO_AJELLT( a, b, c, d, e, f, g, h, i )                                                  \
    CALLOOPOOOOPO( AJELLT, ajellt, a, b, c, d, e, f, g, h, i )
extern void DEFSSPSSSSPS( AJELLT, ajellt, const char *, STRING_SIZE, const char *, STRING_SIZE,
                          const ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                          STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                          const ASTERINTEGER *, const char *, STRING_SIZE );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_MESH_H */
#endif
