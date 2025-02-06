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

#ifndef ASTER_FORT_DS_H_
#define ASTER_FORT_DS_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALLO_ALCART( a, b, c, d ) CALLOOOO( ALCART, alcart, a, b, c, d )
void DEFSSSS( ALCART, alcart, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_CELCES_WRAP( a, b, c ) CALLOOO( CELCES_WRAP, celces_wrap, a, b, c )
void DEFSSS( CELCES_WRAP, celces_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALL_CARCES( a, b, c, d, e, f, g ) CALLOOOOOOP( CARCES, carces, a, b, c, d, e, f, g )
void DEFSSSSSSP( CARCES, carces, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, ASTERINTEGER * );

#define CALLO_CESCEL_WRAP( a, b, c, d, e, f, g, h, i, j )                                          \
    CALLOOOOOPOOOP( CESCEL_WRAP, cescel_wrap, a, b, c, d, e, f, g, h, i, j )
void DEFSSSSSPSSSP( CESCEL_WRAP, cescel_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_CESCNS( a, b, c, d, e, f ) CALLOOOOOP( CESCNS, cescns, a, b, c, d, e, f )
void DEFSSSSSP( CESCNS, cescns, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALL_CESCRE_WRAP( a, b, c, d, e, f, g, h, i, j, k )                                        \
    CALLSSSSSPSPPPP( CESCRE_WRAP, cescre_wrap, a, b, c, d, e, f, g, h, i, j, k )
#define CALLO_CESCRE_WRAP( a, b, c, d, e, f, g, h, i, j, k )                                       \
    CALLOOOOOPOPPPP( CESCRE_WRAP, cescre_wrap, a, b, c, d, e, f, g, h, i, j, k )
void DEFSSSSSPSPPPP( CESCRE_WRAP, cescre_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                     const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                     STRING_SIZE, ASTERINTEGER *, const char *, STRING_SIZE, ASTERINTEGER *,
                     ASTERINTEGER *, ASTERINTEGER *, ASTERLOGICAL * );

#define CALL_CESRED_WRAP( a, b, c, d, e, g, h )                                                    \
    CALLSPPPSSS( CESRED_WRAP, cesred_wrap, a, b, c, d, e, g, h )
#define CALLO_CESRED_WRAP( a, b, c, d, e, g, h )                                                   \
    CALLOPPPOOO( CESRED_WRAP, cesred_wrap, a, b, c, d, e, g, h )
void DEFSPPPSSS( CESRED_WRAP, cesred_wrap, const char *, STRING_SIZE, ASTERINTEGER *,
                 ASTERINTEGER *, ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_CNOCNS_WRAP( a, b, c ) CALLOOO( CNOCNS_WRAP, cnocns_wrap, a, b, c )
void DEFSSS( CNOCNS_WRAP, cnocns_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALL_CNSCRE_WRAP( a, b, c, d, e, g, h )                                                    \
    CALLSSPSSSP( CNSCRE_WRAP, cnscre_wrap, a, b, c, d, e, g, h )
#define CALLO_CNSCRE_WRAP( a, b, c, d, e, g, h )                                                   \
    CALLOOPOOOP( CNSCRE_WRAP, cnscre_wrap, a, b, c, d, e, g, h )
void DEFSSPSSSP( CNSCRE_WRAP, cnscre_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, ASTERLOGICAL * );

#define CALLO_CNSCNO_WRAP( a, b, c, d, e, g, h )                                                   \
    CALLOOOOOOP( CNSCNO_WRAP, cnscno_wrap, a, b, c, d, e, g, h )
void DEFSSSSSSP( CNSCNO_WRAP, cnscno_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, ASTERINTEGER * );

#define CALL_CLEAN_JEVEUX_MEMORY() CALL0( CLEANJEVEUXMEMORY, cleanjeveuxmemory )
void DEF0( CLEANJEVEUXMEMORY, cleanjeveuxmemory );

#define CALL_DETMAT() CALL0( DETMAT, detmat )
extern void DEF0( DETMAT, detmat );

#define CALL_DETMATRIX( a ) CALLS( DETMATRIX, detmatrix, a )
#define CALLO_DETMATRIX( a ) CALLO( DETMATRIX, detmatrix, a )
extern void DEFS( DETMATRIX, detmatrix, const char *, STRING_SIZE );

#define CALLO_DELETE_MATRIX( a, b ) CALLOO( DELETE_MATRIX, delete_matrix, a, b )
void DEFSS( DELETE_MATRIX, delete_matrix, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_DELETE_CACHED_OBJECTS() CALL0( DELETECACHEDOBJECTS, deletecachedobjects )
void DEF0( DELETECACHEDOBJECTS, deletecachedobjects );

#define CALL_DELETE_TEMPORARY_OBJECTS() CALL0( DELETETEMPORARYOBJECTS, deletetemporaryobjects )
void DEF0( DELETETEMPORARYOBJECTS, deletetemporaryobjects );

#define CALLO_EXISD( a, b, c ) CALLOOP( EXISD, exisd, a, b, c )
void DEFSSP( EXISD, exisd, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const ASTERINTEGER * );

#define CALLO_MATR_ASSE_SET_VALUES( a, b, c, d, e )                                                \
    CALLOPPPP( MATR_ASSE_SET_VALUES, matr_asse_set_values, a, b, c, d, e )
extern void DEFSPPPP( MATR_ASSE_SET_VALUES, matr_asse_set_values, const char *, STRING_SIZE,
                      const ASTERINTEGER *, const ASTERINTEGER *, const ASTERINTEGER *,
                      const ASTERDOUBLE * );

#define CALLO_MATR_ASSE_SCALE( a, b, c ) CALLOPP( CALLO_MATR_ASSE_SCALE, matr_asse_scale, a, b, c )
extern void DEFSPP( MATR_ASSE_SCALE, matr_asse_scale, const char *, STRING_SIZE,
                    const ASTERDOUBLE *, const ASTERDOUBLE * );

#define CALLO_MATR_ASSE_TRANSPOSE( a ) CALLO( MATR_ASSE_TRANSPOSE, matr_asse_transpose, a )
void DEFS( MATR_ASSE_TRANSPOSE, matr_asse_transpose, const char *, STRING_SIZE );

#define CALLO_MATR_ASSE_TRANSPOSE_CONJUGATE( a )                                                   \
    CALLO( MATR_ASSE_TRANSPOSE_CONJUGATE, matr_asse_transpose_conjugate, a )
void DEFS( MATR_ASSE_TRANSPOSE_CONJUGATE, matr_asse_transpose_conjugate, const char *,
           STRING_SIZE );

#define CALLO_MATR_ASSE_PRINT( a, b, c ) CALLOPO( MATIMP, matimp, a, b, c )
void DEFSPS( MATIMP, matimp, const char *, STRING_SIZE, const ASTERINTEGER *, const char *,
             STRING_SIZE );

#define CALLO_POSDDL( a, b, c, d, e, f ) CALLOOOOPP( POSDDL, posddl, a, b, c, d, e, f )
void DEFSSSSPP( POSDDL, posddl, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                STRING_SIZE, const char *, STRING_SIZE, const ASTERINTEGER *,
                const ASTERINTEGER * );

#define CALL_EXLIM1( a, b, c, d, e ) CALLPPOOO( EXLIM1, exlim1, a, b, c, d, e )
void DEFPPSSS( EXLIM1, exlim1, const ASTERINTEGER *, const ASTERINTEGER *, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_EXLIM2( a, b, c, d, e ) CALLPPOOO( EXLIM2, exlim2, a, b, c, d, e )
void DEFPPSSS( EXLIM2, exlim2, const ASTERINTEGER *, const ASTERINTEGER *, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_DETRSD( a, b ) CALLOO( DETRSD, detrsd, a, b )
extern void DEFSS( DETRSD, detrsd, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_LGTLGR( a, b, c ) CALLSSS( LGTLGR, lgtlgr, a, b, c )
#define CALLO_LGTLGR( a, b, c ) CALLOOO( LGTLGR, lgtlgr, a, b, c )
extern void DEFSSS( LGTLGR, lgtlgr, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE );

#define CALLO_COPISD( a, b, c, d ) CALLOOOO( COPISD, copisd, a, b, c, d )
extern void DEFSSSS( COPISD, copisd, const char *, STRING_SIZE, const char *, STRING_SIZE,
                     const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ADALIG_WRAP( a ) CALLO( ADALIG_WRAP, adalig_wrap, a )
extern void DEFS( ADALIG_WRAP, adalig_wrap, const char *, STRING_SIZE );

#define CALLO_INITEL( a, b ) CALLOP( INITEL, initel, a, b )
extern void DEFSP( INITEL, initel, const char *, STRING_SIZE, ASTERLOGICAL * );

#define CALLO_CORMGI( a, b ) CALLOO( CORMGI, cormgi, a, b )
extern void DEFSS( CORMGI, cormgi, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_RGNDAS_WRAP( a, b, c, d, e, f )                                                      \
    CALLOPOOOO( RGNDAS_WRAP, rgndas_wrap, a, b, c, d, e, f )
void DEFSPSSSS( RGNDAS_WRAP, rgndas_wrap, const char *, STRING_SIZE, const ASTERINTEGER *,
                const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE );

#define CALL_NUMEDDL_GET_COMPONENTS( a, b, c, d, e, f )                                            \
    CALLSSPPSP( NUMEDDL_GET_COMPONENTS, numeddl_get_components, a, b, c, d, e, f )
extern void DEFSSPPSP( NUMEDDL_GET_COMPONENTS, numeddl_get_components, const char *, STRING_SIZE,
                       const char *, STRING_SIZE, const ASTERINTEGER *, ASTERINTEGER *, char *,
                       STRING_SIZE, const ASTERINTEGER * );

#define CALLO_NUMEDDL_GET_COMPONENT_NAME( a, b, c )                                                \
    CALLOPO( NUMEDDL_GET_COMPONENT_NAME, numeddl_get_component_name, a, b, c )
extern void DEFSPS( NUMEDDL_GET_COMPONENT_NAME, numeddl_get_component_name, const char *,
                    STRING_SIZE, const ASTERINTEGER *, const char *, STRING_SIZE );

#define CALLO_NOCARTC( a, b, c, d, e, f, g, h, i )                                                 \
    CALLOPPOOPOPO( NOCART_C, nocart_c, a, b, c, d, e, f, g, h, i )
void DEFSPPSSPSPS( NOCART_C, nocart_c, const char *, STRING_SIZE, const ASTERINTEGER *,
                   const ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
                   const ASTERINTEGER *, const char *, STRING_SIZE, const ASTERINTEGER *,
                   const char *, STRING_SIZE );

#ifdef ASTER_STRLEN_AT_END
#define CALL_RSACCH( nomsd, numch, nomch, nbord, liord, nbcmp, liscmp )                            \
    F_FUNC( RSACCH, rsacch )                                                                       \
    ( nomsd, numch, nomch, nbord, liord, nbcmp, liscmp, strlen( nomsd ), 16, 8 )
#else
#define CALL_RSACCH( nomsd, numch, nomch, nbord, liord, nbcmp, liscmp )                            \
    F_FUNC( RSACCH, rsacch )                                                                       \
    ( nomsd, strlen( nomsd ), numch, nomch, 16, nbord, liord, nbcmp, liscmp, 8 )
#endif
extern void DEFSPSPPPS( RSACCH, rsacch, char *, STRING_SIZE, ASTERINTEGER *, char *, STRING_SIZE,
                        ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, char *, STRING_SIZE );

#define CALL_RSACPA( nomsd, numva, icode, nomva, ctype, ival, rval, kval, ier )                    \
    CALLSPPSPPPSP( RSACPA, rsacpa, nomsd, numva, icode, nomva, ctype, ival, rval, kval, ier )
extern void DEFSPPSPPPSP( RSACPA, rsacpa, char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *,
                          char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *, ASTERDOUBLE *,
                          char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_RSADPA_ZK8_WRAP( a, b, c, d, e )                                                     \
    CALLOPOOO( RSADPA_ZK8_WRAP, rsadpa_zk8_wrap, a, b, c, d, e )
extern void DEFSPSSS( RSADPA_ZK8_WRAP, rsadpa_zk8_wrap, const char *, STRING_SIZE, ASTERINTEGER *,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                      STRING_SIZE );

#define CALLO_RSADPA_ZK16_WRAP( a, b, c, d, e )                                                    \
    CALLOPOOO( RSADPA_ZK16_WRAP, rsadpa_zk16_wrap, a, b, c, d, e )
extern void DEFSPSSS( RSADPA_ZK16_WRAP, rsadpa_zk16_wrap, const char *, STRING_SIZE, ASTERINTEGER *,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                      STRING_SIZE );

#define CALLO_RSADPA_ZR_WRAP( a, b, c, d ) CALLOPPO( RSADPA_ZR_WRAP, rsadpa_zr_wrap, a, b, c, d )
extern void DEFSPPS( RSADPA_ZR_WRAP, rsadpa_zr_wrap, const char *, STRING_SIZE, ASTERINTEGER *,
                     ASTERDOUBLE *, const char *, STRING_SIZE );

#define CALLO_RSADPA_ZK24_WRAP( a, b, c, d, e )                                                    \
    CALLOPOOO( RSADPA_ZK24_WRAP, rsadpa_zk24_wrap, a, b, c, d, e )
extern void DEFSPSSS( RSADPA_ZK24_WRAP, rsadpa_zk24_wrap, const char *, STRING_SIZE, ASTERINTEGER *,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                      STRING_SIZE );

#define CALLO_RSAGSD( a, b ) CALLOP( RSAGSD, rsagsd, a, b )
void DEFSP( RSAGSD, rsagsd, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_RSCRSD( a, b, c, d ) CALLOOOP( RSCRSD, rscrsd, a, b, c, d )
void DEFSSSP( RSCRSD, rscrsd, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const ASTERINTEGER * );

#define CALLO_RSRUSD( a, b ) CALLOP( RSRUSD, rsrusd, a, b )
void DEFSP( RSRUSD, rsrusd, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_RSINFO( a, b ) CALLOP( RSINFO, rsinfo, a, b )
extern void DEFSP( RSINFO, rsinfo, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_VTGPLD( a, b, c, d, e, f ) CALLOPOOOO( VTGPLD, vtgpld, a, b, c, d, e, f )
extern void DEFSPSSSS( VTGPLD, vtgpld, const char *, STRING_SIZE, ASTERDOUBLE *, const char *,
                       STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                       const char *, STRING_SIZE );

#define CALLO_APLCPGN( a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p )                            \
    CALLOOOOPPPPPPPPPOOO( APLCPGN, aplcpgn, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p )
extern void DEFSSSSPPPPPPPPPSSS( APLCPGN, aplcpgn, const char *, STRING_SIZE, const char *,
                                 STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                                 ASTERDOUBLE *, ASTERDOUBLE *, ASTERINTEGER *, ASTERINTEGER *,
                                 ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *,
                                 ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                                 STRING_SIZE, const char *, STRING_SIZE );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_DS_H */
#endif
