/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org             */
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

#ifndef ASTER_FORT_CALCUL_H_
#define ASTER_FORT_CALCUL_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALLO_ASASVE( a, b, c, d ) CALLOOOO( ASASVE, asasve, a, b, c, d )
void DEFSSSS( ASASVE, asasve, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ASCOVA( a, b, c, d, e, f, g, h )                                                     \
    CALLOOOOPOOO( ASCOVA, ascova, a, b, c, d, e, f, g, h )
void DEFSSSSPSSS( ASCOVA, ascova, const char *, STRING_SIZE, const char *, STRING_SIZE,
                  const char *, STRING_SIZE, const char *, STRING_SIZE, const ASTERDOUBLE *,
                  const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ASCAVC_WRAP( a, b, c, d, e, f, g )                                                   \
    CALLOOOOPOO( ASCAVC_WRAP, ascavc_wrap, a, b, c, d, e, f, g )
void DEFSSSSPSS( ASCAVC_WRAP, ascavc_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const char *, STRING_SIZE, const ASTERDOUBLE *,
                 const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_ASMATR( a, b, c, d, e, f, g, h, i )                                                   \
    CALLPSSSSSSPS( ASMATR, asmatr, a, b, c, d, e, f, g, h, i )
void DEFPSSSSSSPS( ASMATR, asmatr, ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                   STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                   STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                   STRING_SIZE );

#define CALL_ASSVEC( a, b, c, d, e, f, g)                                                          \
    CALLSSPSPSP( ASSVEC_WRAP, assvec_wrap, a, b, c, d, e, f, g )
void DEFSSPSPSP( ASSVEC_WRAP, assvec_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE,
                 ASTERDOUBLE *, const char *, STRING_SIZE,
                 ASTERINTEGER *);

#define CALL_ASSVECWITHMASK( a, b, c, d, e, f, g, h, i)                                          \
    CALLSSPSPSPSP( ASSVECWITHMASK, assvecwithmask, a, b, c, d, e, f, g, h, i )
void DEFSSPSPSPSP( ASSVECWITHMASK, assvecwithmask, const char *, STRING_SIZE,
                   const char *, STRING_SIZE,
                   ASTERINTEGER *, const char *, STRING_SIZE,
                   ASTERDOUBLE *, const char *, STRING_SIZE,
                   ASTERINTEGER *, const char *, STRING_SIZE, const ASTERLOGICAL * );

#define CALLO_AP_ASSEMBLY_VECTOR( a ) CALLO( AP_ASSEMBLY_VECTOR, ap_assembly_vector, a )
void DEFS( AP_ASSEMBLY_VECTOR, ap_assembly_vector, const char *, STRING_SIZE );

#define CALLO_CACHVC( a, b, c, d, e, f, g, h, i, j, k, l )                                         \
    CALLOOOOOOOOPPPP( CACHVC, cachvc, a, b, c, d, e, f, g, h, i, j, k, l )
void DEFSSSSSSSSPPPP( CACHVC, cachvc, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                      STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *,
                      ASTERINTEGER * );

#define CALLO_COMPSTRESSFIELD( a, b, c, d, e, f, g, h, i, j )                                      \
    CALLOOOOOOPPPP( COMPSTRESSFIELD, compstressfield, a, b, c, d, e, f, g, h, i, j)
void DEFSSSSSSPPPP( COMPSTRESSFIELD, compstressfield, const char *, STRING_SIZE,
                      const char *, STRING_SIZE,const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE,
                      ASTERLOGICAL *, ASTERLOGICAL *, ASTERINTEGER *, ASTERDOUBLE * );

#define CALLO_CONLAG( a, b ) CALLOP( CONLAG, conlag, a, b )
void DEFSP( CONLAG, conlag, const char *, STRING_SIZE, ASTERDOUBLE * );

#define CALLO_CORICHWRITE( a, b ) CALLOP( CORICHWRITE, corichwrite, a, b)
void DEFSP( CORICHWRITE, corichwrite,
            const char *, STRING_SIZE, ASTERINTEGER *);

#define CALLO_CRESOL_WRAP( a, b, c ) CALLOOO( CRESOL_WRAP, cresol_wrap, a, b, c )
void DEFSSS( CRESOL_WRAP, cresol_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALLO_MATRIX_FACTOR( a, b, c, d, e, f, g )                                                 \
    CALLOOPOOPP( MATRIX_FACTOR, matrix_factor, a, b, c, d, e, f, g )
void DEFSSPSSPP( MATRIX_FACTOR, matrix_factor, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, ASTERINTEGER * );

#define CALLO_NMDOCH_WRAP( a, b, c, d ) CALLOPOO( NMDOCH_WRAP, nmdoch_wrap, a, b, c, d )
void DEFSPSS( NMDOCH_WRAP, nmdoch_wrap, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NTDOCH_WRAP( a, b, c, d ) CALLOPOO( NTDOCH_WRAP, ntdoch_wrap, a, b, c, d )
void DEFSPSS( NTDOCH_WRAP, ntdoch_wrap, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NUMCIMA( a, b, c, d )                                                      \
    CALLOOOO( NUMCIMA, numcima, a, b, c, d)
void DEFSSSS( NUMCIMA, numcima, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NUMERO_WRAP( a, b, c, d, e, f )                                                      \
    CALLOOOOOO( NUMERO_WRAP, numero_wrap, a, b, c, d, e, f )
void DEFSSSSSS( NUMERO_WRAP, numero_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE );

#define CALLO_NUME_DDL_MATR( a, b, c ) CALLOOP( NUME_DDL_MATR, nume_ddl_matr, a, b, c )
void DEFSSP( NUME_DDL_MATR, nume_ddl_matr, const char *, STRING_SIZE, const char *, STRING_SIZE,
             ASTERINTEGER * );

#define CALLO_RESOUD( a, b, c, d, e, f, g, h, i, j, k, l, m, n )                                   \
    CALLOOOOPOOOPPOPPP( RESOUD, resoud, a, b, c, d, e, f, g, h, i, j, k, l, m, n )
void DEFSSSSPSSSPPSPPP( RESOUD, resoud, const char *, STRING_SIZE, const char *, STRING_SIZE,
                        const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *,
                        const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                        STRING_SIZE, ASTERDOUBLE *, ASTERDOUBLE *, const char *, STRING_SIZE,
                        ASTERLOGICAL *, ASTERINTEGER *, ASTERINTEGER * );

#define CALLO_VECHME_WRAP( a, b, c, d, e, f, g, h, i, l )                                          \
    CALLOOOOPOOOOO( VECHME_WRAP, vechme_wrap, a, b, c, d, e, f, g, h, i, l )
void DEFSSSSPSSSSS( VECHME_WRAP, vechme_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE, const char *, STRING_SIZE, const ASTERDOUBLE *,
                    const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                    const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VEDIME( a, b, c, d, e, f ) CALLOOOPOO( VEDIME, vedime, a, b, c, d, e, f )
void DEFSSSPSS( VEDIME, vedime, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                STRING_SIZE, ASTERDOUBLE *, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VEBTLA( a, b, c, d, e, f, g ) CALLOOOOOOO( VEBTLA, vebtla, a, b, c, d, e, f, g )
void DEFSSSSSSS( VEBTLA, vebtla, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VEBUME( a, b, c, d, e, f ) CALLOOOOPO( VEBUME, vebume, a, b, c, d, e, f )
void DEFSSSSPS( VEBUME, vebume, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE,
                const ASTERDOUBLE *, const char *, STRING_SIZE );

#define CALLO_VELAME( a, b, c, d, e ) CALLOOOOO( VELAME, velame, a, b, c, d, e )
void DEFSSSSS( VELAME, velame, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VRCINS_WRAP( a, b, c, d, e, f )                                                      \
    CALLOOOPOO( VRCINS_WRAP, vrcins_wrap, a, b, c, d, e, f )
void DEFSSSPSS( VRCINS_WRAP, vrcins_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE,
                const char *, STRING_SIZE );

#define CALLO_VRCREF( a, b, c, d ) CALLOOOO( VRCREF, vrcref, a, b, c, d )
void DEFSSSS( VRCREF, vrcref, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VTCREB_WRAP( a, b, c, d ) CALLOOOO( VTCREB_WRAP, vtcreb_wrap, a, b, c, d )
void DEFSSSS( VTCREB_WRAP, vtcreb_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NMDOCC( a, b, c, d, e, f, g ) CALLOOPPOOP( NMDOCC, nmdocc, a, b, c, d, e, f, g )
void DEFSSPPSSP( NMDOCC, nmdocc, const char *, STRING_SIZE, const char *, STRING_SIZE,
                ASTERLOGICAL *, ASTERLOGICAL *, const char *, STRING_SIZE,
                const char *, STRING_SIZE, ASTERLOGICAL * );

#define CALLO_NMDOCR( a, b, c, d ) CALLOOPO( NMDOCR, nmdocr, a, b, c, d )
void DEFSSPS( NMDOCR, nmdocr, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERLOGICAL *,
            const char *, STRING_SIZE );

#define CALLO_NMDOCM( a, b, c ) CALLOOO( NMDOCM, nmdocm, a, b, c )
void DEFSSS( NMDOCM, nmdocm, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALLO_AFVARC( a, b, c ) CALLOOO( AFVARC, afvarc, a, b, c )
void DEFSSS( AFVARC, afvarc, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALLO_CMTREF( a, b ) CALLOO( CMTREF, cmtref, a, b )
void DEFSS( CMTREF, cmtref, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_CESVAR( a, b, c, d ) CALLOOOO( CESVAR, cesvar, a, b, c, d )
void DEFSSSS( CESVAR, cesvar, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE);

#define CALLO_ALCHML( a, b, c, d, e, f, g ) CALLOOOOOPO( ALCHML, alchml, a, b, c, d, e, f, g )
void DEFSSSSSPS( ALCHML, alchml, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE );

#define CALLO_RSAGSD( a, b ) CALLOP( RSAGSD, rsagsd, a, b)
void DEFSP( RSAGSD, rsagsd, const char *, STRING_SIZE, ASTERINTEGER *);

#define CALLO_CALCUL( a, b, c, d, e, f, g, h, i, j, k )                                           \
    CALLSSSPSSPSSSS( CALCUL_CWRAP, calcul_cwrap, a, b, c, d, e, f, g, h, i, j, k )
void DEFSSSPSSPSSSS( CALCUL_CWRAP, calcul_cwrap, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE,
              ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
              ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_CHECKSUPERELEMENT( a, b )                                                           \
    CALLOO( CHECKSUPERELEMENT, checksuperelement, a, b )
void DEFSS( CHECKSUPERELEMENT, checksuperelement, const char *, STRING_SIZE,
            const char *, STRING_SIZE);

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_CALCUL_H */
#endif
