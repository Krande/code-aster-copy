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

#ifndef ASTER_FORT_SUPERV_H_
#define ASTER_FORT_SUPERV_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALL_EXPASS( a ) CALLP( EXPASS, expass, a )
extern void DEFP( EXPASS, expass, ASTERINTEGER * );

#define CALL_EXECOP( a ) CALLP( EXECOP, execop, a )
extern void DEFP( EXECOP, execop, ASTERINTEGER * );

#define CALL_OP9999( a ) CALLP( OP9999, op9999, a )
extern void DEFP( OP9999, op9999, ASTERINTEGER * );

#define CALL_IMPERS() CALL0( IMPERS, impers )
extern void DEF0( IMPERS, impers );

#define CALL_IBMAIN() CALL0( IBMAIN, ibmain )
extern void DEF0( IBMAIN, ibmain );

#define CALL_CMD_CTXT_ENTER() CALL0( CMD_CTXT_ENTER, cmd_ctxt_enter )
extern void DEF0( CMD_CTXT_ENTER, cmd_ctxt_enter );

#define CALL_CMD_CTXT_EXIT() CALL0( CMD_CTXT_EXIT, cmd_ctxt_exit )
extern void DEF0( CMD_CTXT_EXIT, cmd_ctxt_exit );

#define CALL_ONERRF( a, b, c ) CALLSSP( ONERRF, onerrf, a, b, c )
extern void DEFSSP( ONERRF, onerrf, char *, STRING_SIZE, _OUT char *, STRING_SIZE,
                    _OUT ASTERINTEGER * );

#define CALL_GCNCON( a, b ) CALLSS( GCNCON, gcncon, a, b )
extern void DEFSS( GCNCON, gcncon, char *, STRING_SIZE, char *, STRING_SIZE );

#define CALLO_GNOMSD( a, b, c, d ) CALLOOPP( GNOMSD, gnomsd, a, b, c, d )
void DEFSSPP( GNOMSD, gnomsd, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *,
              ASTERINTEGER * );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_SUPERV_H */
#endif
