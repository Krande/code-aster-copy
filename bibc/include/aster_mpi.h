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

/* person_in_charge: mathieu.courtois at edf.fr */

#ifndef ASTER_MPI_H_
#define ASTER_MPI_H_

#include "aster.h"
#ifdef ASTER_HAVE_MPI
#include "mpi.h"
#else
#define MPI_COMM_WORLD 0
#define MPI_COMM_NULL -12345
#define MPI_Fint int
#define MPI_Comm int
#define MPI_Comm_c2f( a ) a
#define MPI_Comm_f2c( a ) a
#define MPI_Request int
#define MPI_Request_c2f( a ) a
#define MPI_Request_f2c( a ) a
#define MPI_Op int
#define MPI_Op_c2f( a ) a
#define MPI_Op_f2c( a ) a
#define MPI_Datatype int
#define MPI_DOUBLE_PRECISION 1
#define MPI_DOUBLE_COMPLEX 2
#define MPI_INTEGER8 3
#define MPI_INTEGER4 4
#define MPI_CHAR 5
#endif

/*
 *  Structure to store the communicator tree.
 * - id : the current communicator
 * - parent : its parent aster_comm_t communicator
 * - level: 0 for MPI_COMM_WORLD, +1 at each split
 * - childs: child communicators
 * - nbchild: number of childs
 * - name: for nicer print and debug
 */
#define MAX_CHILDS 2000
#define NAME_LENGTH 16

#ifdef __cplusplus
extern "C" {
#endif
typedef struct aster_comm_t aster_comm_t;

struct aster_comm_t {
    MPI_Comm id;
    aster_comm_t *parent;
    int level;
    aster_comm_t *childs[MAX_CHILDS];
    int nbchild;
    char name[NAME_LENGTH];
};

// #define ASTER_DEBUG_UNITTEST_MPI

/*
 *   PUBLIC FUNCTIONS
 *
 */

extern void aster_mpi_init( const MPI_Fint );

extern aster_comm_t *aster_get_comm_world();
extern aster_comm_t *aster_get_current_comm();
extern void aster_set_current_comm( aster_comm_t * );
extern void aster_get_mpi_info( aster_comm_t *, int *, int * );
extern aster_comm_t *aster_split_comm( aster_comm_t *, int, int, char * );
extern void aster_free_comm( aster_comm_t * );
extern int aster_set_mpi_barrier( aster_comm_t * );
extern int aster_mpi_bcast( void *, int, MPI_Datatype, int, aster_comm_t * );
extern int aster_mpi_allreduce( void *, void *, int, MPI_Datatype, MPI_Op, aster_comm_t * );
extern int aster_mpi_gather( void *, int, MPI_Datatype, void *, int, MPI_Datatype, int,
                             aster_comm_t * );
extern int aster_mpi_gatherv( void *, int, MPI_Datatype, void *, int *, int *, MPI_Datatype, int,
                              aster_comm_t * );
extern int aster_mpi_allgather( void *, int, MPI_Datatype, void *, int, MPI_Datatype,
                                aster_comm_t * );
extern int aster_mpi_allgatherv( void *, int, MPI_Datatype, void *, int *, int *, MPI_Datatype,
                                 aster_comm_t * );
extern int aster_mpi_send( void *, int, MPI_Datatype, int, int, aster_comm_t * );
extern int aster_mpi_recv( void *, int, MPI_Datatype, int, int, aster_comm_t * );

extern int aster_mpi_sendrecv( void *, int, MPI_Datatype, int, int, void *, int, MPI_Datatype, int,
                               int, aster_comm_t * );

extern void DEFSP( ASMPI_COMM, asmpi_comm, const char *, STRING_SIZE, MPI_Fint * );
extern void DEFPPPSP( ASMPI_SPLIT_COMM, asmpi_split_comm, MPI_Fint *, MPI_Fint *, MPI_Fint *,
                      const char *, STRING_SIZE, MPI_Fint * );
extern void DEFPPP( ASMPI_INFO_WRAP, asmpi_info_wrap, MPI_Fint *, MPI_Fint *, MPI_Fint * );
extern void terminate( void );

#define CALL_ASABRT( a ) CALLP( ASABRT, asabrt, a )
extern void DEFP( ASABRT, asabrt, _IN ASTERINTEGER * );

#define CALLO_ASMPI_COMM( a, b ) CALLOP( ASMPI_COMM, asmpi_comm, a, b )
#define CALL_ASMPI_SPLIT_COMM( a, b, c, d, e )                                                     \
    CALLPPPOP( ASMPI_SPLIT_COMM, asmpi_split_comm, a, b, c, d, e );
#define CALL_ASMPI_INFO( a, b, c ) CALLPPP( ASMPI_INFO_WRAP, asmpi_info_wrap, a, b, c )

/* AS_MPICHECK checks the return code of an MPI function */
#define AS_MPICHECK( code )                                                                        \
    do {                                                                                           \
        int ret = ( code );                                                                        \
        if ( ( ret ) == MPI_ERR_COMM ) {                                                           \
            DEBUG_LOC;                                                                             \
            DBGV( "Invalid communicator: %s", #code );                                             \
            INTERRUPT( 17 );                                                                       \
        } else if ( ( ret ) == MPI_ERR_TYPE ) {                                                    \
            DEBUG_LOC;                                                                             \
            DBGV( "Invalid datatype argument: %s", #code );                                        \
            INTERRUPT( 17 );                                                                       \
        } else if ( ( ret ) == MPI_ERR_COUNT ) {                                                   \
            DEBUG_LOC;                                                                             \
            DBGV( "Invalid count argument: %s", #code );                                           \
            INTERRUPT( 17 );                                                                       \
        } else if ( ( ret ) == MPI_ERR_TAG ) {                                                     \
            DEBUG_LOC;                                                                             \
            DBGV( "Invalid tag argument: %s", #code );                                             \
            INTERRUPT( 17 );                                                                       \
        } else if ( ( ret ) == MPI_ERR_RANK ) {                                                    \
            DEBUG_LOC;                                                                             \
            DBGV( "Invalid source or destination rank: %s", #code );                               \
            INTERRUPT( 17 );                                                                       \
        } else if ( ( ret ) != MPI_SUCCESS ) {                                                     \
            DEBUG_LOC;                                                                             \
            DBGV( "Unknown error: %s", #code );                                                    \
            INTERRUPT( 17 );                                                                       \
        }                                                                                          \
    } while ( 0 );

/*
 *   PRIVATE FUNCTIONS
 *
 */
extern void errhdlr_func( MPI_Comm *, int *, ... );
extern aster_comm_t *_search_id( aster_comm_t *, MPI_Comm * );
aster_comm_t *get_node_by_id( MPI_Comm * );
#ifdef ASTER_DEBUG_UNITTEST_MPI
extern void _unittest_aster_mpi();
#endif

#ifdef __cplusplus
}
#endif

#endif
