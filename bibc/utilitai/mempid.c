/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
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

/* person_in_charge: j-pierre.lefebvre at edf.fr */
#include "aster.h"

#if defined ASTER_PLATFORM_LINUX
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#elif defined ASTER_PLATFORM_FREEBSD
#include <err.h>
#include <fcntl.h>
#include <kvm.h>
#include <sys/param.h>
#include <sys/sysctl.h>
#include <sys/user.h>

#elif defined ASTER_PLATFORM_MINGW
#include <windows.h>
#
#include <psapi.h>
#endif

/*
 * This function returns the memory consumption of the current process.
 * It returns VmSize and VmPeak.
 */

/*
 * OS X does not support retrieving memory consumptions through /proc or kvm library
 */
#ifdef ASTER_PLATFORM_DARWIN
#define ASTER_DISABLE_MEMORY_STATS
#endif

/*
 * If ASTER_DISABLE_MEMORY_STATS, it always return 0
 */
#ifdef ASTER_DISABLE_MEMORY_STATS
ASTERINTEGER DEFP( MEMPID, mempid, ASTERINTEGER *val ) {
    val[0] = 0;
    val[1] = 0;
    return 0;
}

#elif defined ASTER_PLATFORM_LINUX
ASTERINTEGER get_memory_usage_cgroup_v2( pid_t, ASTERINTEGER * );
ASTERINTEGER get_memory_usage_proc( pid_t, ASTERINTEGER * );

ASTERINTEGER DEFP( MEMPID, mempid, ASTERINTEGER *val ) {
    pid_t pid = getpid();
#if defined ASTER_ENABLE_CGROUP_STATS
#define get_memory_stats get_memory_usage_cgroup_v2
#elif defined ASTER_ENABLE_PROC_STATUS
#define get_memory_stats get_memory_usage_proc
#else
#error ASTER_ENABLE_[CGROUP_STATS or PROC_STATUS] is required, or ASTER_DISABLE_MEMORY_STATS
#endif
    return get_memory_stats( pid, val );
}

#elif defined ASTER_PLATFORM_FREEBSD
/*
 * FreeBSD and some others without /proc ?
 */
ASTERINTEGER DEFP( MEMPID, mempid, ASTERINTEGER *val ) {
    pid_t getpid( void );

    pid_t pid = getpid();

#define B2K( x ) ( ( x ) >> 10 ) /* bytes to kbytes */

    char errbuf[_POSIX2_LINE_MAX];
    struct kinfo_proc *kp;
    kvm_t *kd;
    int count;
    kd = kvm_openfiles( NULL, "/dev/null", NULL, O_RDONLY, errbuf );
    if ( kd == NULL )
        errx( 1, "kvm_openfiles: %s", errbuf );

    kp = kvm_getprocs( kd, KERN_PROC_PID, pid, &count );
    if ( kp == NULL ) {
        (void)fprintf( stderr, "kvm_getprocs: %s", kvm_geterr( kd ) );
        kvm_close( kd );
        return -1;
    }
    kvm_close( kd );

    /* VmSize */
    val[0] = B2K( (uintmax_t)kp->ki_size );
    /* VmPeak - not defined in /compat/linux/proc/<PID>/status */
    val[1] = -1;
    return 0;
}

#elif defined ASTER_PLATFORM_MINGW
ASTERINTEGER DEFP( MEMPID, mempid, ASTERINTEGER *val ) {
    PROCESS_MEMORY_COUNTERS pmc; // PROCESS_MEMORY_COUNTERS_EX is the same but with one additional
                                 // field the PrivateUsage one.
    GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof( pmc ) );
    /* VmSize */
    val[0] = (ASTERINTEGER)pmc.WorkingSetSize / 1024;
    /* VmPeak */
    val[1] = (ASTERINTEGER)pmc.PeakWorkingSetSize / 1024;
    if ( val[1] == 0 )
        val[1] = -1;
    /* VmRSS */
    // val[3] = (ASTERINTEGER)pmc.PrivateUsage/1024;
    /* VmStk */
    return 0; // stack size on windows ?
}
#endif

#ifdef ASTER_PLATFORM_LINUX
/**
 * Read the memory use by the process from /proc/<PID>/status (default).
 * Returns -1 in case of error, 0 if it succeeds.
 */
ASTERINTEGER get_memory_usage_proc( pid_t pid, ASTERINTEGER *val ) {
    char status_path[256];
    FILE *status_file;
    long memory_kb = -1;
    char line[256];

    snprintf( status_path, sizeof( status_path ), "/proc/%d/status", pid );
    status_file = fopen( status_path, "r" );
    if ( !status_file ) {
        perror( "ERROR: fopen /proc/<PID>/status" );
        return -1;
    }

    while ( fgets( line, sizeof( line ), status_file ) ) {
        // VmRSS or VmSize
        if ( sscanf( line, "VmSize: %ld kB", &memory_kb ) == 1 ) {
            val[0] = memory_kb;
        }
        if ( sscanf( line, "VmPeak: %ld kB", &memory_kb ) == 1 ) {
            val[1] = memory_kb;
        }
    }
    fclose( status_file );
    return 0;
}

/**
 * Read the memory use by the process from cgroup v2.
 * Returns -1 in case of error, 0 if it succeeds.
 */
ASTERINTEGER get_memory_usage_cgroup_v2( pid_t pid, ASTERINTEGER *val ) {
    char cgroup_path[256] = { 0 };
    char memory_path[512] = { 0 };
    FILE *cgroup_file = NULL;
    FILE *memory_file = NULL;
    long memory_bytes;

    snprintf( cgroup_path, sizeof( cgroup_path ), "/proc/%d/cgroup", pid );
    cgroup_file = fopen( cgroup_path, "r" );
    if ( !cgroup_file ) {
        perror( "ERROR: fopen /proc/<PID>/cgroup" );
        return -1;
    }

    // Lire la première ligne (format: 0::/<chemin_du_cgroup>)
    if ( !fgets( cgroup_path, sizeof( cgroup_path ), cgroup_file ) ) {
        perror( "ERROR: read first line of /proc/<PID>/cgroup, fgets" );
        fclose( cgroup_file );
        return -1;
    }
    fclose( cgroup_file );

    cgroup_path[strcspn( cgroup_path, "\n" )] = '\0';
    // Extraire le chemin du cgroup (après "0::/")
    char *cgroup_rel_path = strchr( cgroup_path, '/' );
    if ( !cgroup_rel_path ) {
        fprintf( stderr, "Format inattendu dans /proc/%d/cgroup\n", pid );
        return -1;
    }

    // memory.current
    snprintf( memory_path, sizeof( memory_path ), "/sys/fs/cgroup%s/memory.current",
              cgroup_rel_path );
    memory_file = fopen( memory_path, "r" );
    if ( !memory_file ) {
        perror( "ERROR: fopen memory.current" );
        return -1;
    }
    if ( fscanf( memory_file, "%ld", &memory_bytes ) != 1 ) {
        perror( "ERROR: fscanf memory.current" );
        fclose( memory_file );
        return -1;
    }
    fclose( memory_file );
    val[0] = memory_bytes / 1024;

    // memory.peak
    snprintf( memory_path, sizeof( memory_path ), "/sys/fs/cgroup%s/memory.peak", cgroup_rel_path );
    memory_file = fopen( memory_path, "r" );
    if ( !memory_file ) {
        perror( "ERROR: fopen memory.current" );
        return -1;
    }
    if ( fscanf( memory_file, "%ld", &memory_bytes ) != 1 ) {
        perror( "ERROR: fscanf memory.current" );
        fclose( memory_file );
        return -1;
    }
    fclose( memory_file );
    val[1] = memory_bytes / 1024;

    return 0;
}
#endif
