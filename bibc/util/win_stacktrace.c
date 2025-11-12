/* --------------------------------------------------------------------
// Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
// This file is part of code_aster.
//
// code_aster is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// code_aster is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
// --------------------------------------------------------------------
*/

#ifdef ASTER_PLATFORM_MSVC64

#include <windows.h>
#include <dbghelp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma comment(lib, "dbghelp.lib")

#define MAX_STACK_FRAMES 64
#define MAX_SYMBOL_NAME 512

static int stacktrace_initialized = 0;

// Print stack trace to file
void win_print_stacktrace(FILE* fp) {
    if (!fp) fp = stderr;

    HANDLE process = GetCurrentProcess();
    HANDLE thread = GetCurrentThread();

    void* stack[MAX_STACK_FRAMES];
    WORD frames = CaptureStackBackTrace(0, MAX_STACK_FRAMES, stack, NULL);

    fprintf(fp, "\nC/C++ Stack trace:\n");
    fprintf(fp, "==================\n");

    // Try to initialize if not already done (might already be done by Intel Fortran)
    // We don't care if this fails - symbols might already be initialized
    if (!stacktrace_initialized) {
        SymSetOptions(SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS | SYMOPT_LOAD_LINES);

        // Note: With /Z7 debug format, debug info is embedded in the executable/DLL
        // Source files are only needed if user wants to see source code in debugger
        // Symbol resolution (function names, line numbers) works without source files

        if (SymInitialize(process, NULL, TRUE) || GetLastError() == ERROR_INVALID_PARAMETER) {
            stacktrace_initialized = 1;
        }
    }

    SYMBOL_INFO* symbol = (SYMBOL_INFO*)calloc(sizeof(SYMBOL_INFO) + MAX_SYMBOL_NAME, 1);
    if (symbol) {
        symbol->MaxNameLen = MAX_SYMBOL_NAME - 1;
        symbol->SizeOfStruct = sizeof(SYMBOL_INFO);

        for (WORD i = 0; i < frames; i++) {
            DWORD64 address = (DWORD64)(stack[i]);

            // Get symbol name
            if (SymFromAddr(process, address, 0, symbol)) {
                // Get line number information
                IMAGEHLP_LINE64 line = { 0 };
                line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);
                DWORD displacement = 0;

                if (SymGetLineFromAddr64(process, address, &displacement, &line)) {
                    fprintf(fp, "[%2d] %-40s  %s:%lu\n",
                            i, symbol->Name, line.FileName, line.LineNumber);
                } else {
                    fprintf(fp, "[%2d] %-40s  (no line info)\n", i, symbol->Name);
                }
            } else {
                // If symbol not found, at least print the address
                fprintf(fp, "[%2d] 0x%016llX  (no symbol)\n", i, address);
            }
        }

        free(symbol);
    }

    fprintf(fp, "\n");
    fflush(fp);
}

// Vectored exception handler for crashes (higher priority than Python's handler)
static LONG WINAPI vectored_exception_handler(EXCEPTION_POINTERS* exceptionInfo) {
    DWORD exceptionCode = exceptionInfo->ExceptionRecord->ExceptionCode;

    // Only handle access violations and other serious exceptions
    if (exceptionCode == EXCEPTION_ACCESS_VIOLATION ||
        exceptionCode == EXCEPTION_ARRAY_BOUNDS_EXCEEDED ||
        exceptionCode == EXCEPTION_DATATYPE_MISALIGNMENT ||
        exceptionCode == EXCEPTION_FLT_DENORMAL_OPERAND ||
        exceptionCode == EXCEPTION_FLT_DIVIDE_BY_ZERO ||
        exceptionCode == EXCEPTION_FLT_INEXACT_RESULT ||
        exceptionCode == EXCEPTION_FLT_INVALID_OPERATION ||
        exceptionCode == EXCEPTION_FLT_OVERFLOW ||
        exceptionCode == EXCEPTION_FLT_STACK_CHECK ||
        exceptionCode == EXCEPTION_FLT_UNDERFLOW ||
        exceptionCode == EXCEPTION_ILLEGAL_INSTRUCTION ||
        exceptionCode == EXCEPTION_IN_PAGE_ERROR ||
        exceptionCode == EXCEPTION_INT_DIVIDE_BY_ZERO ||
        exceptionCode == EXCEPTION_INT_OVERFLOW ||
        exceptionCode == EXCEPTION_INVALID_DISPOSITION ||
        exceptionCode == EXCEPTION_NONCONTINUABLE_EXCEPTION ||
        exceptionCode == EXCEPTION_PRIV_INSTRUCTION ||
        exceptionCode == EXCEPTION_STACK_OVERFLOW) {

        fprintf(stderr, "\n");
        fprintf(stderr, "====================================\n");
        fprintf(stderr, "Windows Exception Caught!\n");
        fprintf(stderr, "Exception Code: 0x%08lX\n", exceptionCode);
        fprintf(stderr, "Exception Address: 0x%p\n", exceptionInfo->ExceptionRecord->ExceptionAddress);
        fprintf(stderr, "====================================\n");
        fflush(stderr);

        // Print our detailed C/C++ stack trace with line numbers
        // This works for all code (C/C++/Fortran) and shows file:line for C/C++
        win_print_stacktrace(stderr);

        // Note: Intel Fortran's traceback will also be printed by its runtime handler
        // That traceback shows better info for Fortran code, but "Unknown" for C/C++
        // Together, these two tracebacks give complete information:
        // - Our trace: Full info for C/C++, function names only for Fortran
        // - Intel trace: Full info for Fortran, "Unknown" for C/C++

        fprintf(stderr, "\nNote: Intel Fortran traceback will follow (if Fortran code involved)\n");
        fprintf(stderr, "      It shows Fortran source but C/C++ shows as 'Unknown'\n");
        fprintf(stderr, "      The above C/C++ trace shows the missing C/C++ details.\n\n");
        fflush(stderr);

        // Don't handle it, let it continue to other handlers (including Python's and Intel's)
        // This way we get our C/C++ stack trace, Intel's Fortran trace, and Python's exception
    }

    return EXCEPTION_CONTINUE_SEARCH;
}

static PVOID vectored_handler = NULL;

// Install exception handler
void win_install_crash_handler() {
    // Don't initialize symbols here - let Intel Fortran do it if needed
    // We'll initialize on-demand when we actually need to print a stack trace

    // Use vectored exception handler which has higher priority than structured exception handling
    if (!vectored_handler) {
        vectored_handler = AddVectoredExceptionHandler(1, vectored_exception_handler);
    }
}


// Cleanup
void win_cleanup_stacktrace() {
    if (vectored_handler) {
        RemoveVectoredExceptionHandler(vectored_handler);
        vectored_handler = NULL;
    }
    // Don't call SymCleanup here - let the process cleanup handle it
    // This prevents conflicts with Intel Fortran runtime which may also use DbgHelp
    if (stacktrace_initialized) {
        stacktrace_initialized = 0;
    }
}

#else

// Stub functions for non-Windows platforms
#include <stdio.h>

void win_print_stacktrace(FILE* fp) {
    (void)fp;
}

void win_install_crash_handler() {
}

int win_is_crash_handler_installed() {
    return 0;
}

void win_cleanup_stacktrace() {
}

#endif

