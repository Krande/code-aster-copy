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

#ifndef WIN_STACKTRACE_H
#define WIN_STACKTRACE_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Print stack trace to file (defaults to stderr if NULL)
void win_print_stacktrace(FILE* fp);

// Install Windows exception handler for crash reporting
void win_install_crash_handler();


// Cleanup symbol handler
void win_cleanup_stacktrace();

#ifdef __cplusplus
}
#endif

#endif // WIN_STACKTRACE_H

