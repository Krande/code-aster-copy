! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

! This module proposes logging-like routines.
! It uses some global variables some performance reasons to keep some indicators
! in cache.

module logging_module
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"

    integer(kind=8), parameter :: LOGLEVEL_SIZE = 3
    ! generic logger
    integer(kind=8), parameter :: LOGLEVEL = 1
    ! for memory (jeveux) traces
    integer(kind=8), parameter :: LOGLEVEL_MEM = 2
    ! for MGIS usages
    integer(kind=8), parameter :: LOGLEVEL_MGIS = 3

    ! levels (same value as in the *logging* Python module + verbose)
    integer(kind=8), parameter :: DEBUG = 10, VERBOSE = 15, INFO = 20
    integer(kind=8), parameter :: WARNING = 30, ERROR = 40, UNSET = -1

    ! current level for each logger
    integer(kind=8) :: level_(LOGLEVEL_SIZE)
    ! environment variables to set each logger
    character(len=24) :: envvar_(LOGLEVEL_SIZE)

    private :: envvar_, level_

contains

!>  Initialization of loggers
    subroutine initialize()
        level_ = UNSET
        envvar_ = [ &
                  "ASTER_LOGLEVEL          ", &
                  "ASTER_LOGLEVEL_MEM      ", &
                  "ASTER_LOGLEVEL_MGIS     " &
                  ]
        call setLevel(LOGLEVEL, INFO)
    end subroutine initialize

!>  Set the level of a logger
    subroutine setLevel(logid, level)
        integer(kind=8), intent(in) :: logid, level
        ASSERT(logid .gt. 0 .and. logid .le. LOGLEVEL_SIZE)
        level_(logid) = level
    end subroutine setLevel

!>  Returns the current level of a logger
    function getLevel(logid)
        integer(kind=8) :: getLevel
        integer(kind=8), intent(in) :: logid

        character(len=8) :: value
        integer(kind=8) :: level
        ASSERT(logid .gt. 0 .and. logid .le. LOGLEVEL_SIZE)
        if (level_(logid) .eq. UNSET) then
            call get_environment_variable(envvar_(logid), value)
            if (value .eq. "DEBUG") then
                level = DEBUG
            elseif (value .eq. "VERBOSE") then
                level = VERBOSE
            elseif (value .eq. "INFO") then
                level = INFO
            elseif (value .eq. "WARNING") then
                level = WARNING
            elseif (value .eq. "ERROR") then
                level = ERROR
            elseif (value .eq. " ") then
                level = INFO
            else
                print *, "WARNING: expecting one of DEBUG, VERBOSE, INFO, WARNING, "&
                   &"ERROR, not: ", value
                level = ERROR
            end if
            call setLevel(logid, level)
        end if
        getLevel = level_(logid)
    end function getLevel

!>  Tell if a logger is enabled at a given level
    function is_enabled(logid, level)
        aster_logical :: is_enabled
        integer(kind=8), intent(in) :: logid, level
        ASSERT(logid .gt. 0 .and. logid .le. LOGLEVEL_SIZE)
        is_enabled = (getLevel(logid) .le. level)
    end function is_enabled

end module
