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

subroutine op9999(options)
    use parameters_module, only: ST_OK
    use allocate_module
    implicit none
    integer(kind=8), intent(in) :: options
#include "asterc/chkmsg.h"
#include "asterc/dllcls.h"
#include "asterc/lcdiscard.h"
#include "asterc/rmfile.h"
#include "asterfort/apetsc.h"
#include "asterfort/asmpi_checkalarm.h"
#include "asterfort/get_jvbasename.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jefini.h"
#include "asterfort/jelibf.h"
#include "asterfort/jemarq.h"
#include "asterfort/jetass.h"
#include "asterfort/jxcopy.h"
#include "asterfort/jxveri.h"
#include "asterfort/ststat.h"
#include "asterfort/uimpba.h"
#include "asterfort/utmess.h"
#include "jeveux.h"

!   Warning: 'options' has not necessarly the same value on all processes
!   but OnlyProc0 must be set everywhere with the same value.
!   Options:
    integer(kind=8), parameter :: SaveBase = 1, Repack = 4, OnlyProc0 = 8, InfoBase = 16
!   InfoResu = 2, Set = 32 not used here
!   - InfoBase:
!       If enabled, list the objects existing in the database
!       (automatically enabled if SaveBase).
!   - SaveBase:
!       If enabled, the objects must be saved properly.
!       Otherwise, the objects can be wiped out (== automatically called at exit).
!   - Repack:
!       Enabled if RETASSAGE="OUI"
!   - OnlyProc0:
!       The objects are only saved on rank #0 if SaveBase.
!   Same values are in 'fin.py'

    character(len=512) :: path
    integer(kind=8) :: iunres, iunmes
    integer(kind=8) :: idx, iret, nbext
    aster_logical :: info_base, close_base

    call jemarq()

    info_base = iand(options, InfoBase) .ne. 0
    close_base = iand(options, SaveBase) .ne. 0

    call ststat(ST_OK)

!   Cleaning in libraries, warnings, errors, mpi...

#ifdef ASTER_HAVE_PETSC
!   Finalize PETSc
    call apetsc('FIN', ' ', ' ', [0.d0], ' ', 0, 0, iret)
#endif
!

!   Free dynamically loaded components
    call dllcls()

!    call lcdiscard(" ")

    if (iand(options, OnlyProc0) .eq. 0) then
!       Check warning messages in parallel
        call asmpi_checkalarm()
    end if

!   Check error messages of type 'E' not followed by 'F' message
    call chkmsg(1, iret)

    if (info_base) then
!       Remove temporay objects from macro-commands
        call jedetc('G', '.', 1)

!       Print the size of objects existing on the GLOBALE database
        iunmes = iunifi('MESSAGE')
        call uimpba('G', iunmes)
    end if

    if (close_base) then
!       Repacking of the GLOBALE database
        if (iand(options, Repack) .ne. 0) then
            call jetass('G')
        end if

!       Free as_allocate object
        call free_slvec(slvec)

!       Call jxveri to check that the execution is ending properly
        call jxveri()
        call jelibf('SAUVE', 'G', 1)
        call jelibf('DETRUIT', 'V', 1)

!       Effective repacking
        if (iand(options, Repack) .ne. 0) then
            call jxcopy('G', 'GLOBALE', 'V', 'VOLATILE', nbext)
            iunres = iunifi('RESULTAT')
            if (iunres .gt. 0) then
                write (iunres, '(A,I2,A)') &
                    ' <I> <FIN> RETASSAGE DE LA BASE "GLOBALE" EFFECTUEE, ', &
                    nbext, ' FICHIER(S) UTILISE(S).'
            end if
        end if
    end if

    call jedema()

!   The diagnosis of the execution is OK thanks to this message.
!   This message is not printed when FIN is automatically called on termination.
    if (options .ne. 0) then
        call utmess('I', 'SUPERVIS2_99')
    end if
    call jefini('NORMAL', close_base)

    if (.not. close_base) then
        idx = 1
        iret = 0
        do while (iret .eq. 0 .and. idx .lt. 99)
            call get_jvbasename("glob", idx, path)
            call rmfile(path, 1, iret)
            idx = idx+1
        end do
        idx = 1
        iret = 0
        do while (iret .eq. 0 .and. idx .lt. 99)
            call get_jvbasename("vola", idx, path)
            call rmfile(path, 1, iret)
            idx = idx+1
        end do
    end if

end subroutine
