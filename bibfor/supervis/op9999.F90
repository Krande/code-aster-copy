! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
    use parameters_module, only : ST_OK
    implicit none
    integer, intent(in) :: options
#include "asterc/chkmsg.h"
#include "asterc/dllcls.h"
#include "asterc/lcdiscard.h"
#include "asterc/rmfile.h"
#include "asterfort/apetsc.h"
#include "asterfort/asmpi_checkalarm.h"
#include "asterfort/assert.h"
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
#include "asterfort/rsinfo.h"
#include "asterfort/ststat.h"
#include "asterfort/uimpba.h"
#include "asterfort/ulexis.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"

!   Warning: 'options' has not necessarly the same value on all processes!
!   Options:
    integer, parameter :: SaveBase = 1, Repack = 4
!   InfoResu = 2, OnlyProc0 = 8: not used here
!   - SaveBase:
!       If enabled, the objects must be saved properly.
!       Otherwise, the objects can be wiped out (== automatically called).
!   - Repack:
!       Enabled if RETASSAGE="OUI"
!   Same values are in 'fin.py'

    integer :: iunres, iunmes
    integer :: i, iret, nbext
    aster_logical :: close_base
    character(len=256) :: fbase

    call jemarq()

    close_base = iand(options, SaveBase) .ne. 0

    call ststat(ST_OK)

!   Cleaning in libraries, warnings, errors, mpi...

#ifdef ASTER_HAVE_PETSC
!   Finalize PETSc
    call apetsc('FIN', ' ', ' ', [0.d0], ' ', 0, 0, iret)
#endif

!   Free dynamically loaded components
    call dllcls()

!    call lcdiscard(" ")

    if ( close_base ) then
!       Check warning messages in parallel
        call asmpi_checkalarm()

!       Check error messages of type 'E' not followed by 'F' message
        call chkmsg(1, iret)

!       Remove temporay objects from macro-commands
        call jedetc('G', '.', 1)

!       Print the size of objects existing on the GLOBALE database
        iunmes = iunifi('MESSAGE')
        call uimpba('G', iunmes)

!       Repacking of the GLOBALE database
        if ( iand(options, Repack) .ne. 0 ) then
            call jetass('G')
        endif

!       Call jxveri to check that the execution is ending properly
        call jxveri()
        call jelibf('SAUVE', 'G', 1)
        call jelibf('DETRUIT', 'V', 1)

!       Effective repacking
        if ( iand(options, Repack) .ne. 0 ) then
            call jxcopy('G', 'GLOBALE', 'V', 'VOLATILE', nbext)
            iunres = iunifi('RESULTAT')
            if (iunres .gt. 0) then
                write(iunres,'(A,I2,A)') &
                    ' <I> <FIN> RETASSAGE DE LA BASE "GLOBALE" EFFECTUEE, ',&
                    nbext, ' FICHIER(S) UTILISE(S).'
            endif
        endif
    endif

    call jedema()

!   The diagnosis of the execution is OK thanks to this message
    if ( close_base ) then
        call utmess('I', 'SUPERVIS2_99')
    endif
    call jefini('NORMAL', close_base)

    if ( .not. close_base ) then
        do i=1,99
            call get_jvbasename('glob', i, fbase)
            call rmfile(fbase, 0, iret)
            if (iret .ne. 0) then
                exit
            endif
            call get_jvbasename('vola', i, fbase)
            call rmfile(fbase, 0, iret)
        end do
    endif

end subroutine
