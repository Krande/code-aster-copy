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
#include "asterc/jdcset.h"
#include "asterc/rmfile.h"
#include "asterfort/assert.h"
#include "asterfort/fin999.h"
#include "asterfort/get_jvbasename.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jefini.h"
#include "asterfort/jeimhd.h"
#include "asterfort/jeliad.h"
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
    integer, parameter :: SaveBase = 1, FormatHdf = 4, Repack = 8
!   InfoResu = 2, OnlyProc0 = 16: not used here
!   - SaveBase:
!       If enabled, the objects must be saved properly.
!       Otherwise, the objects can be wiped out (== automatically called).
!   - FormatHdf:
!       Enabled if FormatHdf="OUI"
!   - Repack:
!       Enabled if RETASSAGE="OUI"
!   Same values are in 'fin.py'

    integer :: nbenre, nboct, iret
    integer :: iunres, iunmes
    integer :: i, nbext
    aster_logical :: close_base
    character(len=80) :: fich
    character(len=256) :: fbase

    call jemarq()

    close_base = iand(options, SaveBase) .ne. 0

    call ststat(ST_OK)

!   Cleaning in libraries, warnings, errors, mpi...
    call fin999()

    if ( close_base ) then
!       Remove temporay objects from macro-commands
        call jedetc('G', '.', 1)

!       Print the size of objects existing on the GLOBALE database
        iunmes = iunifi('MESSAGE')
        call uimpba('G', iunmes)

!       Repacking of the GLOBALE database
        if ( iand(options, Repack) .ne. 0 ) then
            call jetass('G')
            if ( iand(options, FormatHdf) .ne. 0 ) then
                call utmess('A', 'SUPERVIS2_8')
            endif
        endif

!       Save the GLOBALE database in HDF5 format
        if ( iand(options, FormatHdf) .ne. 0 ) then
            fich = 'bhdf.1'
            call jeimhd(fich, 'G')
        endif

    endif

!   Get the location of a specific record to identify the execution
    call jeliad('G', nbenre, nboct)
    call jdcset('jeveux_sysaddr', nboct)

!   Call jxveri to check that the execution is ending properly
    if ( close_base ) then
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
