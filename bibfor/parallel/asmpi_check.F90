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

subroutine asmpi_check(iret)
! person_in_charge: mathieu.courtois at edf.fr
!
!
    use parameters_module
    implicit none
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_wtime.h"
#include "asterc/uttrst.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/asmpi_status.h"
#include "asterfort/gtstat.h"
#include "asterfort/asmpi_stop.h"
#include "asterfort/ststat.h"
#include "asterfort/utmess.h"
    integer(kind=8), intent(out) :: iret
!-----------------------------------------------------------------------
!     FONCTION REALISEE : MPI CHECK ERROR
!       AVANT D'EFFECTUER UNE COMMUNICATION BLOQUANTE
!       ON VERIFIE :
!           - QU'AUCUN PROCESSEUR N'A SIGNALE DE PROBLEMES
!           - QUE TOUS LES PROCESSEURS SONT AU RENDEZ-VOUS
!       EN CAS DE PROBLEME, ON RETOURNE IRET != 0, CAR IL NE FAUT ALORS
!       PAS INITIER DE NOUVELLES COMMUNICATIONS.
!-----------------------------------------------------------------------
#if defined(ASTER_HAVE_MPI) && !defined(ASTER_DISABLE_MPI_CHECK)

#include "mpif.h"
#include "asterc/asmpi_irecv_i4.h"
#include "asterc/asmpi_send_i4.h"
#include "asterc/asmpi_cancel.h"
#include "asterc/asmpi_test.h"

    aster_logical :: lcont
    mpi_int :: term
    integer(kind=8) :: i, nbterm, nbproc, np1, resp0
    mpi_int :: nbpro4, rank, istat, mpicou, wki(1), nbv, ip4
    real(kind=8) :: tres, timeout, t0, tf
    aster_logical, allocatable :: isterm(:)
    mpi_int, allocatable :: diag(:)
    mpi_int, allocatable :: request(:)

!   Current mpi communicator
    call asmpi_comm('GET', mpicou)
    iret = 0
    call asmpi_info(mpicou, rank=rank, size=nbpro4)
    nbproc = to_aster_int(nbpro4)
!   if not started by mpiexec for debugging
    if (nbproc .le. 1) then
        goto 999
    end if
    np1 = nbpro4-1
    nbv = 1

    DEBUG_MPI('mpi_check', rank, nbpro4)

!   On the processor #0
    if (rank == 0) then
        allocate (isterm(nbproc))
        allocate (diag(nbproc))
        allocate (request(nbproc))

        call uttrst(tres)
        timeout = tres*0.2d0
        if (timeout < 0) then
            timeout = 120.
            call utmess('A', 'APPELMPI_94', sr=timeout)
        end if

!       Ask each processor for its status
        do i = 1, np1
            isterm(i) = .false.
            ip4 = i
            DEBUG_MPI('mpi_check', 'irecv from ', ip4)
            call asmpi_irecv_i4(diag(i:i), nbv, ip4, ST_TAG_CHK, mpicou, &
                                request(i))
        end do
!
        nbterm = 0
        t0 = asmpi_wtime()
        lcont = .true.
        do while (lcont)
            do i = 1, np1
                if (.not. isterm(i)) then
                    call asmpi_test(request(i), term)
                    if (term .eq. 1) then
                        nbterm = nbterm+1
                        isterm(i) = .true.
                        if (diag(i) .eq. ST_ER) then
                            call utmess('I', 'APPELMPI_84', si=i)
                            call ststat(ST_ER_OTH)
                        end if
                    end if
                end if
            end do
            lcont = nbterm .lt. np1
!           timeout
            tf = asmpi_wtime()
            if (lcont .and. (tf-t0) .gt. timeout) then
                lcont = .false.
                call utmess('E', 'APPELMPI_97', sr=timeout)
                do i = 1, np1
                    if (.not. isterm(i)) then
                        call utmess('E+', 'APPELMPI_96', si=i)
                        call utmess('E', 'APPELMPI_83', sk='MPI_IRECV')
                        call ststat(ST_UN_OTH)
                    end if
                end do
            end if
        end do

        if (gtstat(ST_ER_PR0)) then
            call utmess('I', 'APPELMPI_84', si=0)
        end if
!       Tell to all processors that answered if it continues or not
        if (gtstat(ST_OK)) then
            istat = ST_OK
        else
            istat = ST_ER
        end if
        do i = 1, np1
            if (isterm(i)) then
                if (istat .ne. ST_OK) then
                    call utmess('I', 'APPELMPI_81', si=i)
                end if
                wki(1) = istat
                ip4 = i
                DEBUG_MPI('mpi_check:send status / to', wki(1), ip4)
                call asmpi_send_i4(wki, nbv, ip4, ST_TAG_CNT, mpicou)
            else
!               cancel those have not answered
                DEBUG_MPI('mpi_check', 'cancel request for proc', i)
                call asmpi_cancel(request(i))
            end if
        end do
!
        if (.not. gtstat(ST_OK)) then
            iret = 1
            if (gtstat(ST_UN_OTH)) then
                call asmpi_stop(1)
            else
                call asmpi_stop(2)
            end if
        end if
        deallocate (isterm)
        deallocate (diag)
        deallocate (request)

!   On all others processors (not #0)
    else
!       Each processor sends ST_OK to the processor #0
        call asmpi_status(ST_OK, resp0)
        if (resp0 .ne. ST_OK) then
            iret = 1
            call utmess('I', 'APPELMPI_80')
            call asmpi_stop(2)
        end if
    end if
999 continue

#else
    iret = 0
#endif

end subroutine asmpi_check
