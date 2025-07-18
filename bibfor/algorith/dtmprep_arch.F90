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

subroutine dtmprep_arch(sd_dtm_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmprep_arch : Prepare the archiving list prior to calculation, parse the
!                user input under the keywords ARCHIVAGE=_F(...
!
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/getfac.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
!
!   -0.2- Local variables
    aster_logical     :: checktfin
    integer(kind=8)           :: iret1, iret2, iparch, nbarch, nbsaves
    integer(kind=8)           :: nbinst, nbocc, i, j, sizearch
    real(kind=8)      :: tinit, tfin, epsi, dt, residue, perarch
    character(len=8)  :: sd_dtm
    character(len=19) :: numarc
!
    real(kind=8), target  :: list(1)
    real(kind=8), pointer :: inst(:) => null()
    real(kind=8), pointer :: archlst(:) => null()
!
!   0 - Initializations
    call jemarq()
    sd_dtm = sd_dtm_
    epsi = r8prem()

    call dtmget(sd_dtm, _INST_FIN, rscal=tfin)

    list = [tfin]
    nbinst = 0
    inst => list
    iparch = 1
    iret2 = 0
!
    call getfac('ARCHIVAGE', nbocc)
    if (nbocc .ne. 0) then
        call getvid('ARCHIVAGE', 'LIST_INST', iocc=1, scal=numarc, nbret=iret1)
        if (iret1 .ne. 0) then
            nullify (inst)
            call jeveuo(numarc//'.VALE', 'L', vr=inst)
            nbinst = size(inst)
        else
            call getvr8('ARCHIVAGE', 'INST', iocc=1, nbval=0, nbret=iret2)
            if (iret2 .ne. 0) then
                nullify (inst)
                nbinst = -iret2
                AS_ALLOCATE(vr=inst, size=nbinst)
                call getvr8('ARCHIVAGE', 'INST', iocc=1, nbval=nbinst, vect=inst)
            end if
        end if
        call getvis('ARCHIVAGE', 'PAS_ARCH', iocc=1, scal=iparch, nbret=iret1)
        if (iret1 .eq. 0) iparch = ismaem()-1
        if (iparch .lt. 0) iparch = ismaem()
    end if
!
    call dtmget(sd_dtm, _INST_INI, rscal=tinit)
    call dtmget(sd_dtm, _DT, rscal=dt)
!
!   --- Check that the list of instants belong to the interval {tinit - tfin}
    sizearch = nbinst+2
    AS_ALLOCATE(vr=archlst, size=sizearch)
    nbarch = 0

    perarch = real(iparch)*dt
    if (perarch .gt. (tfin-tinit)) then
        nbsaves = 1
    else
        nbsaves = 1+nint((tfin-tinit)/perarch)
        nbsaves = 1+int((tfin+epsi*nbsaves-tinit)/perarch)
    end if

    do i = 1, nbinst
        if (inst(i) .lt. (tinit-epsi)) then
            ASSERT(.false.)
        else if (inst(i) .gt. (tfin+epsi)) then
            ASSERT(.false.)
        else
            if (abs(inst(i)-tinit) .gt. epsi) then
                nbarch = nbarch+1
                archlst(nbarch) = inst(i)
                j = int((inst(i)-tinit)/dt)
                residue = mod(inst(i)-tinit, perarch)
                residue = min(residue, abs(inst(i)-perarch))
                if ((residue .gt. (j*epsi)) .and. (abs(residue-perarch) .gt. (j*epsi))) then
                    nbsaves = nbsaves+1
                end if
            end if
        end if
    end do

    checktfin = .false.
    if (nbarch .gt. 0) then
        if (abs(archlst(nbarch)-tfin) .gt. epsi) checktfin = .true.
    end if

    if ((nbarch .eq. 0) .or. (checktfin)) then
        nbarch = nbarch+1
        archlst(nbarch) = tfin
        j = int((tfin-tinit)/dt)
        residue = mod(tfin-tinit, perarch)
        residue = min(residue, abs(tfin-perarch))
        if ((residue .gt. (j*epsi)) .and. (abs(residue-perarch) .gt. (j*epsi))) then
            nbsaves = nbsaves+1
        end if
    end if

!   --- The archived instant list must be monotone, check renforced
    do i = 1, nbarch-1
        if (archlst(i) .ge. archlst(i+1)) then
            ASSERT(.false.)
        end if
    end do
!
    if (iret2 .ne. 0) then
        AS_DEALLOCATE(vr=inst)
    end if
!
    call dtmsav(sd_dtm, _AR_LINST, nbarch, rvect=archlst)
    call dtmsav(sd_dtm, _ARCH_NB, 1, iscal=nbsaves)
    call dtmsav(sd_dtm, _ARCH_PER, 1, iscal=iparch)
!
    AS_DEALLOCATE(vr=archlst)

    call jedema()
end subroutine
