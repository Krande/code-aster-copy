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
! person_in_charge: nicolas.sellenet at edf.fr
!
subroutine impres_component_hpc(nomgd, ntncmp, ncmpvl, ncmpve, indcmp)
!
    implicit none
    character(len=32)  :: nomgd
    character(len=*)  :: ntncmp
    integer(kind=8) :: ncmpvl, ncmpve
    character(len=24)  :: indcmp
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/wkvect.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_allgatherv_char16.h"
!
! --------------------------------------------------------------------------------------------------
!
!     IMPR_RESU components communication
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jnocm1, cmpt, icmp1, vnbcmp(1), nbcmpmax, rang, nbproc
    integer(kind=8) :: jtest, jnocmp, jnocm3, nb_cmp_tot, iproc, jindir, numcmp
    mpi_int, pointer :: v_count(:) => null()
    mpi_int, pointer :: v_displ(:) => null()
    character(len=16), pointer :: v_nomcmp(:) => null()
    character(len=16), pointer :: v_nomcm2(:) => null()
    mpi_int :: mrank, mnbproc, world, one4, taille
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    cmpt = 0
    if (ncmpve .ne. 0) then
        call jeveuo(ntncmp, 'L', jnocm1)
        call jecreo('&&IRCAME.CMPLOC', 'V N K16')
        call jeecra('&&IRCAME.CMPLOC', 'NOMMAX', ncmpve)
        cmpt = 0
        do icmp1 = 1, ncmpve
            call jecroc(jexnom('&&IRCAME.CMPLOC', zk16(jnocm1+icmp1-1)))
        end do
    end if
    one4 = to_mpi_int(1)
    if (nomgd .eq. 'VARI_R') then
        vnbcmp(1) = ncmpve
        call asmpi_comm_vect('MPI_MAX', 'I', 1, vi=vnbcmp)
        nbcmpmax = vnbcmp(1)
    else
        call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmpmax)
    end if
    call asmpi_comm('GET', world)
    call asmpi_info(rank=mrank, size=mnbproc)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(mnbproc)
    call wkvect('&&IRCAME.TEST', 'V V I', nbproc, jtest)
    call wkvect('&&IRCAME.COUNT', 'V V S', nbproc, vi4=v_count)
    call wkvect('&&IRCAME.DISPL', 'V V S', nbproc+1, vi4=v_displ)
    if (ncmpve .ne. 0) then
        call wkvect('&&IRCAME.NOMCMP2', 'V V K16', ncmpve, jnocmp)
        call wkvect('&&IRCAME.NOMCMP3', 'V V K16', ncmpve, jnocm3)
    else
        call wkvect('&&IRCAME.NOMCMP2', 'V V K16', 1, jnocmp)
        call wkvect('&&IRCAME.NOMCMP3', 'V V K16', 1, jnocm3)
    end if
    call asmpi_allgather_i([ncmpve], one4, zi(jtest), one4, world)
    nb_cmp_tot = 0
    do iproc = 0, nbproc-1
        nb_cmp_tot = nb_cmp_tot+zi(jtest+iproc)
        v_count(iproc+1) = to_mpi_int(zi(jtest+iproc))
        v_displ(iproc+2) = to_mpi_int(nb_cmp_tot)
    end do
    do icmp1 = 1, ncmpve
        zk16(jnocmp+icmp1-1) = zk16(jnocm1+icmp1-1)
    end do
    if (nb_cmp_tot .ne. 0) then
        call wkvect('&&IRCAME.NOMFAG', 'V V K16', nb_cmp_tot, vk16=v_nomcmp)
        call wkvect('&&IRCAME.NOMFA2', 'V V K16', nb_cmp_tot, vk16=v_nomcm2)
        call jedetr(ntncmp)
        call wkvect(ntncmp, 'V V K16', nb_cmp_tot, jnocm1)
        call wkvect(indcmp, 'V V I', nb_cmp_tot, jindir)
        taille = to_mpi_int(ncmpve)
        call asmpi_allgatherv_char16(zk16(jnocmp), taille, v_nomcmp, v_count, v_displ, world)
        call asmpi_allgatherv_char16(zk16(jnocm3), taille, v_nomcm2, v_count, v_displ, world)
        call jecreo('&&IRCAME.PTRNOM', 'V N K16')
        call jeecra('&&IRCAME.PTRNOM', 'NOMMAX', nb_cmp_tot)
        cmpt = 0
        do icmp1 = 1, nb_cmp_tot
            call jenonu(jexnom('&&IRCAME.PTRNOM', v_nomcmp(icmp1)), numcmp)
            if (numcmp .eq. 0) then
                call jecroc(jexnom('&&IRCAME.PTRNOM', v_nomcmp(icmp1)))
                cmpt = cmpt+1
                zk16(jnocm1+cmpt-1) = v_nomcmp(icmp1)
                if (ncmpve .ne. 0) then
                    call jenonu(jexnom('&&IRCAME.CMPLOC', v_nomcmp(icmp1)), numcmp)
                    if (numcmp .ne. 0) then
                        zi(jindir+icmp1-1) = cmpt
                    end if
                end if
            end if
        end do
    end if
    call jedetr('&&IRCAME.TEST')
    call jedetr('&&IRCAME.COUNT')
    call jedetr('&&IRCAME.DISPL')
    call jedetr('&&IRCAME.NOMCMP2')
    call jedetr('&&IRCAME.NOMCMP3')
    call jedetr('&&IRCAME.NOMFAG')
    call jedetr('&&IRCAME.NOMFA2')
    call jedetr('&&IRCAME.PTRNOM')
    call jedetr('&&IRCAME.CMPLOC')
    ncmpvl = ncmpve
    ncmpve = cmpt
!
    call jedema()
!
end subroutine
