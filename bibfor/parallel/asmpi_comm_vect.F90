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

subroutine asmpi_comm_vect(optmpi, typsca, nbval, bcrank, vi, &
                           vi4, vr, vc, sci, sci4, &
                           scr, scc)
! person_in_charge: jacques.pellet at edf.fr
!
!
!  FONCTION REALISEE : SUR-COUCHE MPI
!
!  FAIRE UN ECHANGE BCAST/REDUCE/ALL_REDUCE SUR UN VECTEUR FORTRAN
!
! Arguments d'appels
! in optmpi :
!       /'MPI_MAX'  == 'ALLREDUCE + MAX' (interdit pour typsca= 'C')
!       /'MPI_MIN'  == 'ALLREDUCE + MIN' (interdit pour typsca= 'C')
!       /'MPI_SUM'  == 'ALLREDUCE + SUM'
!
!       /'REDUCE'   == 'REDUCE + SUM'      : tous -> 0
!       /'BCAST'    == 'BCAST'             : proc de rang=bcrank -> tous
!
! in    typsca : /'I' /'S' /'R' /'C'
! in    nbval  : longueur du vecteur v* (optionnel, 1 si absent)
! in    bcrank : rang du processus mpi d'ou emane le bcast
!-si nbval > 1:
! inout vi(*)  : vecteur d'entiers a echanger    (si typsca='I')
! inout vi4(*) : vecteur d'entiers a echanger    (si typsca='S')
! inout vr(*)  : vecteur de reels a echanger     (si typsca='R')
! inout vc(*)  : vecteur de complexes a echanger (si typsca='C')
!-si nbval == 1:
! inout sci    : entier a echanger    (si typsca='I')
! inout sci4   : entier a echanger    (si typsca='S')
! inout scr    : reel a echanger      (si typsca='R')
! inout scc    : complexe a echanger  (si typsca='C')
!----------------------------------------------------------------------
    implicit none
!
#include "asterc/asmpi_comm.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_check.h"
#include "asterfort/asmpi_comm_mvect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpu.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: optmpi
    character(len=*), intent(in) :: typsca
    integer(kind=8), intent(in), optional :: nbval
    integer(kind=8), intent(in), optional :: bcrank
    integer(kind=8), intent(inout), optional :: vi(*)
    integer(kind=4), intent(inout), optional :: vi4(*)
    real(kind=8), intent(inout), optional :: vr(*)
    complex(kind=8), intent(inout), optional :: vc(*)
    integer(kind=8), intent(inout), optional :: sci
    integer(kind=4), intent(inout), optional :: sci4
    real(kind=8), intent(inout), optional :: scr
    complex(kind=8), intent(inout), optional :: scc
!
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
#include "asterf_mpi.h"
!
    character(len=1) :: typsc1
    integer(kind=8) :: iret
    mpi_int :: nbpro4, mpicou, proc
! ---------------------------------------------------------------------
    call jemarq()
!
!   -- communicateur mpi de travail :
    call asmpi_comm('GET', mpicou)
!   -- compteur de temps CPU :
    call uttcpu('CPU.CMPI.1', 'DEBUT', ' ')
!
!   -- s'il n'y a qu'un seul proc, il n'y a rien a faire :
    call asmpi_info(mpicou, rank=proc, size=nbpro4)
    if (nbpro4 .eq. 1) goto 999
    DEBUG_MPI('mpi_comm_vect', proc, nbpro4)
!
!
!   -- verification rendez-vous - a voir si on garde car beacuoup de comm
    iret = 1
#ifdef ASTER_DEBUG_MPI
    call asmpi_check(iret)
    if (iret .ne. 0) then
        call utmess('I', 'APPELMPI_83', sk=optmpi)
        goto 999
    end if
#endif
!
!
!   -- Calcul de typsc1 :
!   ---------------------------------------
    typsc1 = typsca
    ASSERT(typsc1 .eq. 'R' .or. typsc1 .eq. 'C' .or. typsc1 .eq. 'I' .or. typsc1 .eq. 'S')
    ASSERT(typsc1 .ne. 'R' .or. present(vr) .or. present(scr))
    ASSERT(typsc1 .ne. 'C' .or. present(vc) .or. present(scc))
    ASSERT(typsc1 .ne. 'I' .or. present(vi) .or. present(sci))
    ASSERT(typsc1 .ne. 'S' .or. present(vi4) .or. present(sci4))
!
!
!   -- communication :
!   ----------------------------
!
    if (present(sci)) then
        call asmpi_comm_mvect(optmpi, typsca, 1, bcrank, sci=sci)
    else if (present(sci4)) then
        call asmpi_comm_mvect(optmpi, typsca, 1, bcrank, sci4=sci4)
    else if (present(scr)) then
        call asmpi_comm_mvect(optmpi, typsca, 1, bcrank, scr=scr)
    else if (present(scc)) then
        call asmpi_comm_mvect(optmpi, typsca, 1, bcrank, scc=scc)
    else if (present(vi)) then
        call asmpi_comm_mvect(optmpi, typsca, nbval, bcrank, vi=vi)
    else if (present(vi4)) then
        call asmpi_comm_mvect(optmpi, typsca, nbval, bcrank, vi4=vi4)
    else if (present(vr)) then
        call asmpi_comm_mvect(optmpi, typsca, nbval, bcrank, vr=vr)
    else if (present(vc)) then
        call asmpi_comm_mvect(optmpi, typsca, nbval, bcrank, vc=vc)
    else
        ASSERT(ASTER_FALSE)
    end if
!
999 continue
    call uttcpu('CPU.CMPI.1', 'FIN', ' ')
    call jedema()
!
#else
    character(len=1) :: kdummy
    integer(kind=8) :: idummy
    integer(kind=4) :: i4dummy
    real(kind=8) :: rdummy
    complex(kind=8) :: cdummy
!
    if (present(nbval) .and. present(vi) .and. present(vr) .and. present(vc) .and. &
        present(bcrank) .and. present(sci) .and. present(scr) .and. present(scc)) then
        kdummy = optmpi(1:1)
        kdummy = typsca(1:1)
        idummy = nbval
        idummy = bcrank
        idummy = vi(1)
        i4dummy = vi4(1)
        rdummy = vr(1)
        cdummy = vc(1)
        idummy = sci
        i4dummy = sci4
        rdummy = scr
        cdummy = scc
    end if
#endif
end subroutine
