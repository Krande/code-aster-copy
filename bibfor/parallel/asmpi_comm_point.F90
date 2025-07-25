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

subroutine asmpi_comm_point(optmpi, typsca, nudest, numess, nbval, &
                            vi, vi4, vr, sci, sci4, &
                            scr)
! person_in_charge: nicolas.sellenet at edf.fr
!
!
!  FONCTION REALISEE : SUR-COUCHE MPI
!
!  COMMUNICATION MPI POINT A POINT D'UN VECTEUR FORTRAN
!
! Arguments d'appels
! in optmpi :
!      /'MPI_SEND' == envoyer un message mpi
!      /'MPI_RECV' == recevoir un message mpi
!
! in typsca : /'I' /'I4' /'R'
! in nudest : numero du processeur d'origine ou destinataire
! in numess : numero mpi du message
! in nbval  : longueur du vecteur vi, vr (optionnel, 1 par défaut)
!-si nbval > 1:
! inout vi(*)  : vecteur d'entiers a echanger (si typsca='I')
! inout vi4(*) : vecteur d'entiers a echanger (si typsca='I4')
! inout vr(*)  : vecteur de reels a echanger  (si typsca='R')
!-si nbval == 1:
! inout sci    : entier a echanger    (si typsca='I')
! inout sci4   : entier 4 a echanger  (si typsca='I4')
! inout scr    : réel a echanger      (si typsca='R')
!----------------------------------------------------------------------
    implicit none
! DECLARATION PARAMETRES D'APPELS
#include "asterf.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/uttcpu.h"
    character(len=*), intent(in) :: optmpi
    character(len=*), intent(in) :: typsca
    integer(kind=8), intent(in) :: nudest
    integer(kind=8), intent(in) :: numess
    integer(kind=8), intent(in), optional :: nbval
    integer(kind=8), intent(inout), optional :: vi(*)
    integer(kind=4), intent(inout), optional :: vi4(*)
    real(kind=8), intent(inout), optional :: vr(*)
    integer(kind=8), intent(inout), optional :: sci
    integer(kind=4), intent(inout), optional :: sci4
    real(kind=8), intent(inout), optional :: scr
!
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
#include "asterc/asmpi_recv_r.h"
#include "asterc/asmpi_recv_i.h"
#include "asterc/asmpi_recv_i4.h"
#include "asterc/asmpi_send_r.h"
#include "asterc/asmpi_send_i.h"
#include "asterc/asmpi_send_i4.h"
! DECLARATION VARIABLES LOCALES
    character(len=2) :: typsc1
    integer(kind=8) :: nbv
    mpi_int :: nbv4, nbpro4, nudes4, numes4
    mpi_int :: mpicou
    aster_logical :: scal
    real(kind=8) :: wkr(1)
    integer(kind=8) :: wki(1)
    integer(kind=4) :: wki4(1)
! ---------------------------------------------------------------------
    call jemarq()
!---- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicou)
! --- COMPTEUR
    call uttcpu('CPU.CMPI.1', 'DEBUT', ' ')
!
!     -- S'IL N'Y A QU'UN SEUL PROC, IL N'Y A RIEN A FAIRE :
    call asmpi_info(mpicou, size=nbpro4)
    if (nbpro4 .eq. 1) goto 999
    DEBUG_MPI('mpi_comm_point', nbpro4, ' ')
!
!     -- SCALAIRE :
!     -------------
    typsc1 = typsca
    scal = present(sci) .or. present(sci4) .or. present(scr)
    if (.not. scal) then
        ASSERT(present(nbval))
        nbv = nbval
    else
        nbv = 1
    end if
    ASSERT(typsc1 .eq. 'I' .or. typsc1 .eq. 'I4' .or. typsc1 .eq. 'R')
    ASSERT(typsc1 .ne. 'I' .or. present(vi) .or. present(sci))
    ASSERT(typsc1 .ne. 'I4' .or. present(vi4) .or. present(sci4))
    ASSERT(typsc1 .ne. 'R' .or. present(vr) .or. present(scr))
    nbv4 = nbv
    nudes4 = nudest
    numes4 = numess
!
    if (optmpi .eq. 'MPI_SEND') then
!     ---------------------------------
        if (scal) then
            if (typsc1 .eq. 'R') then
                wkr(1) = scr
                call asmpi_send_r(wkr, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I') then
                wki(1) = sci
                call asmpi_send_i(wki, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I4') then
                wki4(1) = sci4
                call asmpi_send_i4(wki4, nbv4, nudes4, numes4, mpicou)
            else
                ASSERT(.false.)
            end if
        else
            if (typsc1 .eq. 'R') then
                call asmpi_send_r(vr, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I') then
                call asmpi_send_i(vi, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I4') then
                call asmpi_send_i4(vi4, nbv4, nudes4, numes4, mpicou)
            else
                ASSERT(.false.)
            end if
        end if
    else if (optmpi .eq. 'MPI_RECV') then
!     ---------------------------------
        if (scal) then
            if (typsc1 .eq. 'R ') then
                call asmpi_recv_r(wkr, nbv4, nudes4, numes4, mpicou)
                scr = wkr(1)
            else if (typsc1 .eq. 'I ') then
                call asmpi_recv_i(wki, nbv4, nudes4, numes4, mpicou)
                sci = wki(1)
            else if (typsc1 .eq. 'I4') then
                call asmpi_recv_i4(wki4, nbv4, nudes4, numes4, mpicou)
                sci4 = wki4(1)
            else
                ASSERT(.false.)
            end if
        else
            if (typsc1 .eq. 'R ') then
                call asmpi_recv_r(vr, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I ') then
                call asmpi_recv_i(vi, nbv4, nudes4, numes4, mpicou)
            else if (typsc1 .eq. 'I4') then
                call asmpi_recv_i4(vi4, nbv4, nudes4, numes4, mpicou)
            else
                ASSERT(.false.)
            end if
        end if
    else
        ASSERT(.false.)
    end if
!
999 continue
! --- COMPTEUR
    call uttcpu('CPU.CMPI.1', 'FIN', ' ')
    call jedema()
#else
    character(len=1) :: kdummy
    integer(kind=8) :: idummy
    integer(kind=4) :: i4dummy
    real(kind=8) :: rdummy
!
    if (present(nbval) .and. present(vi) .and. present(vi4) .and. present(vr) .and. &
        present(sci) .and. present(sci4) .and. present(scr)) then
        kdummy = optmpi(1:1)
        kdummy = typsca(1:1)
        idummy = nudest
        idummy = numess
        idummy = nbval
        idummy = vi(1)
        i4dummy = vi4(1)
        rdummy = vr(1)
        idummy = sci
        i4dummy = sci4
        rdummy = scr
    end if
#endif
end subroutine
