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

subroutine dtmforc(sd_dtm_, sd_int_, index, buffdtm, buffint, nlaccnt)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmforc : Calculate the forces at the current step, specified by the argument
!           "index".
!
#include "jeveux.h"
#include "asterfort/dtmcase_coder.h"
#include "asterfort/dtmforc_calcnoli.h"
#include "asterfort/dtmget.h"
#include "asterfort/intget.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/intbuff.h"
#include "asterfort/intinivec.h"
#include "asterfort/mdfext.h"
#include "asterfort/nlget.h"
#include "asterfort/pmavec.h"
#include "asterfort/pmtvec.h"
#include "asterfort/uttcpu.h"
#include "asterfort/uttcpr.h"
#include "asterfort/sdmpic.h"
#include "asterfort/wkvect.h"
!
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_info.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
    integer(kind=8), intent(in) :: index
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer              :: buffint(:)
    integer(kind=8), optional, intent(in) :: nlaccnt

!
!   -0.2- Local variables
    aster_logical    :: instrum
    integer(kind=8)          :: nbmode, ntotex, nbnli, itime
    integer(kind=8)          :: iret, i, nlcase, nlacc
    real(kind=8)     :: temps, dt, r8b, tps(7)
    character(len=7) :: casek7
    character(len=8) :: sd_dtm, sd_int, sd_nl
    character(len=24) :: sd_fextnl_tmp
!
    integer(kind=8), pointer :: idescf(:) => null()
    integer(kind=8), pointer :: liad(:) => null()
    integer(kind=8), pointer :: inumor(:) => null()
    real(kind=8), pointer :: depl0(:) => null()
    real(kind=8), pointer :: vite0(:) => null()
    real(kind=8), pointer :: acce0(:) => null()
    real(kind=8), pointer :: fext0(:) => null()
    real(kind=8), pointer :: depl(:) => null()
    real(kind=8), pointer :: vite(:) => null()
    real(kind=8), pointer :: acce(:) => null()
    real(kind=8), pointer :: fext(:) => null()
    real(kind=8), pointer :: fext_tmp(:) => null()
    real(kind=8), pointer :: phi(:) => null()
!    real(kind=8)    , pointer :: fext_nl(:) => null()
!    real(kind=8)    , pointer :: fext_tgt(:)=> null()
    real(kind=8), pointer :: fadd_nl(:) => null()
    real(kind=8), pointer :: coefm(:) => null()
    character(len=8), pointer :: nomfon(:) => null()
    integer(kind=8), pointer          :: buffnl(:) => null()
    mpi_int :: i_proc, nb_proc, mpicou
    aster_logical :: one_proc
    integer(kind=8) :: nb_nbnli_mpi, nbr_nbnli_mpi, idx_start, idx_end
    integer(kind=8) :: iret_mpi
    character(len=16), pointer :: valk(:) => null()

!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_int = sd_int_
!
    nlacc = 0
    if (present(nlaccnt)) nlacc = nlaccnt
!
    instrum = .false.
    if (instrum) then
        call uttcpu('CPU.DTMFORC', 'INIT', ' ')
        call uttcpu('CPU.DTMFORC', 'DEBUT', ' ')
    end if

!
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
    call dtmget(sd_dtm, _NB_EXC_T, iscal=ntotex, buffer=buffdtm)
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli, buffer=buffdtm)
!

    call intget(sd_int, FORCE_EX, iocc=index, lonvec=iret, buffer=buffint)
    if (iret .eq. 0) then
        call intinivec(sd_int, FORCE_EX, nbmode, iocc=index, vr=fext)
        nullify (buffint)
        call intbuff(sd_int, buffint, level=index)
    else
        call intget(sd_int, FORCE_EX, iocc=index, vr=fext, buffer=buffint)
    end if
!
    call dtmget(sd_dtm, _NL_CASE, iscal=nlcase, buffer=buffdtm)

    if (nlcase .eq. 0) then
!       --- No implicit treatment of chocs or free state (vol), simply use the
!           integration displacements, velocities, and accelerations
        call intget(sd_int, DEPL, iocc=index, vr=depl, buffer=buffint)
        call intget(sd_int, VITE, iocc=index, vr=vite, buffer=buffint)
!       --- The acceleration at the preceding step is required
        call intget(sd_int, ACCE, iocc=1, vr=acce, buffer=buffint)
    else
!       --- Implicit treatment of chocs, project these vectors back to the original
!           basis by calculating : [Phi] x DEPL, [Phi] x VITE, [Phi] x ACCE
        call intget(sd_int, DEPL, iocc=index, vr=depl0, buffer=buffint)
        call intget(sd_int, VITE, iocc=index, vr=vite0, buffer=buffint)
        call intget(sd_int, ACCE, iocc=1, vr=acce0, buffer=buffint)
        call dtmcase_coder(nlcase, casek7)
        call jeveuo(sd_dtm//'.PRJ_BAS.'//casek7, 'E', vr=phi)
        call dtmget(sd_dtm, _IMP_DEPL, vr=depl, buffer=buffdtm)
        call dtmget(sd_dtm, _IMP_VITE, vr=vite, buffer=buffdtm)
        call dtmget(sd_dtm, _IMP_ACCE, vr=acce, buffer=buffdtm)
        call pmavec('ZERO', nbmode, phi, depl0, depl)
        call pmavec('ZERO', nbmode, phi, vite0, vite)
        call pmavec('ZERO', nbmode, phi, acce0, acce)
        nullify (fext)
        call dtmget(sd_dtm, _IMP_FEXT, vr=fext, buffer=buffdtm)
        call intget(sd_int, FORCE_EX, iocc=index, vr=fext0, buffer=buffint)
    end if

    fext(:) = 0.d0

    call intget(sd_int, INDEX, iocc=index, iscal=itime, buffer=buffint)
    call intget(sd_int, TIME, iocc=index, rscal=temps, buffer=buffint)
    call intget(sd_int, STEP, iocc=index, rscal=dt, buffer=buffint)
!
    if (ntotex .ne. 0) then
        call dtmget(sd_dtm, _DESC_FRC, vi=idescf, buffer=buffdtm)
        call dtmget(sd_dtm, _FUNC_NAM, vk8=nomfon, buffer=buffdtm)
        call dtmget(sd_dtm, _COEF_MLT, vr=coefm, buffer=buffdtm)
        call dtmget(sd_dtm, _ADRES_VC, vi=liad, buffer=buffdtm)
        call dtmget(sd_dtm, _N_ORD_VC, vi=inumor, buffer=buffdtm)
        call mdfext(temps, r8b, nbmode, ntotex, idescf, &
                    nomfon, coefm, liad, inumor, 1, &
                    fext)
    end if
!

!        print *, "-----bfore FEXT----"
!        print *, fext
!        print *, "-------------"

    if (nbnli .ne. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl, buffer=buffdtm)
        call dtmget(sd_dtm, _NL_BUFFER, vi=buffnl, buffer=buffdtm)

!        call nlget (sd_nl,  _F_TOT_WK , vr=fext_nl , buffer=buffnl)
!        call nlget (sd_nl,  _F_TAN_WK , vr=fext_tgt, buffer=buffnl)
!        call nlget(sd_nl, _FEXT_MPI, vr=fext_tmp,savejv=sd_fextnl_tmp)
!        write (6,*) "OBJET JEVEUX MPI" , sd_fextnl_tmp

        !
        ! ----- Mpi informations
        !
        one_proc = .false.
        call asmpi_comm('GET', mpicou)
        call asmpi_info(mpicou, rank=i_proc, size=nb_proc)
        if (nb_proc .eq. 1 .or. nbnli .le. nb_proc) one_proc = .true.

        if (one_proc) then
            idx_start = 1
            idx_end = nbnli
            call dtmforc_calcnoli(sd_dtm, sd_nl, buffdtm, buffnl, &
                                  temps, dt, depl, vite, fext, &
                                  nbmode, nlacc, nlcase, idx_start, idx_end)
        else
            ! ----- MPI initialisation
            ! ----- Temporary non linear force
            call nlget(sd_nl, _FEXT_MPI, savejv=sd_fextnl_tmp, vr=fext_tmp)
            fext_tmp(:) = 0.d0
            call jeexin(sd_fextnl_tmp(1:20)//'.MPI', iret_mpi)

            if (iret_mpi .eq. 0) then
                call wkvect(sd_fextnl_tmp(1:20)//'.MPI', 'V V K16', 1, vk16=valk)
            else
                call jeveuo(sd_fextnl_tmp(1:20)//'.MPI', 'E', vk16=valk)
            end if
            valk(1) = 'MPI_INCOMPLET'
            nb_nbnli_mpi = int(nbnli/nb_proc)
            nbr_nbnli_mpi = nbnli-nb_nbnli_mpi*nb_proc
            idx_start = 1+(i_proc)*nb_nbnli_mpi
            idx_end = idx_start+nb_nbnli_mpi-1+(nbr_nbnli_mpi*int((i_proc+1)/nb_proc))
            call dtmforc_calcnoli(sd_dtm, sd_nl, buffdtm, buffnl, &
                                  temps, dt, depl, vite, fext_tmp, &
                                  nbmode, nlacc, nlcase, idx_start, idx_end)
            !
            ! - All-reduce _NL_FEXT_MPI
            !
            call sdmpic('&&OP29NL', sd_fextnl_tmp)
            fext = fext+fext_tmp
        end if

        if (nlcase .ne. 0) then
!           --- Implicit treatment of chocs, project the force to the new basis
!               by calculating : [Phi]_t x FORC
            call pmtvec('ZERO', nbmode, phi, fext, fext0)
            if (nlacc .eq. 0) then
                call dtmget(sd_dtm, _F_NL_ADD, vr=fadd_nl, buffer=buffdtm)
                do i = 1, nbmode
                    fext0(i) = fext0(i)+fadd_nl(i)
                end do
            end if
        end if
    end if

    if (instrum) then
        call uttcpu('CPU.DTMFORC', 'FIN', ' ')
        call uttcpr('CPU.DTMFORC', 7, tps)
        write (*, *) "ELPSD DTMFORC : ", tps(7)
    end if

end subroutine
