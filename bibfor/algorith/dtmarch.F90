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
!
subroutine dtmarch(sd_dtm_, sd_int_, buffdtm, buffint)
    use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmarch : Archive the current step.
!
#include "jeveux.h"
#include "blas/dcopy.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dtmallo.h"
#include "asterfort/dtmcase_coder.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmsav.h"
#include "asterfort/intget.h"
#include "asterfort/jelira.h"
#include "asterfort/jgetptc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nlget.h"
#include "asterfort/pmavec.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer :: buffint(:)
!
!   -0.2- Local variables
    integer(kind=8) :: ipas, nbmode, nbsauv, nbnoli, nbvint
    integer(kind=8) :: ndec
    integer(kind=8) :: ind, iret, nlcase
    integer(kind=8) :: index, iarch_sd
    real(kind=8) :: t, dt
    character(len=7) :: casek7, intk7
    character(len=8) :: sd_dtm, sd_int, sd_nl
    character(len=16) :: nomres16
    character(len=24) :: nomres
    type(c_ptr) :: pc
!
    integer(kind=8), pointer :: vindx(:) => null()
    integer(kind=8), pointer :: allocs(:) => null()
    integer(kind=8), pointer :: isto(:) => null()
    integer(kind=8), pointer :: iorsto(:) => null()
    real(kind=8), pointer :: depgen(:) => null()
    real(kind=8), pointer :: vitgen(:) => null()
    real(kind=8), pointer :: accgen(:) => null()
    real(kind=8), pointer :: depl0(:) => null()
    real(kind=8), pointer :: vite0(:) => null()
    real(kind=8), pointer :: acce0(:) => null()
    real(kind=8), pointer :: phi(:) => null()
    real(kind=8), pointer :: temsto(:) => null()
    real(kind=8), pointer :: passto(:) => null()
    real(kind=8), pointer :: depsto(:) => null()
    real(kind=8), pointer :: vitsto(:) => null()
    real(kind=8), pointer :: accsto(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    real(kind=8), pointer :: vintsto(:) => null()
    integer(kind=8), pointer :: buffnl(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
!   0 - Initializations
!
    sd_dtm = sd_dtm_
    sd_int = sd_int_
!
    call dtmget(sd_dtm, _ARCH_STO, vi=isto, buffer=buffdtm)
    call intget(sd_int, IND_ARCH, iscal=index, buffer=buffint)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
!
    call dtmget(sd_dtm, _CALC_SD, kscal=nomres, buffer=buffdtm)
    call dtmget(sd_dtm, _IARCH_SD, iscal=iarch_sd)
    if (iarch_sd .gt. 0) then
        call codent(iarch_sd, 'D0', intk7)
        nomres16 = nomres(1:8)//'.'//intk7
        call jelira(nomres16//'   .ORDR', 'LONMAX', nbsauv)
        if (isto(1) .ge. (nbsauv)) then
            iarch_sd = iarch_sd+1
            call dtmsav(sd_dtm, _IARCH_SD, 1, iscal=iarch_sd)
            call dtmallo(sd_dtm)
        end if
    else
        call dtmget(sd_dtm, _ARCH_NB, iscal=nbsauv, buffer=buffdtm)
        if (isto(1) .ge. (nbsauv)) then
            ASSERT(.false.)
        end if
    end if
!
    call intget(sd_int, STEP, iocc=index, rscal=dt, buffer=buffint)
    call intget(sd_int, INDEX, iocc=index, iscal=ipas, buffer=buffint)
    call intget(sd_int, TIME, iocc=index, rscal=t, buffer=buffint)
!
    call dtmget(sd_dtm, _NL_CASE, iscal=nlcase, buffer=buffdtm)
!
    call intget(sd_int, DEPL, iocc=index, lonvec=iret, buffer=buffint)
    if (iret .ne. 0) then
        if (nlcase .eq. 0) then
            call intget(sd_int, DEPL, iocc=index, vr=depgen, buffer=buffint)
            call intget(sd_int, VITE, iocc=index, vr=vitgen, buffer=buffint)
            call intget(sd_int, ACCE, iocc=index, vr=accgen, buffer=buffint)
        else
!           --- Implicit treatment of chocs, project the acceleration to the
!               previous basis : [Phi] x ACCE
!               Same treatment is done for displacement and velocity
!
            call dtmcase_coder(nlcase, casek7)
            call jeveuo(sd_dtm//'.PRJ_BAS.'//casek7, 'E', vr=phi)
            call intget(sd_int, DEPL, iocc=index, vr=depl0, buffer=buffint)
            call intget(sd_int, VITE, iocc=index, vr=vite0, buffer=buffint)
            call intget(sd_int, ACCE, iocc=index, vr=acce0, buffer=buffint)
            call dtmget(sd_dtm, _IMP_DEPL, vr=depgen, buffer=buffdtm)
            call dtmget(sd_dtm, _IMP_VITE, vr=vitgen, buffer=buffdtm)
            call dtmget(sd_dtm, _IMP_ACCE, vr=accgen, buffer=buffdtm)
            call pmavec('ZERO', nbmode, phi, depl0, depgen)
            call pmavec('ZERO', nbmode, phi, vite0, vitgen)
            call pmavec('ZERO', nbmode, phi, acce0, accgen)
        end if
    else
        ASSERT(.false.)
    end if
!
!
    call dtmget(sd_dtm, _IND_ALOC, vi=allocs, buffer=buffdtm)
!
!   [jordr, jdisc, jptem, jdepl , jvite, jacce,
!    jfcho, jdcho, jvcho, jadcho, jredc, jredd,
!    jrevc, jrevv, jvint                       ]
!
    call jgetptc(allocs(1)+isto(1), pc, vi=zi(1))
    call c_f_pointer(pc, iorsto, [1])
!
    call jgetptc(allocs(2)+isto(1), pc, vr=zr(1))
    call c_f_pointer(pc, temsto, [1])
!
    call jgetptc(allocs(3)+isto(1), pc, vr=zr(1))
    call c_f_pointer(pc, passto, [1])
!
    ind = nbmode*isto(1)
!
    call jgetptc(allocs(4)+ind, pc, vr=zr(1))
    call c_f_pointer(pc, depsto, [nbmode])
    call jgetptc(allocs(5)+ind, pc, vr=zr(1))
    call c_f_pointer(pc, vitsto, [nbmode])
    call jgetptc(allocs(6)+ind, pc, vr=zr(1))
    call c_f_pointer(pc, accsto, [nbmode])
!
!   Obligatory information (time, displacement, velocity, etc.)
    iorsto(1) = ipas
    temsto(1) = t
    passto(1) = dt
    b_n = to_blas_int(nbmode)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depgen, b_incx, depsto, b_incy)
    b_n = to_blas_int(nbmode)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vitgen, b_incx, vitsto, b_incy)
    b_n = to_blas_int(nbmode)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, accgen, b_incx, accsto, b_incy)
!
!   Optional information (nonlinearities)
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnoli, buffer=buffdtm)
    if (nbnoli .gt. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl, buffer=buffdtm)
        call dtmget(sd_dtm, _NL_BUFFER, vi=buffnl, buffer=buffdtm)
!
        call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx, buffer=buffnl)
        nbvint = vindx(nbnoli+1)-1
!
        ndec = nbvint*isto(1)
        call jgetptc(allocs(7)+ndec, pc, vr=zr(1))
        call c_f_pointer(pc, vintsto, [nbvint])
!
        call nlget(sd_nl, _INTERNAL_VARS, vr=vint, buffer=buffnl)
!
        b_n = to_blas_int(nbvint)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vint, b_incx, vintsto, b_incy)
    end if
!
    isto(1) = isto(1)+1
!
end subroutine
