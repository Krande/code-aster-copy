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
subroutine dtminit(sd_dtm_, sd_int_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtminit : Initialize a transitory calculation over a modal basis.
!
#include "jeveux.h"
#include "asterfort/dtmacce.h"
#include "asterfort/dtmarch.h"
#include "asterfort/dtmbuff.h"
#include "asterfort/dtmforeq.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtminivec.h"
#include "asterfort/dtmsav.h"
#include "asterfort/dtmupmat.h"
#include "asterfort/intbuff.h"
#include "asterfort/intinivec.h"
#include "asterfort/intsav.h"
#include "asterfort/mdinit.h"
#include "asterfort/nlget.h"
#include "asterfort/nlinivec.h"
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
!
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
!
!   -0.2- Local variables
    aster_logical :: reuse
    integer(kind=8) :: nbrede, nbrevi, nbsauv, nbnli, nbpal
    integer(kind=8) :: nbmode, iret, jmass, nbvint, jchor
    integer(kind=8) :: nltreat, append, appendind, i
!
    real(kind=8) :: tinit, dt
    character(len=8) :: sd_dtm, sd_int, nomres, basemo, sd_nl
    character(len=24) :: matasm, solver
!
    integer(kind=8), pointer :: vindx(:) => null()
    integer(kind=8), pointer :: matdesc(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    real(kind=8), pointer :: depl(:) => null()
    real(kind=8), pointer :: vite(:) => null()
    real(kind=8), pointer :: acce(:) => null()
    real(kind=8), pointer :: fext(:) => null()
    integer(kind=8), pointer :: buffdtm(:) => null()
    integer(kind=8), pointer :: buffint(:) => null()
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_int = sd_int_
!
!   1 - Retrieval of the necessary information
!
    call dtmget(sd_dtm, _ARCH_NB, iscal=nbsauv)
    call dtmget(sd_dtm, _CALC_SD, kscal=nomres)
    call dtmget(sd_dtm, _INST_INI, rscal=tinit)
    call dtmget(sd_dtm, _BASE_MOD, kscal=basemo)
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode)
    call dtmget(sd_dtm, _DT, rscal=dt)
!
    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)
    nbrede = 0
    nbrevi = 0
    if (nbnli .gt. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)
        call nlget(sd_nl, _NB_PALIE, iscal=nbpal)
        call nlget(sd_nl, _NB_REL_FX, iscal=nbrede)
        call nlget(sd_nl, _NB_REL_FV, iscal=nbrevi)
!
        call nlinivec(sd_nl, _F_TOT_WK, nbmode)
        call nlinivec(sd_nl, _F_TAN_WK, nbmode)
    end if
!
    call intsav(sd_int, TIME, 1, iocc=1, rscal=tinit)
    call intsav(sd_int, INDEX, 1, iocc=1, iscal=0)
    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt)
!
!   --- Update matrices (Gyroscropy)
!                                            /--- set to 1 <=> updates matrices for ind=1
    call intsav(sd_int, IND_ARCH, 1, iscal=1)
!
    call dtmbuff(sd_dtm, buffdtm)
    call intbuff(sd_int, buffint)
    call dtmupmat(sd_dtm, sd_int, buffdtm, buffint)
!
!   --- Factorize the mass M-matrix
    call dtmget(sd_dtm, _MAT_DESC, vi=matdesc)
    if (matdesc(1) .ne. 0) then
        matasm = zk24(zi(matdesc(1)+1)) (1:19)
        call dtmget(sd_dtm, _SOLVER, savejv=solver)
        call preres(solver, 'V', iret, '&&DTMCAL.MATPRE', matasm, &
                    iret, -9999)
    end if
!
    call dtmget(sd_dtm, _MASS_FUL, lonvec=iret)
    if (iret .gt. 0) then
!       --- Copy the full mass matrix
        call intinivec(sd_int, MASS_FUL, nbmode*nbmode, iocc=1, address=jmass)
        call dtmget(sd_dtm, _MASS_FUL, rvect=zr(jmass))
!       --- Copy the factorized mass matrix
        call intinivec(sd_int, MASS_FAC, nbmode*nbmode, iocc=1, address=jmass)
        call dtmget(sd_dtm, _MASS_FAC, rvect=zr(jmass))
    end if
    call intinivec(sd_int, MASS_DIA, nbmode, iocc=1, address=jmass)
    call dtmget(sd_dtm, _MASS_DIA, rvect=zr(jmass))
!
!   --- Initialize displacements, velocities, and accelerations
    call intinivec(sd_int, DEPL, nbmode, iocc=1, vr=depl)
    call intinivec(sd_int, VITE, nbmode, iocc=1, vr=vite)
    call intinivec(sd_int, ACCE, nbmode, iocc=1, vr=acce)
    call intinivec(sd_int, FORCE_EX, nbmode, iocc=1, vr=fext)
!
!   --- Call mdinit to retrieve the initial conditions into depl, vite, and acce.
!       in addition to internal variables, initial time, etc.
    nbvint = 0
    if (nbnli .gt. 0) then
        call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx)
        nbvint = vindx(nbnli+1)-1
        call nlinivec(sd_nl, _INTERNAL_VARS, nbvint, vr=vint)
    end if
!
    call mdinit(basemo, nbmode, nbnli, depl, vite, &
                vint, iret, tinit, reprise=reuse, accgen=acce, &
                index=appendind)
!
    call dtmget(sd_dtm, _APPND_SD, iscal=append)
    if (append .ne. 0) then
        call dtmsav(sd_dtm, _APPND_SD, 1, iscal=appendind-1)
    end if
!
    if (iret .ne. 0) call utmess('F', 'ALGORITH5_24')
!
    if (nbnli .ne. 0) then
!   --- Non linear case with implicit treatment of choc-type non-linearities
!       initialize the displacement to the user given value
        call dtmget(sd_dtm, _NL_TREAT, iscal=nltreat)
    end if
!
!
!   --- Copy the retrieved internal variables to the work object zr(jchor)
    if (nbvint .gt. 0) then
        call dtminivec(sd_dtm, _NL_SAVES, nbvint, address=jchor)
        do i = 1, nbvint
            zr(jchor+i-1) = vint(i)
        end do
    end if
!
!   Calculate the initial acceleration
!
!   --- If no reuse, then acce0 must be calculated based on depl0, vite0, and the
!       forces equilibrium equation (using dtmacce). To do so, delete acce0 for now.
    call dtminivec(sd_dtm, _ACC_WORK, nbmode, address=iret)
    nullify (buffdtm)
    nullify (buffint)
    call dtmbuff(sd_dtm, buffdtm)
    call intbuff(sd_int, buffint)
!
    if (.not. (reuse)) then
        call dtmacce(sd_dtm, sd_int, 1, buffdtm, buffint)
    else
        call dtmforeq(sd_dtm, sd_int, 1, buffdtm, buffint)
    end if
!
!   Initialize archiving indices and archive the initial state
    call dtmsav(sd_dtm, _ARCH_STO, 4, ivect=[0, 0, 0, 0])
    call dtmarch(sd_dtm, sd_int, buffdtm, buffint)
!
end subroutine
