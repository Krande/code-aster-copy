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
subroutine dtmforc_dvis(nl_ind, sd_dtm_, sd_nl_, buffdtm, buffnl, &
                        time, step, depl, vite, fext)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmforc_dvis : Calculates the Generalized Zener discret element's force
!                at the current step (t)
!
!       nl_ind           : nonlinearity index (for sd_nl access)
!       sd_dtm_, buffdtm : dtm data structure and its buffer
!       sd_nl_ , buffnl  : nl  data structure and its buffer
!       time, step       : current t and integration dt
!       depl, vite       : structural modal displacement and velocity at "t"
!       fext             : projected total non-linear force
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dtmget.h"
#include "asterfort/fointe.h"
#include "asterfort/gloloc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/locglo.h"
#include "asterfort/rk5adp.h"
#include "asterfort/nlget.h"
#include "asterfort/tophys.h"
#include "asterfort/tophys_ms.h"
#include "asterfort/togene.h"
#include "asterfort/utmess.h"
#include "asterfort/zengen.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in) :: nl_ind
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer :: buffnl(:)
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: step
    real(kind=8), pointer :: depl(:)
    real(kind=8), pointer :: vite(:)
    real(kind=8), pointer :: fext(:)
!
!   -0.2- Local variables
    aster_logical :: multi_support
    integer(kind=8) :: i, iex, nbexci, ier, nbno
    integer(kind=8) :: ino, nbdecp, iret, start, finish
    real(kind=8) :: sina, cosa, sinb, cosb, sing
    real(kind=8) :: cosg, depglo(3), vitglo(3), deploc(6), vitloc(6)
    real(kind=8) :: dvitlo(3), flocal(3), errmax, fgloba(3), y0(4)
    real(kind=8) :: dy0(4), ldcpar(5), resu(8)
    integer(kind=8) :: ldcpai(1)
    character(len=8) :: ldcpac(1)
    character(len=8) :: sd_dtm, sd_nl, monmot, obst_typ
    character(len=19) :: nomres
!
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer :: origob(:) => null()
    real(kind=8), pointer :: coedep(:) => null()
    real(kind=8), pointer :: coevit(:) => null()
    real(kind=8), pointer :: psidel(:) => null()
    real(kind=8), pointer :: psidel1(:) => null()
    real(kind=8), pointer :: psidel2(:) => null()
    real(kind=8), pointer :: sincos_angle_a(:) => null()
    real(kind=8), pointer :: sincos_angle_b(:) => null()
    real(kind=8), pointer :: sincos_angle_g(:) => null()
    real(kind=8), pointer :: sign_dyz(:) => null()
    real(kind=8), pointer :: dplmod(:) => null()
    real(kind=8), pointer :: dplmod1(:) => null()
    real(kind=8), pointer :: dplmod2(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    character(len=8), pointer :: nofdep(:) => null()
    character(len=8), pointer :: nofvit(:) => null()
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    call nlget(sd_nl, _INTERNAL_VARS, vr=vint, buffer=buffnl)
    call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx, buffer=buffnl)
    start = vindx(nl_ind)
!
    deploc(1:6) = 0.d0
!
    call dtmget(sd_dtm, _MULTI_AP, kscal=monmot, buffer=buffdtm)
    multi_support = monmot(1:3) .eq. 'OUI'
    if (multi_support) then
        call dtmget(sd_dtm, _CALC_SD, kscal=nomres, buffer=buffdtm)
        call dtmget(sd_dtm, _NB_EXC_T, iscal=nbexci, buffer=buffdtm)
!
        call jeveuo(nomres//'.FDEP', 'L', vk8=nofdep)
        call jeveuo(nomres//'.FVIT', 'L', vk8=nofvit)
!
        AS_ALLOCATE(vr=coedep, size=nbexci)
        AS_ALLOCATE(vr=coevit, size=nbexci)
        do iex = 1, nbexci
            coedep(iex) = 0.d0
            coevit(iex) = 0.d0
            if (nofdep(iex) .ne. ' ') then
                call fointe('F', nofdep(iex), 1, ['INST'], [time], &
                            coedep(iex), ier)
            end if
            if (nofvit(iex) .ne. ' ') then
                call fointe('F', nofvit(iex), 1, ['INST'], [time], &
                            coevit(iex), ier)
            end if
        end do
    end if
!
    call nlget(sd_nl, _COOR_ORIGIN_OBSTACLE, iocc=nl_ind, vr=origob, buffer=buffnl)
    call nlget(sd_nl, _SINCOS_ANGLE_A, iocc=nl_ind, vr=sincos_angle_a, buffer=buffnl)
    call nlget(sd_nl, _SINCOS_ANGLE_B, iocc=nl_ind, vr=sincos_angle_b, buffer=buffnl)
    call nlget(sd_nl, _SINCOS_ANGLE_G, iocc=nl_ind, vr=sincos_angle_g, buffer=buffnl)
    call nlget(sd_nl, _SIGN_DYZ, iocc=nl_ind, vr=sign_dyz, buffer=buffnl)
    sina = sincos_angle_a(1)
    cosa = sincos_angle_a(2)
    sinb = sincos_angle_b(1)
    cosb = sincos_angle_b(2)
    sing = sincos_angle_g(1)
    cosg = sincos_angle_g(2)
!
    nbno = 1
    call nlget(sd_nl, _OBST_TYP, iocc=nl_ind, kscal=obst_typ, buffer=buffnl)
    call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=nl_ind, vr=dplmod1, buffer=buffnl)
    if (multi_support) call nlget(sd_nl, _PSI_DELT_NO1, vr=psidel1, buffer=buffnl)
!
    if (obst_typ(1:2) .eq. 'BI') then
        nbno = 2
        call nlget(sd_nl, _MODAL_DEPL_NO2, iocc=nl_ind, vr=dplmod2, buffer=buffnl)
        if (multi_support) call nlget(sd_nl, _PSI_DELT_NO2, vr=psidel2, buffer=buffnl)
    end if
!
    do ino = 1, nbno
!       --- Point toward the modal displacement for the concerned node / 1 or 2 /
        dplmod => dplmod1
        if (multi_support) psidel => psidel1
        if (ino .eq. 2) then
            dplmod => dplmod2
            if (multi_support) psidel => psidel2
        end if
!
!       --- Conversion of generalized displacements/velocities
!           back to the physical (global) basis
        if (multi_support) then
            call tophys_ms(dplmod, psidel, coedep, depl, depglo)
            call tophys_ms(dplmod, psidel, coevit, vite, vitglo)
        else
            call tophys(dplmod, depl, depglo)
            call tophys(dplmod, vite, vitglo)
        end if
!
!   --- Conversion of these vectors to the local basis
        call gloloc(depglo, origob, sina, cosa, sinb, &
                    cosb, sing, cosg, deploc(1+(ino-1)*3))
        call gloloc(vitglo, [0.d0, 0.d0, 0.d0], sina, cosa, sinb, &
                    cosb, sing, cosg, vitloc(1+(ino-1)*3))
    end do
!
    if (nbno .eq. 2) then
        do i = 1, 3
            dvitlo(i) = vitloc(i)-vitloc(3+i)
        end do
    else
        do i = 1, 3
            dvitlo(i) = vitloc(i)
        end do
    end if
!
!   -------------------------------------------------------------------------------------
!
!   --- Non linear behaviour along the local x-axis
!
!       System equations     : 1      2         3     4
!                       yy   : sigma, epsivisq, epsi, puiss
!
!   --- Internal variables   : 1      2         3     4
!                       vari : sigma, epsivisq, epsi, puiss
!
!   --- Retrieve the internal variables
    y0(1) = vint(start+7)
    y0(2) = vint(start+8)
    y0(3) = vint(start+9)
    y0(4) = vint(start+10)
!
    resu(1:8) = 0.d0
!
!   --- At initialization, the step is set to zero, there is no need to integrate
    if (abs(step) .le. r8prem()) then
        do i = 1, 4
            resu(i) = y0(i)
        end do
    else
!
!       --- Physical (behavior) parameters
        call nlget(sd_nl, _DISVISC_K1, iocc=nl_ind, rscal=ldcpar(1), buffer=buffnl)
        call nlget(sd_nl, _DISVISC_K2, iocc=nl_ind, rscal=ldcpar(2), buffer=buffnl)
        call nlget(sd_nl, _DISVISC_K3, iocc=nl_ind, rscal=ldcpar(3), buffer=buffnl)
        call nlget(sd_nl, _DISVISC_C, iocc=nl_ind, rscal=ldcpar(4), buffer=buffnl)
        call nlget(sd_nl, _DISVISC_A, iocc=nl_ind, rscal=ldcpar(5), buffer=buffnl)
!
!       --- Numerical parameters
        call nlget(sd_nl, _MAX_INTE, iocc=nl_ind, iscal=nbdecp, buffer=buffnl)
        call nlget(sd_nl, _RES_INTE, iocc=nl_ind, rscal=errmax, buffer=buffnl)
!
!       --- Velocity (along the local x-axis)
        dy0(1:4) = 0.d0
        dy0(3) = -dvitlo(1)
!
!       --- Runge-Kutta 5/4 integration from t -> t+dt
        iret = 0
        call rk5adp(4, ldcpar, ldcpai, ldcpac, time, &
                    step, nbdecp, errmax, y0, dy0, &
                    zengen, resu, iret)
        if (iret .ne. 0) then
            call utmess('A', 'DISCRETS_42', si=nbdecp, sr=errmax)
        end if
    end if
!
!   --- Retrieve the axial force (along the local x-axis)
    flocal(1:3) = 0.d0
    flocal(1) = resu(1)
!
!   --- Conversion to the global (physical) reference
    call locglo(flocal, sina, cosa, sinb, cosb, &
                sing, cosg, fgloba)
!
!   --- Generalized force on the first node
    call togene(dplmod1, fgloba, fext)
!   --- Generalized force on the second node
    if (nbno .eq. 2) then
        call togene(dplmod2, fgloba, fext, coef=-1.d0)
    end if
!
!   -------------------------------------------------------------------------------------------
!   --- Internal variables, storage
    if (step .gt. 0.d0) then
!
        finish = vindx(nl_ind+1)
        ASSERT((finish-start) .eq. NBVARINT_DVIS)
!
!       --- Local displacement of node 1
        vint(start) = deploc(1)
        vint(start+1) = deploc(2)
        vint(start+2) = deploc(3)
!
!       --- Local displacement of node 2
        vint(start+3) = deploc(4)
        vint(start+4) = deploc(5)
        vint(start+5) = deploc(6)
!
!       --- Axial deformation velocity
        vint(start+6) = dvitlo(1)
!
!       --- Generalized Zener discrete axial force
        vint(start+7) = resu(1)
!
!       --- Deformation velocity (Viscous)
        vint(start+8) = resu(2)
!
!       --- Deformation
        vint(start+9) = resu(3)
!
!       --- Power
        vint(start+10) = resu(4)
    end if
!
    if (multi_support) then
        AS_DEALLOCATE(vr=coedep)
        AS_DEALLOCATE(vr=coevit)
    end if
!
end subroutine
