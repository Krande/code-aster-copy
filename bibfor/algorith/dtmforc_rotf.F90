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

subroutine dtmforc_rotf(nl_ind, sd_dtm_, sd_nl_, buffdtm, buffnl, &
                        time, depl, fext)
    implicit none
!
!
! dtmforc_rotf : Calculates the CRACKED_ROTOR forces at the current step (t)
!
!       nl_ind           : nonlinearity index (for sd_nl access)
!       sd_dtm_, buffdtm : dtm data structure and its buffer
!       sd_nl_ , buffnl  : nl  data structure and its buffer
!       time             : current instant t
!       depl             : structural modal displacement and velocity at "t"
!       fext             : projected total non-linear force
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8rddg.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dtmget.h"
#include "asterfort/fointe.h"
#include "asterfort/gloloc.h"
#include "asterfort/locglo.h"
#include "asterfort/nlget.h"
#include "asterfort/tophys.h"
#include "asterfort/togene.h"

!
!   -0.1- Input/output arguments
    integer(kind=8), intent(in)  :: nl_ind
    character(len=*), intent(in)  :: sd_dtm_
    character(len=*), intent(in)  :: sd_nl_
    integer(kind=8), pointer  :: buffdtm(:)
    integer(kind=8), pointer  :: buffnl(:)
    real(kind=8), intent(in)  :: time
    real(kind=8), pointer  :: depl(:)
    real(kind=8), pointer :: fext(:)
!
!   -0.2- Local variables
    integer(kind=8)           :: i, ier, nbno, nbmode, start
    integer(kind=8)           :: finish
    real(kind=8)      :: sina, cosa, sinb, cosb, sing
    real(kind=8)      :: cosg, depglo1(3), depglo2(3), drl(3), drg(3)
    real(kind=8)      :: origob(3), eps, ml(3), mg(3), angini
    real(kind=8)      :: angrot, fkphi, dfkphi, phi, vrotat

    character(len=3)  :: vitvar
    character(len=8)  :: sd_dtm, sd_nl, fk, dfk, foncp
!
    integer(kind=8), pointer :: vindx(:) => null()
    real(kind=8), pointer :: sincos_angle_a(:) => null()
    real(kind=8), pointer :: sincos_angle_b(:) => null()
    real(kind=8), pointer :: dplmod1(:) => null()
    real(kind=8), pointer :: dplmod2(:) => null()
    real(kind=8), pointer :: vint(:) => null()
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_nl = sd_nl_
!
    call nlget(sd_nl, _INTERNAL_VARS, vr=vint, buffer=buffnl)
    call nlget(sd_nl, _INTERNAL_VARS_INDEX, vi=vindx, buffer=buffnl)
    start = vindx(nl_ind)
!
    eps = r8prem()
!
    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
!
    call dtmget(sd_dtm, _VITE_VAR, kscal=vitvar, buffer=buffdtm)
    if (vitvar .eq. 'OUI') then
        call nlget(sd_nl, _ANG_ROTA, iocc=nl_ind, kscal=foncp, buffer=buffnl)
        call fointe('F ', foncp, 1, ['INST'], [time], &
                    angrot, ier)
    else
        call nlget(sd_nl, _ANG_INIT, iocc=nl_ind, rscal=angini, buffer=buffnl)
        call dtmget(sd_dtm, _V_ROT, rscal=vrotat, buffer=buffdtm)
        angrot = angini+vrotat*time
    end if
!
    sing = sin(angrot)
    cosg = cos(angrot)
!
!
    call nlget(sd_nl, _SINCOS_ANGLE_A, iocc=nl_ind, vr=sincos_angle_a, buffer=buffnl)
    call nlget(sd_nl, _SINCOS_ANGLE_B, iocc=nl_ind, vr=sincos_angle_b, buffer=buffnl)
!
    sina = sincos_angle_a(1)
    cosa = sincos_angle_a(2)
    sinb = sincos_angle_b(1)
    cosb = sincos_angle_b(2)
!
    call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=nl_ind, vr=dplmod1, buffer=buffnl)
    call nlget(sd_nl, _MODAL_DEPL_NO2, iocc=nl_ind, vr=dplmod2, buffer=buffnl)

    nbno = 2

    origob(1) = 0.d0
    origob(2) = 0.d0
    origob(3) = 0.d0

    call tophys(dplmod1, depl, depglo1)
    call tophys(dplmod2, depl, depglo2)

    do i = 1, 3
        drg(i) = depglo2(i)-depglo1(i)
    end do

    call gloloc(drg, origob, sina, cosa, sinb, &
                cosb, sing, cosg, drl)

    phi = atan2(drl(2), drl(3))
    if (phi .lt. 0.d0) phi = r8depi()+phi

!
    call nlget(sd_nl, _ROTR_FK, iocc=nl_ind, kscal=fk, buffer=buffnl)
    call nlget(sd_nl, _ROTR_DFK, iocc=nl_ind, kscal=dfk, buffer=buffnl)

    call fointe('F', fk, 1, ['ABSC'], [phi], &
                fkphi, ier)
    call fointe('F', dfk, 1, ['ABSC'], [phi], &
                dfkphi, ier)
!
    ml(1) = 0.d0
    ml(2) = fkphi*drl(2)+0.5d0*dfkphi*drl(3)
    ml(3) = fkphi*drl(3)-0.5d0*dfkphi*drl(2)
!
!
    call locglo(ml, sina, cosa, sinb, cosb, &
                sing, cosg, mg)
!
!       --- Generalized force on the first node
    call togene(dplmod1, mg, fext)
!       --- Generalized force on the second node
    call togene(dplmod2, mg, fext, coef=-1.d0)
!
! --------------------------------------------------------------------------------------------------
!   --- Internal variables, storage
!
    finish = vindx(nl_ind+1)
    ASSERT((finish-start) .eq. NBVARINT_ROTF)

!   --- Angle, in degrees
    vint(start) = phi*r8rddg()

!   --- Force (local)
    vint(start+1) = ml(2)
    vint(start+2) = ml(3)

end subroutine
