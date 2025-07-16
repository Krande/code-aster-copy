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
subroutine nmprdc(ds_algopara, nume_dof, disp_prev, sddisc, nume_inst, &
                  incr_esti, disp_esti)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/copisd.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsinch.h"
#include "asterfort/utmess.h"
#include "asterfort/chamnoIsSame.h"
#include "asterfort/vtcopy.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: disp_prev
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: incr_esti
    character(len=19), intent(in) :: disp_esti
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm - Euler prediction
!
! DEPL_CALCULE option
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_algopara      : datastructure for algorithm parameters
! In  nume_dof         : name of numbering (NUME_DDL)
! In  disp_prev        : previous displacement (T-)
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current time step
! In  incr_esti        : name of increment estimation field
! In  disp_esti        : name of displacement estimation field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nb_equa, iret
    real(kind=8) :: time
    character(len=19) :: disp_extr
    character(len=8) :: result_extr
    real(kind=8), pointer :: v_disp_esti(:) => null()
    real(kind=8), pointer :: v_disp_prev(:) => null()
    real(kind=8), pointer :: v_incr_esti(:) => null()
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: prec = 1.0d-10
    character(len=8), parameter :: crit = 'ABSOLU'
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... PAR DEPL. CALCULE'
    end if
!
! - Initializations
!
    call dismoi('NB_EQUA', nume_dof, 'NUME_DDL', repi=nb_equa)
    time = diinst(sddisc, nume_inst)
!
! - Get results datastructure for PREDICTION='DEPL_CALCULE
!
    result_extr = ds_algopara%result_prev_disp
!
! - Get displacement in results datastructure
!
    disp_extr = '&&NMPRDC.DEPEST'
    call rsinch(result_extr, 'DEPL', 'INST', time, disp_extr, &
                'EXCLU', 'EXCLU', 0, 'V', prec, crit, iret)
    if (iret .gt. 0) then
        call utmess('F', 'MECANONLINE2_27', sk=result_extr, sr=time)
    end if

! - Copy displacement
    if (nume_inst .eq. 1) then
        call vtcopy(disp_extr, disp_esti, iret)
        if (iret .ne. 0) then
            call utmess('F', 'MECANONLINE2_29')
        end if
    else
        call chamnoIsSame(disp_extr, disp_esti, iret)
        if (iret .gt. 0) then
            call utmess('F', 'MECANONLINE2_28', sr=time)
        else
            call copisd('CHAMP_GD', 'V', disp_extr, disp_esti)
        end if
    end if
!
! - Compute increment: incr_esti = disp_esti - disp_prev
!
    call jeveuo(disp_esti(1:19)//'.VALE', 'L', vr=v_disp_esti)
    call jeveuo(disp_prev(1:19)//'.VALE', 'L', vr=v_disp_prev)
    call jeveuo(incr_esti(1:19)//'.VALE', 'E', vr=v_incr_esti)
    b_n = to_blas_int(nb_equa)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, v_disp_esti, b_incx, v_incr_esti, b_incy)
    b_n = to_blas_int(nb_equa)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, v_disp_prev, b_incx, v_incr_esti, &
               b_incy)
!
end subroutine
