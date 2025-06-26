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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romAlgoMGS(nb_mode, nb_equa, syst_type, field_iden, base, &
                      vr_mode_in, vr_mode_out, vc_mode_in, vc_mode_out)
!
    implicit none
!
#include "asterfort/assert.h"
#include "blas/zdotc.h"
#include "blas/ddot.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
!
    integer(kind=8), intent(in) :: nb_mode, nb_equa
    character(len=1), intent(in) :: syst_type
    character(len=8), intent(in) :: base
    character(len=24), intent(in) :: field_iden
    real(kind=8), pointer, optional :: vr_mode_in(:), vr_mode_out(:)
    complex(kind=8), pointer, optional :: vc_mode_in(:), vc_mode_out(:)
!
! --------------------------------------------------------------------------------------------------
!
! Greedy algorithm
!
! Orthogonalization the basis with algorithme MGS
! (PS : We do not normalize the basis here and we suppose that base is already normalized !)
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_mode             : number of empiric modes
! In  nb_equa             : number of equations (length of empiric mode)
! In  syst_type           : global type of system (real or complex)
! In  field_iden          : identificator of field (name in results datastructure)
! In  vr_mode_in          : pointer to mode before orthogonalization (real)
! In  vr_mode_out         : pointer to mode after orthogonalization (real)
! In  vc_mode_in          : pointer to mode before orthogonalization (complex)
! In  vc_mode_out         : pointer to mode after orthogonalization (complex)
!
! --------------------------------------------------------------------------------------------------
!
    complex(kind=8) :: term_c
    real(kind=8) :: term_r
    real(kind=8), pointer :: vr_mode(:) => null()
    complex(kind=8), pointer :: vc_mode(:) => null()
    character(len=19) :: mode
    integer(kind=8) :: i_mode, iret
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    term_c = dcmplx(0.d0, 0.d0)
    term_r = 0.d0
!
! - Orthogonalization the basis with algorithme MGS
!
    if (syst_type .eq. 'R') then
        vr_mode_out(1:nb_equa) = vr_mode_in(1:nb_equa)
        do i_mode = 1, nb_mode
            call rsexch(' ', base, field_iden, i_mode, mode, &
                        iret)
            call jeveuo(mode(1:19)//'.VALE', 'L', vr=vr_mode)
            b_n = to_blas_int(nb_equa)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            term_r = ddot(b_n, vr_mode, b_incx, vr_mode_in, b_incy)
            vr_mode_out(1:nb_equa) = vr_mode_out(1:nb_equa)-term_r*vr_mode(1:nb_equa)
        end do
    else if (syst_type .eq. 'C') then
        vc_mode_out(1:nb_equa) = vc_mode_in(1:nb_equa)
        do i_mode = 1, nb_mode
            call rsexch(' ', base, field_iden, i_mode, mode, &
                        iret)
            call jeveuo(mode(1:19)//'.VALE', 'L', vc=vc_mode)
            b_n = to_blas_int(nb_equa)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            term_c = zdotc(b_n, vc_mode, b_incx, vc_mode_in, b_incy)
            vc_mode_out(1:nb_equa) = vc_mode_out(1:nb_equa)-term_c*vc_mode(1:nb_equa)
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
