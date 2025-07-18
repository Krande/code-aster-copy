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
subroutine mmGetAlgo(l_large_slip, ndexfr, jeusup, lambds, &
                     ialgoc, ialgof, i_reso_fric, i_reso_geom, &
                     l_pena_cont, l_pena_fric, &
                     lambds_prev_, jeu_prev_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
!
    aster_logical, intent(out) :: l_large_slip
    integer(kind=8), intent(out) :: ndexfr
    real(kind=8), intent(out) :: jeusup
    real(kind=8), intent(out) :: lambds
    integer(kind=8), intent(out) :: ialgoc, i_reso_fric, ialgof, i_reso_geom
    aster_logical, intent(out) :: l_pena_cont, l_pena_fric
    real(kind=8), optional, intent(out) :: lambds_prev_, jeu_prev_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Get algorithms
!
! --------------------------------------------------------------------------------------------------
!
! Out l_large_slip     : flag for GRAND_GLISSEMENT
! Out ndexfr           : integer for EXCL_FROT_* keyword
! Out jeusup           : gap from DIST_ESCL/DIST_MAIT
! Out lambds           : contact pressure (fixed trigger)
! Out ialgoc           : formulation for contact
!                        1 - Standard
!                        3 - Penalization
!                        5 - LAC
! Out ialgof           : formulation for friction
!                        1 - Standard
!                        3 - Penalization
! Out i_reso_fric      : algorithm for friction
! Out i_reso_geom      : algorithm for geometry
! Out l_pena_cont      : flag for penalized contact
! Out l_pena_fric      : flag for penalized friction
! Out lambds_prev      : contact pressure from previous iteration (fixed trigger)
! Out jeu_prev         : gap from previous iteration
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jpcf
    real(kind=8) :: lambds_prev, jeu_prev
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PCONFR', 'L', jpcf)
!
    lambds = zr(jpcf-1+13)
    jeusup = zr(jpcf-1+14)
    ialgoc = nint(zr(jpcf-1+15))
    i_reso_fric = nint(zr(jpcf-1+17))
    ialgof = nint(zr(jpcf-1+18))
    ndexfr = nint(zr(jpcf-1+21))
    i_reso_geom = nint(zr(jpcf-1+25))
    l_large_slip = nint(zr(jpcf-1+48)) .eq. 1
    lambds_prev = zr(jpcf-1+26)
    jeu_prev = zr(jpcf-1+29)
    l_pena_cont = (ialgoc .eq. 3) .or. nint(zr(jpcf-1+45)) .eq. 4
    l_pena_fric = (ialgof .eq. 3) .or. nint(zr(jpcf-1+46)) .eq. 4
    if (present(lambds_prev_)) lambds_prev_ = lambds_prev
    if (present(jeu_prev_)) jeu_prev_ = jeu_prev
!
end subroutine
