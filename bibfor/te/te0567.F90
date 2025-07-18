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
subroutine te0567(nomopt, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jevech.h"
#include "asterfort/lcelem.h"
#include "asterfort/lcstco.h"
#include "asterfort/lcgeominit.h"
#include "asterfort/lcgeog.h"
#include "asterfort/lclaze.h"
#include "asterfort/mmmtdb.h"
#include "asterfort/lcmatr.h"
!
    character(len=16), intent(in) :: nomopt, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Option: RIGI_CONT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, ij
    integer(kind=8) :: nb_node_slav, nb_node_mast, nb_lagr, nb_dof
    integer(kind=8) :: indi_lagc(10)
    integer(kind=8) :: elem_dime
    integer(kind=8) :: jmatt
    real(kind=8) :: lagrc_curr, gap_curr, gapi
    integer(kind=8) :: indi_cont, nmcp, i_reso_geom
    aster_logical :: l_norm_smooth
    aster_logical :: l_axis, debug, l_upda_jaco
    character(len=8) :: elem_slav_code, elem_mast_code
    real(kind=8) :: elem_mast_coor(27), elem_slav_coor(27)
    real(kind=8) :: elem_mast_init(27), elem_slav_init(27)
    real(kind=8) :: elem_mast_coop(27), elem_slav_coop(27)
    real(kind=8) :: poin_inte_sl(16)
    real(kind=8) :: poin_inte_ma(16)
    integer(kind=8) :: nb_poin_inte
    character(len=8) :: elga_fami_slav, elga_fami_mast
    real(kind=8) :: mmat(55, 55)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    mmat(1:55, 1:55) = 0.d0
    elem_mast_coor(1:27) = 0.d0
    elem_mast_coop(1:27) = 0.d0
    elem_slav_coor(1:27) = 0.d0
    elem_slav_coop(1:27) = 0.d0
    debug = ASTER_FALSE
    ASSERT(nomopt .eq. 'RIGI_CONT')
!
! - Get informations about contact element
!
    call lcelem(nomte, elem_dime, l_axis, &
                nb_dof, nb_lagr, indi_lagc, &
                elem_slav_code, elga_fami_slav, nb_node_slav, &
                elem_mast_code, elga_fami_mast, nb_node_mast)
    ASSERT(nb_dof .le. 55)
    ASSERT(elga_fami_slav .eq. elga_fami_mast)
!
! - Get indicators
!
    call lcstco(l_upda_jaco, l_norm_smooth, i_reso_geom, &
                lagrc_curr, gap_curr, &
                indi_cont, &
                gapi, nmcp, &
                nb_poin_inte, poin_inte_sl, poin_inte_ma)
!
! - Get initial coordinates
!
    call lcgeominit(elem_dime, &
                    nb_node_slav, nb_node_mast, &
                    elem_mast_init, elem_slav_init)
!
! - Compute updated geometry
!
    call lcgeog(elem_dime, i_reso_geom, &
                nb_lagr, indi_lagc, &
                nb_node_slav, nb_node_mast, &
                elem_mast_init, elem_slav_init, &
                elem_mast_coor, elem_slav_coor)
!
! - Compute matrix
!
    if (indi_cont .eq. 1) then
        call lcmatr(elem_dime, l_axis, l_upda_jaco, l_norm_smooth, &
                    nb_lagr, indi_lagc, elga_fami_slav, &
                    nb_node_slav, elem_slav_code, elem_slav_init, elem_slav_coor, &
                    nb_node_mast, elem_mast_code, elem_mast_init, elem_mast_coor, &
                    nb_poin_inte, poin_inte_sl, poin_inte_ma, &
                    mmat)
    elseif (indi_cont .eq. 0) then
        call lclaze(elem_dime, nb_lagr, nb_node_slav, indi_lagc, mmat)
    else
!
    end if
!
! - Write (symmetric matrix)
!
    call jevech('PMATUUR', 'E', jmatt)
    do j = 1, nb_dof
        do i = 1, j
            ij = (j-1)*j/2+i
            zr(jmatt+ij-1) = mmat(i, j)
            if (debug) then
                call mmmtdb(mmat(i, j), 'IJ', i, j)
            end if
        end do
    end do
!
    call jedema()
end subroutine
