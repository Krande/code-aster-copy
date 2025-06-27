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

subroutine char_impo_bloc(nomg, istype_bloc, cmp_nb, cmp_name, cmp_index, &
                          vale_real, vale_cplx, vale_fonc)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
!
!
    character(len=8), intent(in) :: nomg
    aster_logical, intent(in) :: istype_bloc(3)
    character(len=8), intent(out) :: cmp_name(39)
    integer(kind=8), intent(out) :: cmp_index(39)
    integer(kind=8), intent(out) :: cmp_nb
    real(kind=8), intent(out) :: vale_real
    character(len=8), intent(out) :: vale_fonc
    complex(kind=8), intent(out):: vale_cplx
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_CHAR_MECA
!
! Prepare data for BLOCAGE in DDL_IMPO
!
! --------------------------------------------------------------------------------------------------
!
! In  nomg      : name of <GRANDEUR>
! In  istype_bloc : true if type of BLOCAGE exits
!!!!!!TYPE of BLOCAGE : DEPLACEMENT ROTATION TUYAU_FOURIER
! Out cmp_nb    : number of components
! Out cmp_name  : components name
! Out cmp_index : components index in <GRANDEUR>
! Out vale_real : value of all components (if real)
! Out vale_cplx : value of all components (if complex)
! Out vale_fonc : value of all components (if function)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: look_name_depla(3)
    character(len=8) :: look_name_rota(3)
    character(len=8) :: look_name_fourier(33)
    integer(kind=8) :: i_exis, i_cmp, jnom
!
    data look_name_depla/'DX', 'DY', 'DZ'/
    data look_name_rota/'DRX', 'DRY', 'DRZ'/
    data look_name_fourier/'UI2', 'UI3', 'VI2', 'VI3', 'WI2', 'WI3', 'UO2', &
        'UO3', 'VO2', 'VO3', 'WO2', 'WO3', 'UI4', 'UI5', &
        'VI4', 'VI5', 'WI4', 'WI5', 'UO4', 'UO5', 'VO4', &
        'VO5', 'WO4', 'WO5', 'UI6', 'UO6', 'VI6', &
        'VO6', 'WI6', 'WO6', 'WO', 'WI1', 'WO1'/
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ASSERT(istype_bloc(1) .or. istype_bloc(2) .or. istype_bloc(3))
    vale_real = 0.d0
    vale_cplx = (0.d0, 0.d0)
    vale_fonc = '&FOZERO'
    cmp_nb = 0
!
! - Information about <GRANDEUR>
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomg), 'L', jnom)

! - DOF index in catalogue

!!  --TYPE DEPLACEMENT
    if (istype_bloc(1)) then
        do i_cmp = 1, 3
            i_exis = indik8(zk8(jnom), look_name_depla(i_cmp), 1, 3)
            if (i_exis .ne. 0) then
                cmp_nb = cmp_nb+1
                cmp_index(cmp_nb) = i_exis
                cmp_name(cmp_nb) = look_name_depla(i_cmp)
            end if
        end do
    end if

!!  --TYPE ROTATION
    if (istype_bloc(2)) then
        do i_cmp = 1, 3
            i_exis = indik8(zk8(jnom), look_name_rota(i_cmp), 1, 6)
            if (i_exis .ne. 0) then
                cmp_nb = cmp_nb+1
                cmp_index(cmp_nb) = i_exis
                cmp_name(cmp_nb) = look_name_rota(i_cmp)
            end if
        end do
    end if

!!  --TYPE FOURIER
    if (istype_bloc(3)) then

        do i_cmp = 1, 33
            i_exis = indik8(zk8(jnom), look_name_fourier(i_cmp), 1, 120)
            if (i_exis .ne. 0) then
                cmp_nb = cmp_nb+1
                cmp_index(cmp_nb) = i_exis
                cmp_name(cmp_nb) = look_name_fourier(i_cmp)
            end if
        end do
    end if
!
    call jedema()
end subroutine
