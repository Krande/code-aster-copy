! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine load_neut_evol(nb_type_neumz, type_calc, time_curr, load_name, load_type_ligr, &
                          load_opti_r, load_para_r, load_obje, nb_obje)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/load_neut_data.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsinch.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer, intent(in) :: nb_type_neumz
    character(len=4), intent(in) :: type_calc
    real(kind=8), intent(in) :: time_curr
    character(len=8), intent(in) :: load_name
    character(len=6), intent(out) :: load_type_ligr
    character(len=16), intent(out) :: load_opti_r
    character(len=8), intent(out) :: load_para_r(2)
    character(len=19), intent(out) :: load_obje(2)
    integer, intent(out) :: nb_obje
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation - Thermic
!
! Get information about load for EVOL_CHAR (Neumann)
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_type_neumz    : maximum number of Neumann load type
! In  type_calc        : type of option to compute
!                        '2MBR' for second member (vector)
!                        'RESI' for residual (vector)
!                        'MRIG' for rigidity (matrix)
!                        'MTAN' for tangent matrix
! In  time_curr        : current time
! In  load_name        : name of current load
! Out load_type_ligr   : type of LIGREL for current load
! Out load_opti_r      : option for real parameter
! Out load_para_r      : name of parameterS (real)
! Out load_obje        : name of objectS (cart in AFFE_CHAR_THER)
! Out nb_obje          : number of objects
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter(nb_type_neum=11)
!
    integer :: i_type_neum, i_type_echa, iret, nb_cham
    character(len=24) :: load_keyw, evol_obje
    character(len=19) :: load_name_evol(2)
    character(len=16) :: type_sd
    character(len=8) :: evol_char
    character(len=8), pointer :: p_object(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    i_type_echa = 0
    ASSERT(nb_type_neumz .eq. nb_type_neum)
!
! - Get EVOL_CHAR
!
    evol_obje = load_name//'.CHTH.EVOL.CHAR'
    call jeexin(evol_obje, iret)
    ASSERT(iret .gt. 0)
    call jeveuo(evol_obje, 'L', vk8=p_object)
    evol_char = p_object(1)
!
! - Checks
!
    call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nb_cham)
    ASSERT(nb_cham .gt. 0)
    call gettco(evol_char, type_sd)
    ASSERT(type_sd .eq. 'EVOL_CHAR')
!
! - Identify type of load (ECHANGE with COEF_H and T_EXT or FLUX_REP with FLUN)
!
    load_name_evol(1) = '&&NTDEPR.FLURE'
    call rsinch(evol_char, 'FLUN', 'INST', time_curr, load_name_evol(1), &
                'EXCLU', 'EXCLU', 0, 'V', iret)
    if (iret .gt. 2) then

        load_name_evol(1) = '&&NTDEPR.COEFH'
        load_name_evol(2) = '&&NTDEPR.T_EXT'
!
! -     Get exterior temperature
!
        call rsinch(evol_char, 'COEF_H', 'INST', time_curr, load_name_evol(1), &
                    'EXCLU', 'EXCLU', 0, 'V', iret)
        if (iret .gt. 2) then
            call utmess('F', 'CHARGES3_12', sk=evol_char, sr=time_curr)
        end if
!
! -     Get exchange coefficient
!
        call rsinch(evol_char, 'T_EXT', 'INST', time_curr, load_name_evol(2), &
                    'EXCLU', 'EXCLU', 0, 'V', iret)
        if (iret .gt. 2) then
            call utmess('F', 'CHARGES3_12', sk=evol_char, sr=time_curr)
        end if
!
! -     Identify ECHANGE load
!
        do i_type_neum = 1, nb_type_neum
            call load_neut_data(i_type_neum, nb_type_neum, &
                                load_keyw_=load_keyw)
            if (load_keyw .eq. 'ECHANGE') then
                i_type_echa = i_type_neum
            end if
        end do
    else
!
! -     Identify FLUX_REP_NORM load
!
        do i_type_neum = 1, nb_type_neum
            call load_neut_data(i_type_neum, nb_type_neum, &
                                load_keyw_=load_keyw)
            if (load_keyw .eq. 'FLUX_REP_XYZ') then
                i_type_echa = i_type_neum
            end if
        end do
    end if

    ASSERT(i_type_echa .ne. 0)
!
! - Get data for select load
!
    call load_neut_data(i_type_echa, nb_type_neum, type_calc, &
                        load_type_ligr, &
                        load_opti_r_=load_opti_r, &
                        load_para_r_=load_para_r, &
                        load_obje_=load_obje, nb_obje_=nb_obje)
!
! - Output objects
!
    ASSERT(nb_obje .gt. 0 .and. nb_obje .le. 2)

    load_obje(1) = load_name_evol(1)
    if (nb_obje .eq. 2) load_obje(2) = load_name_evol(2)
!
end subroutine
