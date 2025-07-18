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
subroutine thmGetBehaviour(compor, ds_thm)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/THM_type.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: compor(*)
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get parameters for behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  compor           : name of comportment definition (field)
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    ds_thm%ds_behaviour%defo = compor(DEFO)
    ds_thm%ds_behaviour%rela_thmc = compor(THMC_NAME)
    ds_thm%ds_behaviour%rela_ther = compor(THER_NAME)
    ds_thm%ds_behaviour%rela_hydr = compor(HYDR_NAME)
    ds_thm%ds_behaviour%rela_meca = compor(MECA_NAME)
    ds_thm%ds_behaviour%extern_addr = compor(MGIS_ADDR)
    read (compor(THMC_NVAR), '(I16)') ds_thm%ds_behaviour%nb_vari_thmc
    read (compor(THER_NVAR), '(I16)') ds_thm%ds_behaviour%nb_vari_ther
    read (compor(HYDR_NVAR), '(I16)') ds_thm%ds_behaviour%nb_vari_hydr
    read (compor(MECA_NVAR), '(I16)') ds_thm%ds_behaviour%nb_vari_meca
    read (compor(THMC_NUME), '(I16)') ds_thm%ds_behaviour%nume_thmc
    read (compor(MECA_NUME), '(I16)') ds_thm%ds_behaviour%nume_meca
    read (compor(THER_NUME), '(I16)') ds_thm%ds_behaviour%nume_ther
    read (compor(HYDR_NUME), '(I16)') ds_thm%ds_behaviour%nume_hydr

    read (compor(NVAR), '(I16)') ds_thm%ds_behaviour%nb_vari
!
! - For coupling law
!
    ds_thm%ds_behaviour%nb_pres = 0
    if (ds_thm%ds_behaviour%nume_thmc .eq. NO_LAW) then
        ds_thm%ds_behaviour%nb_pres = 0
    elseif (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_SATU .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ_ATM .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE) then
        ds_thm%ds_behaviour%nb_pres = 1
    else
        ds_thm%ds_behaviour%nb_pres = 2
    end if
!
    if (ds_thm%ds_behaviour%nb_pres .ge. 1) then
        ds_thm%ds_behaviour%nb_phase(1) = 1
        if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE_GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE) then
            ds_thm%ds_behaviour%nb_phase(1) = 2
        end if
        if (ds_thm%ds_behaviour%nb_pres .eq. 2) then
            ds_thm%ds_behaviour%nb_phase(2) = 1
            if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE .or. &
                ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ) then
                ds_thm%ds_behaviour%nb_phase(2) = 2
            end if
        end if
    end if
!
    if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_SATU .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. GAZ) then
        ds_thm%ds_behaviour%satur_type = SATURATED
    elseif (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE) then
        ds_thm%ds_behaviour%satur_type = SATURATED_SPEC
    elseif (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE_GAZ .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ_ATM .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE .or. &
            ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ) then
        ds_thm%ds_behaviour%satur_type = UNSATURATED
    else
        ASSERT(.false.)
    end if
!
    if (ds_thm%ds_behaviour%nume_thmc .ne. GAZ .and. &
        ds_thm%ds_behaviour%nume_thmc .ne. NO_LAW) then
        ds_thm%ds_material%l_liquid = ASTER_TRUE
    end if
    if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. GAZ .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE_GAZ .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ_ATM .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ) then
        ds_thm%ds_material%l_gaz = ASTER_TRUE
    end if
    if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_VAPE_GAZ .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE) then
        ds_thm%ds_material%l_steam = ASTER_TRUE
    end if
    if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ_VAPE .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_AD_GAZ) then
        ds_thm%ds_material%l_ad = ASTER_TRUE
    end if
    ds_thm%ds_material%l_r_gaz = ASTER_TRUE
    if (ds_thm%ds_behaviour%nume_thmc .eq. LIQU_SATU .or. &
        ds_thm%ds_behaviour%nume_thmc .eq. LIQU_GAZ_ATM) then
        ds_thm%ds_material%l_r_gaz = ASTER_FALSE
    end if
!
! - For thermic
!
    ds_thm%ds_behaviour%l_temp = ASTER_TRUE
    if (ds_thm%ds_behaviour%nume_ther .eq. NO_LAW) then
        ds_thm%ds_behaviour%l_temp = ASTER_FALSE
    end if
!
! - For mechanic
!
    ds_thm%ds_behaviour%l_meca = ASTER_TRUE
    if (ds_thm%ds_behaviour%rela_meca .eq. 'VIDE') then
        ds_thm%ds_behaviour%l_meca = ASTER_FALSE
    end if
!
! - For stress measures
!
    if (ds_thm%ds_behaviour%rela_meca .eq. 'GonfElas' .or. &
        ds_thm%ds_behaviour%rela_meca .eq. 'Barcelone') then
        ds_thm%ds_behaviour%l_stress_bishop = ASTER_FALSE
    end if
!

end subroutine
