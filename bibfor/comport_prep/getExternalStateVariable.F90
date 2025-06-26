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
subroutine getExternalStateVariable(rela_comp, rela_code_py, &
                                    l_mfront_offi, l_mfront_proto, &
                                    extern_addr, variExteCode)
!
    use NonLin_Datastructure_type
    use Behaviour_module
!
    implicit none
!
#include "asterc/lcextevari.h"
#include "asterc/lcinfo.h"
#include "asterc/mgis_get_esvs.h"
#include "asterc/mgis_get_number_of_esvs.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/iscode.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: rela_comp, rela_code_py
    aster_logical, intent(in) :: l_mfront_offi, l_mfront_proto
    character(len=16), intent(in) :: extern_addr
    integer(kind=8), intent(out) :: variExteCode(2)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get external states variables
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  rela_code_py     : coded comportment for RELATION (coding in Python)
! In  l_mfront_proto   : .true. if MFront prototype
! In  l_mfront_offi    : .true. if MFront official
! In  extern_addr          : pointer to the MGIS Behaviour
! Out variExteCode     : coded integers for external state variable
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_exte, i_exte, idummy1, idummy2, i_exte_list
    integer(kind=8), parameter :: nb_exte_list = 31
    character(len=64) :: name_exte(ESVA_EXTE_NBMAXI)
    character(len=8) :: varc_aster
    integer(kind=8) :: tabcod(60)
    character(len=8), parameter :: name_varc(nb_exte_list) = (/ &
                                   'ELTSIZE1', 'COORGA  ', &
                                   'GRADVELO', 'HYGR    ', 'NEUT1   ', &
                                   'NEUT2   ', 'TEMP    ', 'DTX     ', &
                                   'DTY     ', 'DTZ     ', 'X       ', &
                                   'Y       ', 'Z       ', 'SECH    ', &
                                   'HYDR    ', 'CORR    ', 'IRRA    ', &
                                   'EPSAXX  ', 'EPSAYY  ', 'EPSAZZ  ', &
                                   'EPSAXY  ', 'EPSAXZ  ', 'EPSAYZ  ', &
                                   'PFERRITE', 'PPERLITE', 'PBAINITE', &
                                   'PMARTENS', 'ALPHPUR ', 'ALPHBET ', &
                                   'TIME    ', 'TEMPREFE'/)
    aster_logical, parameter :: l_allow_mfront(nb_exte_list) = (/.true., .false., &
                                                                 .false., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true., .true., &
                                                                 .true., .true./)
!
! --------------------------------------------------------------------------------------------------
!
    variExteCode = 0

! - Get names of external state variables
    nb_exte = 0
    name_exte = ' '
    if (l_mfront_proto .or. l_mfront_offi) then
        call mgis_get_number_of_esvs(extern_addr, nb_exte)
        ASSERT(nb_exte .le. ESVA_EXTE_NBMAXI)
        call mgis_get_esvs(extern_addr, name_exte)
    else
        call lcinfo(rela_code_py, idummy1, idummy2, nb_exte)
        ASSERT(nb_exte .le. ESVA_EXTE_NBMAXI)
        call lcextevari(rela_code_py, nb_exte, name_exte)
    end if

! - Print
    if (nb_exte .gt. 0) then
        call utmess('I', 'COMPOR4_21', si=nb_exte, sk=rela_comp)
        do i_exte = 1, nb_exte
            call utmess('I', 'COMPOR4_22', si=i_exte, sk=name_exte(i_exte))
        end do
    end if

! - Coding
    tabcod = 0
    do i_exte = 1, nb_exte
        do i_exte_list = 1, nb_exte_list
            varc_aster = getAsterVariableName(name_exte(i_exte))
            if (varc_aster .eq. name_varc(i_exte_list)) then
                tabcod(i_exte_list) = 1
                if (.not. l_allow_mfront(i_exte_list) .and. &
                    (l_mfront_proto .or. l_mfront_offi)) then
                    call utmess('I', 'COMPOR2_25', sk=name_exte(i_exte))
                    tabcod(i_exte_list) = 0
                end if
            end if
        end do
    end do
    call iscode(tabcod, variExteCode, 60)
!
end subroutine
