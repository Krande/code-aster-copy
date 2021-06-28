! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine getExternalStateVariable(rela_comp    , rela_code_py,&
                                    l_mfront_offi, l_mfront_proto,&
                                    cptr_nbvarext, cptr_namevarext,&
                                    variExteCode)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "asterf_types.h"
#include "asterc/lcextevari.h"
#include "asterc/lcinfo.h"
#include "asterc/mfront_get_external_state_variable.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/iscode.h"
#include "asterfort/utmess.h"
!
character(len=16), intent(in) :: rela_comp, rela_code_py
aster_logical, intent(in) :: l_mfront_offi, l_mfront_proto
integer, intent(in) :: cptr_nbvarext, cptr_namevarext
integer, intent(out) :: variExteCode(2)
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
! In  cptr_nbvarext    : pointer to number of external state variable
! In  cptr_namevarext  : pointer to name of external state variable
! Out variExteCode     : coded integers for external state variable
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_exte, i_exte, idummy1, idummy2, i_exte_list
    integer, parameter :: nb_exte_list = 32
    character(len=8) :: name_exte(8)

    integer :: tabcod(60)
    character(len=16), parameter :: name_varc(nb_exte_list)  = (/'ELTSIZE1','ELTSIZE2','COORGA  ',&
                                                                 'GRADVELO','HYGR    ','NEUT1   ',&
                                                                 'NEUT2   ','TEMP    ','DTX     ',&
                                                                 'DTY     ','DTZ     ','X       ',&
                                                                 'Y       ','Z       ','SECH    ',&
                                                                 'HYDR    ','CORR    ','IRRA    ',&
                                                                 'EPSAXX  ','EPSAYY  ','EPSAZZ  ',&
                                                                 'EPSAXY  ','EPSAXZ  ','EPSAYZ  ',&
                                                                 'PFERRITE','PPERLITE','PBAINITE',&
                                                                 'PMARTENS','ALPHPUR ','ALPHBET ',&
                                                                 'TIME    ','TEMPREFE'/)
    aster_logical, parameter :: l_allow_mfront(nb_exte_list) = (/.true.    ,.false.   ,.false.   ,&
                                                                 .false.   ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true.    ,.true.    ,&
                                                                 .true.    ,.true./)
!
! --------------------------------------------------------------------------------------------------
!
    variExteCode = 0

! - Get names of external state variables
    nb_exte = 0
    name_exte = ' '
    if (l_mfront_proto .or. l_mfront_offi) then
        call mfront_get_external_state_variable(cptr_nbvarext, cptr_namevarext,&
                                                name_exte    , nb_exte)
        ASSERT(nb_exte .le. 8)
    else
        call lcinfo(rela_code_py, idummy1, idummy2, nb_exte)
        ASSERT(nb_exte .le. 8)
        call lcextevari(rela_code_py, nb_exte, name_exte)
    endif

! - Print
    if (nb_exte .gt. 0) then
        call utmess('I', 'COMPOR4_21', si = nb_exte, sk = rela_comp)
        do i_exte = 1, nb_exte
            call utmess('I', 'COMPOR4_22', si = i_exte, sk = name_exte(i_exte))
        end do 
    endif

! - Coding
    tabcod = 0
    do i_exte = 1, nb_exte
        do i_exte_list = 1, nb_exte_list
            if (name_exte(i_exte) .eq. name_varc(i_exte_list)) then
                tabcod(i_exte_list) = 1
                if (.not. l_allow_mfront(i_exte_list) .and.&
                    (l_mfront_proto .or. l_mfront_offi)) then
                    call utmess('I', 'COMPOR2_25', sk = name_exte(i_exte))
                    tabcod(i_exte_list) = 0
                endif
            endif
        end do
    end do 
    call iscode(tabcod, variextecode, 60)
!
end subroutine
