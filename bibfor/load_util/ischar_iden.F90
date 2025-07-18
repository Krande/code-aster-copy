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
function ischar_iden(v_load_info, i_load, nb_load, load_type_1, load_type_2, load_name)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
!
    aster_logical :: ischar_iden
    integer(kind=8), pointer :: v_load_info(:)
    integer(kind=8), intent(in) :: i_load
    integer(kind=8), intent(in) :: nb_load
    character(len=4), intent(in) :: load_type_1
    character(len=4), intent(in) :: load_type_2
    character(len=24), optional, intent(in) :: load_name
!
! --------------------------------------------------------------------------------------------------
!
! List of loads - Utility
!
! Return type of load - Identification
!
! --------------------------------------------------------------------------------------------------
!
! In  v_load_info    : vector of loads info
! In  i_load         : index in list of loads
! In  nb_load        : total number of loads
! In  load_type_1    : first level of type
!                'DIRI' - DIRICHLET
!                'NEUM' - NEUMANN
! In  load_type_2    : second level of type
! -> For Dirichlet loads
!                'DUAL' - AFFE_CHAR_MECA
!                'ELIM' - AFFE_CHAR_CINE
!                'DIDI' - Differential
!                'SUIV' - Undead load
!                '    ' - All types
! -> For Neumann loads
!                'ONDE' - ONDE PLANE
!                'SIGM' - SIGMA_INTERNE
!                'TARD' - ELEMENTS TARDIFS
!                'SUIV' - Undead load
!                '    ' - All types
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: load_nume_diri, load_nume_neum, jafci, ier
    aster_logical :: ldiri, lelim, ldual, ldidi, lneum
    aster_logical :: londe, lsigm, lelem, lsuiv, lpilo
!
! --------------------------------------------------------------------------------------------------
!
    ischar_iden = ASTER_FALSE
    lelim = ASTER_FALSE
    ldual = ASTER_FALSE
    ldiri = ASTER_FALSE
    ldidi = ASTER_FALSE
    lneum = ASTER_FALSE
    londe = ASTER_FALSE
    lsigm = ASTER_FALSE
    lsuiv = ASTER_FALSE
    lelem = ASTER_FALSE
    lpilo = ASTER_FALSE
!
    load_nume_diri = v_load_info(i_load+1)
    load_nume_neum = v_load_info(i_load+nb_load+1)

    if ((load_nume_diri .eq. -1) .or. (load_nume_diri .eq. -2) .or. (load_nume_diri .eq. -3)) then
        if (present(load_name)) then
            call jeexin(load_name(1:19)//'.AFCI', ier)
            if (ier .ne. 0) then
                call jeveuo(load_name(1:19)//'.AFCI', 'L', jafci)
                if (zi(jafci-1+1) .gt. 0) then
                    ldiri = ASTER_TRUE
                    lelim = ASTER_TRUE
                end if
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (load_nume_diri .eq. 1) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
    else if (load_nume_diri .eq. 2) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
    else if (load_nume_diri .eq. 3) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
    else if (load_nume_diri .eq. 4) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
        lsuiv = ASTER_TRUE
    else if (load_nume_diri .eq. 5) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
        lpilo = ASTER_TRUE
    else if (load_nume_diri .eq. 6) then
        ldiri = ASTER_TRUE
        ldual = ASTER_TRUE
        lpilo = ASTER_TRUE
    else if (load_nume_diri .eq. 0) then
        if (load_nume_neum .eq. 1) then
            lneum = ASTER_TRUE
        else if (load_nume_neum .eq. 2) then
            lneum = ASTER_TRUE
        else if (load_nume_neum .eq. 3) then
            lneum = ASTER_TRUE
        else if (load_nume_neum .eq. 4) then
            lneum = ASTER_TRUE
            lsuiv = ASTER_TRUE
        else if (load_nume_neum .eq. 5) then
            lneum = ASTER_TRUE
            lpilo = ASTER_TRUE
        else if (load_nume_neum .eq. 6) then
            lneum = ASTER_TRUE
            londe = ASTER_TRUE
        else if (load_nume_neum .eq. 8) then
            lneum = ASTER_TRUE
            lpilo = ASTER_TRUE
        else if (load_nume_neum .eq. 9) then
            lneum = ASTER_TRUE
            lpilo = ASTER_TRUE
            lsuiv = ASTER_TRUE
        else if (load_nume_neum .eq. 10) then
            lelem = ASTER_TRUE
        else if (load_nume_neum .eq. 11) then
            lneum = ASTER_TRUE
            lpilo = ASTER_TRUE
            lsuiv = ASTER_TRUE
        else if (load_nume_neum .eq. 20) then
            lneum = ASTER_TRUE
        else if (load_nume_neum .eq. 55) then
            lneum = ASTER_TRUE
            lsigm = ASTER_TRUE
        else if (load_nume_neum .eq. 0) then
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
    if (ldiri) then
        if (v_load_info(3*nb_load+i_load+3) .eq. 1) then
            ldidi = ASTER_TRUE
        end if
    end if
    if (load_type_1 .eq. 'DIRI') then
        if (ldiri) then
            if (load_type_2 .eq. 'DUAL') then
                ischar_iden = ldual
            else if (load_type_2 .eq. 'ELIM') then
                ischar_iden = lelim
            else if (load_type_2 .eq. 'DIDI') then
                ischar_iden = ldidi
            else if (load_type_2 .eq. 'SUIV') then
                ischar_iden = lsuiv
            else if (load_type_2 .eq. 'PILO') then
                ischar_iden = lpilo
            else if (load_type_2 .eq. '    ') then
                ischar_iden = ldiri
            else
                write (6, *) 'SOUTYP: ', load_type_2
                ASSERT(ASTER_FALSE)
            end if
        else if (lneum) then
            ischar_iden = ASTER_FALSE
        end if
    else if (load_type_1 .eq. 'NEUM') then
        if (lneum) then
            if (load_type_2 .eq. 'ONDE') then
                ischar_iden = londe
            else if (load_type_2 .eq. 'SIGM') then
                ischar_iden = lsigm
            else if (load_type_2 .eq. 'TARD') then
                ischar_iden = lelem
            else if (load_type_2 .eq. 'SUIV') then
                ischar_iden = lsuiv
            else if (load_type_2 .eq. '    ') then
                ischar_iden = lneum
            else if (load_type_2 .eq. 'PILO') then
                ischar_iden = lpilo
            else
                write (6, *) 'SOUTYP: ', load_type_2
                ASSERT(ASTER_FALSE)
            end if
        else if (ldiri) then
            ischar_iden = ASTER_FALSE
        end if
    else
        write (6, *) 'TYPCHA: ', load_type_1
        ASSERT(ASTER_FALSE)
    end if
end function
