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
subroutine char_excl_keyw(keywordfact, keywordexcl, n_keyexcl, n_suffix, list_suffix)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    character(len=16), intent(in) :: keywordfact
    character(len=24), intent(in) :: keywordexcl
    integer(kind=8), intent(out) :: n_keyexcl
    integer(kind=8), intent(in), optional :: n_suffix
    character(len=8), optional, intent(in) :: list_suffix(*)
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Create list of excluded keywords for using in char_read_keyw
!
! --------------------------------------------------------------------------------------------------
!
! In  n_suffix    : number of sufixes excluded for topoaster_logical keywords
! In  list_suffix : list of sufixes excluded for topoaster_logical keywords
! In  keywordfact : factor keyword
! In  keywordexcl : name of JEVEUX object for excluded keywords
! Out n_keyexcl   : number of excluded keywords
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n_keyexcl_affe
    parameter(n_keyexcl_affe=8)
    character(len=24) :: excl_affe(n_keyexcl_affe)
    integer(kind=8) :: leng_affe(n_keyexcl_affe)
!
    integer(kind=8) :: i_keyw, i_suffix
    character(len=24) :: keyword
    character(len=8) :: suffix
    character(len=24), pointer :: p_keywordexcl(:) => null()
!
    data excl_affe/'GROUP_MA', 'MAILLE', 'GROUP_NO', 'NOEUD', &
        'SANS_GROUP_MA', 'SANS_MAILLE', 'SANS_GROUP_NO', 'SANS_NOEUD'/
    data leng_affe/8, 6, 8, 5, &
        13, 11, 13, 10/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    n_keyexcl = 1

! - Global affectation keywords - Count
    if (present(n_suffix)) then
        ASSERT(present(list_suffix))
        do i_suffix = 1, n_suffix
            suffix = list_suffix(i_suffix)
            do i_keyw = 1, n_keyexcl_affe
                keyword = excl_affe(i_keyw) (1:leng_affe(i_keyw))//suffix
                n_keyexcl = n_keyexcl+1
            end do
        end do
    else
        do i_keyw = 1, n_keyexcl_affe
            keyword = excl_affe(i_keyw)
            n_keyexcl = n_keyexcl+1
        end do
    end if

! - Other keywords - Count
    if (keywordfact .eq. 'FACE_IMPO') then
        n_keyexcl = n_keyexcl+3
    else if (keywordfact .eq. 'ARETE_IMPO') then
        n_keyexcl = n_keyexcl+1
    else if (keywordfact .eq. 'LIAISON_OBLIQUE') then
        n_keyexcl = n_keyexcl+1
    else if (keywordfact .eq. 'DDL_IMPO') then
        n_keyexcl = n_keyexcl+1
    else if (keywordfact .eq. 'TEMP_IMPO') then
! ----- Nothing else components
    else if (keywordfact .eq. 'SECH_IMPO') then
! ----- Nothing else components
    else if (keywordfact .eq. 'PRES_IMPO') then
! ----- Nothing else components
    else if (keywordfact .eq. 'DDL_POUTRE') then
        n_keyexcl = n_keyexcl+4
    else if (keywordfact .eq. 'LIAISON_SOLIDE') then
        n_keyexcl = n_keyexcl+3
    else
        ASSERT(ASTER_FALSE)
    end if

! - Create excluded keyword object
    if (n_keyexcl .ne. 0) then
!
! ----- Allocate keyword object
!
        call wkvect(keywordexcl, 'V V K24', n_keyexcl, vk24=p_keywordexcl)
!
! ----- Global affectation keywords - Affect
!
        n_keyexcl = 1
        p_keywordexcl(n_keyexcl) = 'TOUT'
        if (.not. present(n_suffix)) then
            do i_keyw = 1, n_keyexcl_affe
                keyword = excl_affe(i_keyw)
                n_keyexcl = n_keyexcl+1
                p_keywordexcl(n_keyexcl) = keyword
            end do
        else
            ASSERT(present(list_suffix))
            do i_suffix = 1, n_suffix
                suffix = list_suffix(i_suffix)
                do i_keyw = 1, n_keyexcl_affe
                    keyword = excl_affe(i_keyw) (1:leng_affe(i_keyw))//suffix
                    n_keyexcl = n_keyexcl+1
                    p_keywordexcl(n_keyexcl) = keyword
                end do
            end do
        end if
!
! ----- Other keywords - Affect
!
        if (keywordfact .eq. 'FACE_IMPO') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DNOR'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DRNOR'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DTAN'
        else if (keywordfact .eq. 'ARETE_IMPO') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DTAN'
        else if (keywordfact .eq. 'ARETE_IMPO') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DTAN'
        else if (keywordfact .eq. 'LIAISON_OBLIQUE') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'ANGL_NAUT'
        else if (keywordfact .eq. 'DDL_IMPO') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'BLOCAGE'
        else if (keywordfact .eq. 'TEMP_IMPO') then
! ----- Nothing else components
        else if (keywordfact .eq. 'SECH_IMPO') then
! ----- Nothing else components
        else if (keywordfact .eq. 'PRES_IMPO') then
! ----- Nothing else components
        else if (keywordfact .eq. 'DDL_POUTRE') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'ANGL_VRIL'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'VECT_Y'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'MAILLE_REPE'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'GROUP_MA_REPE'
        else if (keywordfact .eq. 'LIAISON_SOLIDE') then
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'TRAN'
            n_keyexcl = n_keyexcl+1
            p_keywordexcl(n_keyexcl) = 'DIST_MIN'
            n_keyexcl = n_keyexcl+1
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    call jedema()
end subroutine
