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

subroutine char_read_keyw(keywordfact, iocc, val_type, n_keyexcl, keywordexcl, &
                          n_max_keyword, n_keyword, keywordlist, val_nb, val_r, &
                          val_f, val_c)
!
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getmjm.h"
#include "asterfort/assert.h"
#include "asterfort/char_read_val.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
!
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    character(len=4), intent(in) :: val_type
    character(len=24), intent(in) :: keywordexcl
    integer(kind=8), intent(in) :: n_keyexcl
    integer(kind=8), intent(in) :: n_max_keyword
    integer(kind=8), intent(out) :: n_keyword
    character(len=16), intent(out) :: keywordlist(n_max_keyword)
    integer(kind=8), intent(out) :: val_nb(n_max_keyword)
    real(kind=8), intent(out) :: val_r(n_max_keyword)
    character(len=8), intent(out) :: val_f(n_max_keyword)
    complex(kind=8), intent(out) :: val_c(n_max_keyword)
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_CHAR_MECA
!
! Read keywords and their values except someones give in keywordexcl object
!
! --------------------------------------------------------------------------------------------------
!
! In  keywordfact   : factor keyword to read
! In  iocc          : factor keyword index in AFFE_CHAR_MECA
! In  val_type      : type of values (REEL, COMP or FONC) in keyword
! In  keywordexcl   : name of JEVEUX object for excluded keywords
! In  n_keyexcl     : number of excluded keywords
! In  n_max_keyword : maximum number of keywords can been read
! Out n_keyword     : number of keywords has been read
! Out keywordlist   : list of keywords has been read
! Out val_nb        : number of values in keyword
! Out val_r         : values (if real) in keyword
! Out val_f         : names of function (if function) in keyword
! Out val_c         : values (if complex) in keyword
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: keywordread(300)
    character(len=16) :: k16_dummy(300), keyword, val_t_dummy
    integer(kind=8) :: n, i_keyword, i_keyexcl
    aster_logical :: l_excl
    integer(kind=8) :: j_kexcl
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    ASSERT(n_max_keyword .eq. 300)
!
! - Initial number of keywords
!
    call getmjm(keywordfact, iocc, 0, keywordread, k16_dummy, &
                n)
    n_keyword = -n
    ASSERT(n_keyword .le. n_max_keyword)
!
! - Read all keywords
!
    call getmjm(keywordfact, iocc, n_keyword, keywordread, k16_dummy, &
                n)
!
! - Excluding some keywords
!
    if (n_keyexcl .ne. 0) call jeveuo(keywordexcl, 'L', j_kexcl)
    n = n_keyword
    n_keyword = 0
    do i_keyword = 1, n
        l_excl = .false.
        keyword = keywordread(i_keyword)
        do i_keyexcl = 1, n_keyexcl
            if (keyword .eq. zk24(j_kexcl-1+i_keyexcl)) then
                l_excl = .true.
            end if
        end do
        if (.not. l_excl) then
            n_keyword = n_keyword+1
            keywordlist(n_keyword) = keyword
        end if
    end do
!
! - Values for final keywords
!
    do i_keyword = 1, n_keyword
        keyword = keywordlist(i_keyword)
        call char_read_val(keywordfact, iocc, keyword, val_type, val_nb(i_keyword), &
                           val_r(i_keyword), val_f(i_keyword), val_c(i_keyword), val_t_dummy)
        ASSERT(val_nb(i_keyword) .eq. 1)
    end do
!
    call jedema()
end subroutine
