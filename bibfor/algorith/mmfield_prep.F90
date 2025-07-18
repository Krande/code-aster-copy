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

subroutine mmfield_prep(field_in, field_out, &
                        l_update_, field_update_, alpha_, &
                        l_sort_, nb_cmp_, list_cmp_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/vtgpld.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=*), intent(in) :: field_in
    character(len=*), intent(in) :: field_out
    aster_logical, optional, intent(in) :: l_update_
    aster_logical, optional, intent(in) :: l_sort_
    integer(kind=8), optional, intent(in) :: nb_cmp_
    character(len=8), optional, intent(in):: list_cmp_(*)
    character(len=*), optional, intent(in) :: field_update_
    real(kind=8), optional, intent(in) :: alpha_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Convert and update fields
!
! --------------------------------------------------------------------------------------------------
!
! Two possibilities (EXCLUDED !):
! <1> UPDATE FIELDS
!   Update on geometry components (X,Y,Z): field_out(X,Y,Z) <= field_in(...) +
!                                                             alpha * field_update(...)
!   If alpha doesn't exist =>
!   Update on geometry components (X,Y,Z): field_out(X,Y,Z) <= field_in(...)
! <2> SORT FIELDS (select dof)
!   l_sort  => field_out(list_cmp)          <= field_in(...) x list_cmp
!
! In  l_update         : .true. to update on geometry components (X,Y,Z)
! In  l_sort           : .true. to sort on components
! In  nb_cmp           : number of components
! In  list_cmp         : list of components
! In  field_in         : name of input field
! In  field_update     : name of field to use for update
! In  field_out        : name of output field
! In  alpha            : coefficient for update field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_cmp
    aster_logical :: l_sort, l_update
    real(kind=8) :: alpha
    character(len=19) :: field_simple, field_update
!
! --------------------------------------------------------------------------------------------------
!
    field_simple = '&&MMFIELD.SIMPLE'
    l_sort = .false.
    l_update = .false.
    ASSERT(EXCLUS2(l_sort_, l_update_))
    if (present(alpha_)) then
        alpha = alpha_
    else
        alpha = 1.d0
    end if
    if (present(l_sort_)) then
        l_sort = l_sort_
    end if
    if (present(l_update_)) then
        l_update = l_update_
    end if
!
    if (l_update) then
        ASSERT(present(field_update_))
        field_update = field_update_
        if (present(alpha_)) then
            call vtgpld('CUMU', alpha, field_in, field_update, 'V', &
                        field_out)
        else
            call vtgpld('ZERO', alpha, field_in, field_update, 'V', &
                        field_out)
        end if
    end if
    if (l_sort) then
!
! ----- Count number of components to sort
!
        nb_cmp = nb_cmp_
        ASSERT(nb_cmp .gt. 0)
!
! ----- Convert
!
        call cnocns(field_in, 'V', field_simple)
!
! ----- Sort components
!
        call cnsred(field_simple, 0, [0], nb_cmp, list_cmp_, &
                    'V', field_out)
    end if
!
! - Cleaning
!
    call detrsd('CHAM_NO_S', field_simple)
!
end subroutine
