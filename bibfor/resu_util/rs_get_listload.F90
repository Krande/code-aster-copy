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

subroutine rs_get_listload(result_, nume, list_load, iexcit)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterfort/getvid.h"
#include "asterfort/rsadpa.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=*), intent(in) :: result_
    integer(kind=8), intent(in) :: nume
    character(len=*), intent(out) :: list_load
    integer(kind=8), intent(out) :: iexcit
!
! --------------------------------------------------------------------------------------------------
!
! Results datastructure - Utility
!
! Get list of loads at index stored in results datastructure or from command file
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of results datastructure
! In  nume             : index to find in results datastructure
! Out list_load        : name of datastructure for list of loads when came from results
! Out iexcit           : where loads are coming from ?
!                      0 : From results datastructure
!                      1 : From command file
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: result, func_name, valk(4)
    integer(kind=8) :: jv_para, i_load_comm, i_load_resu, n1, n2, vali(2)
    integer(kind=8) :: nb_load_comm, nb_load, nb_load_resu
    character(len=19) :: list_load_r
    character(len=8), pointer :: v_loadc_name(:) => null()
    character(len=8), pointer :: v_loadc_func(:) => null()
    integer(kind=8), pointer :: v_loadr_info(:) => null()
    character(len=24), pointer :: v_loadr_name(:) => null()
    character(len=24), pointer :: v_loadr_func(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    result = result_
    iexcit = 0
    list_load = ' '
    list_load_r = ' '
!
! - Get from command file
!
    nb_load_comm = 0
    if (getexm('EXCIT', 'CHARGE') .eq. 1) then
        call getfac('EXCIT', nb_load_comm)
        if (nb_load_comm .ne. 0) then
            AS_ALLOCATE(vk8=v_loadc_name, size=nb_load_comm)
            AS_ALLOCATE(vk8=v_loadc_func, size=nb_load_comm)
            do i_load_comm = 1, nb_load_comm
                call getvid('EXCIT', 'CHARGE', iocc=i_load_comm, scal=v_loadc_name(i_load_comm), &
                            nbret=n1)
                call getvid('EXCIT', 'FONC_MULT', iocc=i_load_comm, scal=func_name, nbret=n2)
                if (n2 .ne. 0) then
                    v_loadc_func(i_load_comm) = func_name
                end if
            end do
        end if
    end if
    if (getexm(' ', 'CHARGE') .eq. 1) then
        call getvid(' ', 'CHARGE', nbval=0, nbret=nb_load)
        nb_load = -nb_load
        nb_load_comm = max(1, nb_load)
        AS_ALLOCATE(vk8=v_loadc_name, size=nb_load_comm)
        AS_ALLOCATE(vk8=v_loadc_func, size=nb_load_comm)
        call getvid(' ', 'CHARGE', nbval=nb_load, vect=v_loadc_name, nbret=n2)
    end if
!
! - Get from results datastructure
!
    call rsadpa(result, 'L', 1, 'EXCIT', nume, 0, sjv=jv_para)
    list_load_r = zk24(jv_para) (1:19)
!
! - Get objects
!
    nb_load_resu = 0
    if (list_load_r .ne. ' ') then
        call jeveuo(list_load_r(1:19)//'.LCHA', 'L', vk24=v_loadr_name)
        call jeveuo(list_load_r(1:19)//'.INFC', 'L', vi=v_loadr_info)
        call jeveuo(list_load_r(1:19)//'.FCHA', 'L', vk24=v_loadr_func)
        nb_load_resu = v_loadr_info(1)
        list_load = list_load_r
    end if
!
! - Some checks between stored in results datastructure or from command file
!
    if ((nb_load_comm .ne. 0) .and. (nb_load_resu .ne. 0)) then
!
! ----- Total number of loads
!
        if (nb_load_comm .ne. nb_load_resu) then
            vali(1) = nb_load_comm
            vali(2) = nb_load_resu
            call utmess('A', 'RESULT1_65', ni=2, vali=vali)
        end if
!
! ----- Name of loads
!
        do i_load_comm = 1, nb_load_comm
            do i_load_resu = 1, nb_load_resu
                if (v_loadc_name(i_load_comm) .eq. v_loadr_name(i_load_resu) (1:8)) then
                    goto 30
                end if
            end do
            call utmess('A', 'RESULT1_40')
30          continue
        end do
!
! ----- Name of functions
!
        do i_load_comm = 1, nb_load_comm
            do i_load_resu = 1, nb_load_resu
                func_name = v_loadr_func(i_load_resu) (1:8)
                if (func_name(1:2) .eq. '&&') then
                    func_name = ' '
                end if
                if (v_loadc_func(i_load_comm) .eq. func_name) then
                    goto 60
                end if
                if (func_name .eq. ' ') then
                    goto 60
                end if
            end do
            call utmess('A', 'RESULT1_41')
60          continue
        end do
!
! ----- Name of functions/loads
!
        do i_load_comm = 1, nb_load_comm
            do i_load_resu = 1, nb_load_resu
                if (v_loadc_name(i_load_comm) .eq. v_loadr_name(i_load_resu) (1:8)) then
                    func_name = v_loadr_func(i_load_resu) (1:8)
                    if (func_name(1:2) .eq. '&&') then
                        func_name = ' '
                    end if
                    if (v_loadc_func(i_load_comm) .eq. func_name) then
                        goto 95
                    end if
                    valk(1) = v_loadc_name(i_load_comm)
                    valk(2) = v_loadc_func(i_load_comm)
                    valk(3) = v_loadr_name(i_load_resu) (1:8)
                    valk(4) = v_loadr_func(i_load_resu) (1:8)
                    call utmess('A', 'RESULT1_66', nk=4, valk=valk)
                end if
            end do
95          continue
        end do
    end if
!
! - Cleaning
!
    AS_DEALLOCATE(vk8=v_loadc_name)
    AS_DEALLOCATE(vk8=v_loadc_func)
!
! - Where loads are coming from ?
!
    iexcit = 0
    if (nb_load_comm .ne. 0) then
        iexcit = 1
    end if
    if (nb_load_comm .eq. 0 .and. list_load_r(1:1) .eq. ' ') then
        iexcit = 1
    end if
!
end subroutine
