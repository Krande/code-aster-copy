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
subroutine dismca(question_, object_, answeri, answerk_, ierd)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fonbpa.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/dismgd.h"
!
    character(len=*), intent(in) :: question_
    character(len=*), intent(in) :: object_
    integer(kind=8), intent(out) :: answeri, ierd
    character(len=*), intent(out)  :: answerk_
!
! --------------------------------------------------------------------------------------------------
!
! The DISMOI mechanism
!
! Questions for <CARTE> objects
!
! --------------------------------------------------------------------------------------------------
!
! In  question       : question on object
! In  object         : name of object
! Out answeri        : answer when integer
! Out answerk        : answer when string
! Out ierd           : return code for error
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24)   :: question
    character(len=32)   :: answerk
    character(len=19)   :: object, func_name, field
    integer(kind=8), parameter  :: max_para_name = 15
    character(len=8)    :: func_type, para_name(max_para_name), type, nogd
    integer(kind=8)             :: iexi, iret, repi
    integer(kind=8)             :: jvale, i_zone, i_para, nb_zone, ltyp, nb_para
    integer(kind=8), pointer            :: v_desc(:) => null()
    character(len=24), pointer  :: v_prol(:) => null()
    character(len=8), pointer   :: v_noma(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    answerk_ = ' '
    answeri = 0
    ierd = 0
    object = object_
    question = question_
!
! - Object exists ?
!
    call jeexin(object//'.NOMA', iexi)
    if (iexi .eq. 0) then
        ierd = 1
        goto 999
    end if
!
    if (question .eq. 'NOM_MAILLA') then
        call jeveuo(object//'.NOMA', 'L', vk8=v_noma)
        answerk = v_noma(1)
!
    else if (question .eq. 'TYPE_CHAMP') then
        answerk = 'CART'
!
    else if (question .eq. 'TYPE_SUPERVIS') then
        call jeveuo(object//'.DESC', 'L', vi=v_desc)
        call jenuno(jexnum('&CATA.GD.NOMGD', v_desc(1)), nogd)
        answerk = 'CART_'//nogd
!
    else if (question .eq. 'TYPE_SCA') then
        call jeveuo(object//'.DESC', 'L', vi=v_desc)
        call jenuno(jexnum('&CATA.GD.NOMGD', v_desc(1)), nogd)
        call dismgd(question, nogd, repi, answerk, ierd)
!
    else if (question(1:7) .eq. 'NOM_GD ') then
        call jeveuo(object//'.DESC', 'L', vi=v_desc)
        call jenuno(jexnum('&CATA.GD.NOMGD', v_desc(1)), answerk)
!
    else if (question .eq. 'PARA_INST') then
        answerk = ' '
        field = object
        call jeveuo(field//'.VALE', 'L', jvale)
        call jelira(field//'.VALE', 'TYPE', cval=type)
        if (type(1:1) .eq. 'K') then
            call jelira(field//'.VALE', 'LONMAX', nb_zone)
            call jelira(field//'.VALE', 'LTYP', ltyp)
            do i_zone = 1, nb_zone
                if (ltyp .eq. 8) then
                    func_name = zk8(jvale+i_zone-1)
                else if (ltyp .eq. 24) then
                    func_name = zk24(jvale+i_zone-1) (1:19)
                else
                    ASSERT(.false.)
                end if
                if (func_name(1:8) .ne. ' ') then
                    call jeexin(func_name//'.PROL', iret)
                    if (iret .gt. 0) then
                        call jeveuo(func_name//'.PROL', 'L', vk24=v_prol)
                        call fonbpa(func_name, v_prol, func_type, max_para_name, nb_para, para_name)
                        do i_para = 1, nb_para
                            if (para_name(i_para) (1:4) .eq. 'INST') then
                                answerk = 'OUI'
                                goto 999
                            end if
                        end do
                    end if
                end if
            end do
        end if
!
    else if (question .eq. 'PARA_VITE') then
        answerk = ' '
        field = object
        call jeveuo(field//'.VALE', 'L', jvale)
        call jelira(field//'.VALE', 'TYPE', cval=type)
        if (type(1:1) .eq. 'K') then
            call jelira(field//'.VALE', 'LONMAX', nb_zone)
            call jelira(field//'.VALE', 'LTYP', ltyp)
            do i_zone = 1, nb_zone
                if (ltyp .eq. 8) then
                    func_name = zk8(jvale+i_zone-1)
                else if (ltyp .eq. 24) then
                    func_name = zk24(jvale+i_zone-1) (1:19)
                else
                    ASSERT(.false.)
                end if
                if (func_name(1:8) .ne. ' ') then
                    call jeexin(func_name//'.PROL', iret)
                    if (iret .gt. 0) then
                        call jeveuo(func_name//'.PROL', 'L', vk24=v_prol)
                        call fonbpa(func_name, v_prol, func_type, max_para_name, nb_para, para_name)
                        do i_para = 1, nb_para
                            if (para_name(i_para) .eq. 'VITE_X' .or. &
                                para_name(i_para) .eq. 'VITE_Y' .or. &
                                para_name(i_para) .eq. 'VITE_Z') then
                                answerk = 'OUI'
                                goto 999
                            end if
                        end do
                    end if
                end if
            end do
        end if
!
    else if (question .eq. 'PARA_ACCE') then
        answerk = ' '
        field = object
        call jeveuo(field//'.VALE', 'L', jvale)
        call jelira(field//'.VALE', 'TYPE', cval=type)
        if (type(1:1) .eq. 'K') then
            call jelira(field//'.VALE', 'LONMAX', nb_zone)
            call jelira(field//'.VALE', 'LTYP', ltyp)
            do i_zone = 1, nb_zone
                if (ltyp .eq. 8) then
                    func_name = zk8(jvale+i_zone-1)
                else if (ltyp .eq. 24) then
                    func_name = zk24(jvale+i_zone-1) (1:19)
                else
                    ASSERT(.false.)
                end if
                if (func_name(1:8) .ne. ' ') then
                    call jeexin(func_name//'.PROL', iret)
                    if (iret .gt. 0) then
                        call jeveuo(func_name//'.PROL', 'L', vk24=v_prol)
                        call fonbpa(func_name, v_prol, func_type, max_para_name, nb_para, para_name)
                        do i_para = 1, nb_para
                            if (para_name(i_para) .eq. 'ACCE_X' .or. &
                                para_name(i_para) .eq. 'ACCE_Y' .or. &
                                para_name(i_para) .eq. 'ACCE_Z') then
                                answerk = 'OUI'
                                goto 999
                            end if
                        end do
                    end if
                end if
            end do
        end if
    else if (question .eq. 'PARA_TEMP') then
        answerk = ' '
        field = object
        call jeveuo(field//'.VALE', 'L', jvale)
        call jelira(field//'.VALE', 'TYPE', cval=type)
        if (type(1:1) .eq. 'K') then
            call jelira(field//'.VALE', 'LONMAX', nb_zone)
            call jelira(field//'.VALE', 'LTYP', ltyp)
            do i_zone = 1, nb_zone
                if (ltyp .eq. 8) then
                    func_name = zk8(jvale+i_zone-1)
                else if (ltyp .eq. 24) then
                    func_name = zk24(jvale+i_zone-1) (1:19)
                else
                    ASSERT(.false.)
                end if
                if (func_name(1:8) .ne. ' ') then
                    call jeexin(func_name//'.PROL', iret)
                    if (iret .gt. 0) then
                        call jeveuo(func_name//'.PROL', 'L', vk24=v_prol)
                        call fonbpa(func_name, v_prol, func_type, max_para_name, nb_para, para_name)
                        do i_para = 1, nb_para
                            if (para_name(i_para) .eq. 'TEMP') then
                                answerk = 'OUI'
                                goto 999
                            end if
                        end do
                    end if
                end if
            end do
        end if
    else
        ierd = 1
    end if
!
999 continue
!
    answerk_ = answerk
    call jedema()
end subroutine
