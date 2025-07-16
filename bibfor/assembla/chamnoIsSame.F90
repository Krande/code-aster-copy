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
subroutine chamnoIsSame(chamno1_, chamno2_, ier)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/idenob.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
!
    character(len=*), intent(in) :: chamno1_, chamno2_
    integer(kind=8), intent(out) :: ier
!
! --------------------------------------------------------------------------------------------------
!
! CHAMNO - Is the same ?
!
! --------------------------------------------------------------------------------------------------
!
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_same
    integer(kind=8) :: nb_refe1, nb_refe2, nb_ligr1, nb_ligr2, nbno
    character(len=19) :: chamno1, chamno2, numeq1, numeq2, ligrel
    character(len=24) :: refe1, refe2
    character(len=24), pointer :: v_refe1(:) => null()
    character(len=24), pointer :: v_refe2(:) => null()
    integer(kind=8), pointer :: v_nbno(:) => null()
    integer(kind=8) :: iexi
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    ier = 0
    chamno1 = chamno1_
    chamno2 = chamno2_
    nbno = 0
!
! - For REFE objects
!
    refe1 = chamno1//'.REFE'
    refe2 = chamno2//'.REFE'
    call jelira(refe1, 'LONMAX', nb_refe1)
    call jelira(refe2, 'LONMAX', nb_refe2)
    if (nb_refe1 .ne. nb_refe2) then
        ier = ier+abs(nb_refe1-nb_refe2)
    end if
    call jeveuo(refe1, 'L', vk24=v_refe1)
    call jeveuo(refe2, 'L', vk24=v_refe2)
    if (v_refe1(1) .ne. v_refe2(1)) then
        ier = ier+1
    end if
!
! - For NUMEEQUA
!
    numeq1 = v_refe1(2) (1:19)
    numeq2 = v_refe2(2) (1:19)
    l_same = idenob(numeq1//'.DEEQ', numeq2//'.DEEQ')
    if (.not. l_same) then
        ier = ier+1
    end if
    l_same = idenob(numeq1//'.NUEQ', numeq2//'.NUEQ')
    if (.not. l_same) then
        ier = ier+1
    end if
    l_same = idenob(numeq1//'.LILI', numeq2//'.LILI')
    if (l_same) then
        l_same = idenob(numeq1//'.PRNO', numeq2//'.PRNO')
        if (.not. l_same) then
            ier = ier+1
        end if
    else
        nbno = 0
        call jelira(numeq1//'.LILI', 'NOMMAX', nb_ligr1)
        call jelira(numeq2//'.LILI', 'NOMMAX', nb_ligr2)
        if (nb_ligr1 .eq. 2 .and. nb_ligr2 .eq. 1) then
            call jenuno(jexnum(numeq1//'.LILI', 2), ligrel)
            call jeexin(ligrel(1:19)//'.NBNO', iexi)
            if (iexi .gt. 0) then
                call jeveuo(ligrel(1:19)//'.NBNO', 'L', vi=v_nbno)
                nbno = v_nbno(1)
                if (nbno .ne. 0) then
                    ier = ier+1
                end if
            end if
        elseif (nb_ligr2 .eq. 2 .and. nb_ligr1 .eq. 1) then
            call jenuno(jexnum(numeq2//'.LILI', 2), ligrel)
            call jeexin(ligrel(1:19)//'.NBNO', iexi)
            if (iexi .gt. 0) then
                call jeveuo(ligrel(1:19)//'.NBNO', 'L', vi=v_nbno)
                nbno = v_nbno(1)
                if (nbno .ne. 0) then
                    ier = ier+1
                end if
            end if
        else
            ier = ier+1
        end if
    end if
!
    call jedema()
end subroutine
