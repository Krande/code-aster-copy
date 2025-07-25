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

subroutine ccchuc_chamno(field_in_s, field_out_s, nb_node, list_node, nb_cmp, type_comp, &
                         crit, nb_form, name_form, name_gd, nb_cmp_resu, work_out_val, &
                         nb_node_out, ichk)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ccchcf.h"
#include "asterfort/ccchcr.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeundf.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
!
    character(len=19), intent(in) :: field_in_s
    character(len=19), intent(in) :: field_out_s
    integer(kind=8), intent(in) :: nb_node
    character(len=24), intent(in) :: list_node
    integer(kind=8), intent(in) :: nb_cmp
    character(len=16), intent(in) :: type_comp
    character(len=16), intent(in) :: crit
    integer(kind=8), intent(in) :: nb_form
    character(len=8), intent(in) :: name_form(nb_form)
    character(len=8), intent(in) :: name_gd
    integer(kind=8), intent(in) :: nb_cmp_resu
    character(len=24), intent(in) :: work_out_val
    integer(kind=8), intent(out) :: ichk
    integer(kind=8), intent(out) :: nb_node_out
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP - CHAM_UTIL
!
! Compute CHAM_UTIL on <CHAM_NO>
!
! --------------------------------------------------------------------------------------------------
!
! In  field_in_s   : name of <CHAM_NO_S> input field FROM which extract values
! In  field_out_s  : name of <CHAM_NO_S> output field IN which compute values
! In  nb_node      : number of nodes
! In  list_node    : list of nodes
! In  nb_cmp       : number of components in input field
! In  type_comp    : type of computation (CRITERE or FORMULE)
! In  crit         : type of criterion
! In  nb_form      : number of formulas
! In  name_form    : names of formulas
! In  name_gd      : name of <GRANDEUR> of input field
! In  nb_cmp_resu  : number of components in output field
! In  work_out_val : working vector for output field (values)
! Out ichk         : 0 if OK
! Out nb_node_out  : number of nodes in output field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ino, icmp, nb_val_in, i_list_node, iel
    integer(kind=8) :: j_resu, ichk_node
    character(len=19) :: work_val, work_cmp
    integer(kind=8) :: j_val, j_cmp, j_nodein
    integer(kind=8) ::   jchsl
    integer(kind=8) :: jchrl
    character(len=8), pointer :: cnsc(:) => null()
    real(kind=8), pointer :: chrv(:) => null()
    real(kind=8), pointer :: chsv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ichk = -1
    nb_node_out = 0
    work_val = '&&CCCHUC_CHAMNO.VAL'
    work_cmp = '&&CCCHUC_CHAMNO.CMP'
!
! - Access to input field
!
    call jeveuo(field_in_s//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(field_in_s//'.CNSV', 'L', vr=chsv)
    call jeveuo(field_in_s//'.CNSL', 'L', jchsl)
!
! - Access to output field
!
    call jeveuo(field_out_s//'.CNSL', 'E', jchrl)
    call jeveuo(field_out_s//'.CNSV', 'E', vr=chrv)
!
! - Access to output working vector
!
    call jeveuo(work_out_val, 'E', j_resu)
!
! - Create working vectors
!
    call wkvect(work_val, 'V V R', nb_cmp, j_val)
    call wkvect(work_cmp, 'V V K8', nb_cmp, j_cmp)
!
    call jeexin(list_node, i_list_node)
    if (i_list_node .ne. 0) then
        call jeveuo(list_node, 'L', j_nodein)
    end if

    do iel = 1, nb_node

        if (i_list_node .ne. 0) then
            ino = zi(j_nodein-1+iel)
        else
            ino = iel
        end if
!
! ----- Undefine values
!
        call jeundf(work_val)
        call jeundf(work_cmp)
!
! ----- Set values
!
        nb_val_in = 0
        ichk_node = -1
        do icmp = 1, nb_cmp
            if (zl(jchsl-1+(ino-1)*nb_cmp+icmp)) then
                nb_val_in = nb_val_in+1
                zr(j_val-1+nb_val_in) = chsv((ino-1)*nb_cmp+icmp)
                zk8(j_cmp-1+nb_val_in) = cnsc(icmp)
            end if
        end do
!
! ----- Compute result
!
        if (type_comp .eq. 'CRITERE') then
            ASSERT(nb_cmp_resu .eq. 1)
            call ccchcr(crit, name_gd, nb_val_in, zr(j_val), zk8(j_cmp), &
                        nb_cmp_resu, zr(j_resu), ichk_node)
        elseif (type_comp .eq. 'FORMULE') then
            ASSERT(nb_cmp_resu .eq. nb_form)
            call ccchcf(name_form, nb_val_in, zr(j_val), zk8(j_cmp), nb_cmp_resu, &
                        zr(j_resu), ichk_node)
        else
            ASSERT(.false.)
        end if
!
! ----- Copy to output field
!
        if (ichk_node .eq. 0) then
            ichk = 0
            nb_node_out = nb_node_out+1
            do icmp = 1, nb_cmp_resu
                zl(jchrl-1+(ino-1)*nb_cmp_resu+icmp) = .true.
                chrv((ino-1)*nb_cmp_resu+icmp) = zr(j_resu-1+icmp)
            end do
        end if
    end do
!
    call jedetr(work_val)
    call jedetr(work_cmp)
!
    call jedema()
!
end subroutine
