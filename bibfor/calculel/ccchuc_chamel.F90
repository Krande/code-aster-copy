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

subroutine ccchuc_chamel(field_in_s, field_out_s, nb_elem, list_elem, nb_cmp, type_comp, &
                         crit, nb_form, name_form, name_gd, nb_cmp_resu, &
                         work_out_val, work_out_ele, nb_elem_out, ichk)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/ccchcf.h"
#include "asterfort/ccchcr.h"
#include "asterfort/cesexi.h"
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
    integer(kind=8), intent(in) :: nb_elem
    character(len=24), intent(in) :: list_elem
    integer(kind=8), intent(in) :: nb_cmp
    character(len=16), intent(in) :: type_comp
    character(len=16), intent(in) :: crit
    integer(kind=8), intent(in) :: nb_form
    character(len=8), intent(in) :: name_form(nb_form)
    character(len=8), intent(in) :: name_gd
    integer(kind=8), intent(in) :: nb_cmp_resu
    character(len=24), intent(in) :: work_out_val
    character(len=24), intent(in) :: work_out_ele
    integer(kind=8), intent(out) :: ichk
    integer(kind=8), intent(out) :: nb_elem_out
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP - CHAM_UTIL
!
! Compute CHAM_UTIL on <CHAM_ELEM>
!
! --------------------------------------------------------------------------------------------------
!
! In  field_in_s   : name of <CHAM_ELEM_S> input field FROM which extract values
! In  field_out_s  : name of <CHAM_ELEM_S> output field IN which compute values
! In  nb_elem      : number of elements
! In  list_elem    : list of elements
! In  nb_cmp       : number of components in input field
! In  type_comp    : type of computation (CRITERE or FORMULE)
! In  crit         : type of criterion
! In  nb_form      : number of formulas
! In  name_form    : names of formulas
! In  name_gd      : name of <GRANDEUR> of input field
! In  nb_cmp_resu  : number of components in output field
! In  work_out_val : working vector for output field (values)
! In  work_out_ele : working vector for output field (elements)
! Out ichk         : 0 if OK
! Out nb_elem_out  : number of elements in output field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ima, ipt, isp, icmp, iad, nb_val_in
    integer(kind=8) :: nbpt, nbsp, nbcmp, ichk_elem, ichk_sp
    integer(kind=8) :: j_resu, j_elem, i_list_elem, j_elemin, iel
    character(len=24) :: work_val, work_cmp
    integer(kind=8) :: j_val, j_cmp
    integer(kind=8) ::   jchsl, jchsd
    integer(kind=8) :: jchrd, jchrl
    real(kind=8), pointer :: chrv(:) => null()
    real(kind=8), pointer :: chsv(:) => null()
    character(len=8), pointer :: cesc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    ichk = -1
    nb_elem_out = 0
    work_val = '&&CCCHUC_CHAMEL.VAL'
    work_cmp = '&&CCCHUC_CHAMEL.CMP'
!
! - Access to input field
!
    call jeveuo(field_in_s//'.CESL', 'L', jchsl)
    call jeveuo(field_in_s//'.CESV', 'L', vr=chsv)
    call jeveuo(field_in_s//'.CESC', 'L', vk8=cesc)
    call jeveuo(field_in_s//'.CESD', 'L', jchsd)
!
! - Access to output field
!
    call jeveuo(field_out_s//'.CESL', 'E', jchrl)
    call jeveuo(field_out_s//'.CESD', 'E', jchrd)
    call jeveuo(field_out_s//'.CESV', 'E', vr=chrv)
!
! - Access to output working vector
!
    call jeveuo(work_out_val, 'E', j_resu)
!
! - Access to work vector for element in out field
!
    call jeveuo(work_out_ele, 'E', j_elem)
!
! - Create working vectors
!
    call wkvect(work_val, 'V V R', nb_cmp, j_val)
    call wkvect(work_cmp, 'V V K8', nb_cmp, j_cmp)
!
    call jeexin(list_elem, i_list_elem)
    if (i_list_elem .ne. 0) then
        call jeveuo(list_elem, 'L', j_elemin)
    end if

    do iel = 1, nb_elem

        if (i_list_elem .ne. 0) then
            ima = zi(j_elemin-1+iel)
        else
            ima = iel
        end if

        ichk_elem = -1
        nbpt = zi(jchsd-1+5+4*(ima-1)+1)
        nbsp = zi(jchsd-1+5+4*(ima-1)+2)
        nbcmp = zi(jchsd-1+5+4*(ima-1)+3)
        do ipt = 1, nbpt
            do isp = 1, nbsp
!
! ------------- Undefine values
!
                call jeundf(work_val)
                call jeundf(work_cmp)
!
! ------------- Set values
!
                nb_val_in = 0
                do icmp = 1, nbcmp
                    call cesexi('S', jchsd, jchsl, ima, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        nb_val_in = nb_val_in+1
                        zr(j_val-1+nb_val_in) = chsv(iad)
                        zk8(j_cmp-1+nb_val_in) = cesc(icmp)
                    end if
                end do
!
! ------------- Compute result
!
                if (type_comp .eq. 'CRITERE') then
                    ASSERT(nb_cmp_resu .eq. 1)
                    call ccchcr(crit, name_gd, nb_val_in, zr(j_val), zk8(j_cmp), &
                                nb_cmp_resu, zr(j_resu), ichk_sp)
                else if (type_comp .eq. 'FORMULE') then
                    ASSERT(nb_cmp_resu .eq. nb_form)
                    call ccchcf(name_form, nb_val_in, zr(j_val), zk8(j_cmp), nb_cmp_resu, &
                                zr(j_resu), ichk_sp)
                else
                    ASSERT(.false.)
                end if
!
! ------------- Copy to output field
!
                if (ichk_sp .eq. 0) then
                    ichk_elem = 0
                    do icmp = 1, nb_cmp_resu
                        call cesexi('S', jchrd, jchrl, ima, ipt, &
                                    isp, icmp, iad)
                        iad = -iad
                        zl(jchrl-1+iad) = .true.
                        chrv(iad) = zr(j_resu-1+icmp)
                    end do
                end if
            end do
        end do
!
! ----- Add element computed
!
        if (ichk_elem .eq. 0) then
            ichk = 0
            nb_elem_out = nb_elem_out+1
            zi(j_elem-1+nb_elem_out) = ima
        end if
    end do
!
    call jedetr(work_val)
    call jedetr(work_cmp)
!
    call jedema()
!
end subroutine
