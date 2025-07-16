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

subroutine comp_comp_save(mesh, compor, v_info_valk, v_info_vali)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: compor
    character(len=16), intent(in) :: v_info_valk(:)
    integer(kind=8), intent(in) :: v_info_vali(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (AFFE_MATERIAU)
!
! Save informations in COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  compor           : name of <CARTE> COMPOR
! In  v_info_valk      : comportment informations (character)
! In  v_info_vali      : comportment informations (integer)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: list_elem_affe = '&&COMPCOMPSAVE.LIST'
    aster_logical :: l_affe_all
    integer(kind=8) :: nb_elem_affe
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: rela_comp, defo_comp, type_comp, type_cpla, mult_comp, kit_comp(4)
    character(len=16) :: post_iter
    integer(kind=8) :: nb_vari, nume_comp(4), nb_vari_exte, unit_comp, nb_vari_comp(4)
    character(len=16), parameter :: factorKeyword = 'AFFE_COMPOR'
    character(len=16), pointer :: comporValv(:) => null()
    integer(kind=8), pointer :: v_elem_affe(:) => null()
    integer(kind=8), parameter :: nbCmp = COMPOR_SIZE
!
! --------------------------------------------------------------------------------------------------
!
    call getfac(factorKeyword, nbFactorKeyword)

! - Access to map
    call jeveuo(compor//'.VALV', 'E', vk16=comporValv)

! - Read list
    do iFactorKeyword = 1, nbFactorKeyword
        nume_comp = 0
        nb_vari_comp = 0
        nb_vari_exte = v_info_vali(4*(iFactorKeyword-1)+1)
        unit_comp = v_info_vali(4*(iFactorKeyword-1)+2)
        nb_vari = v_info_vali(4*(iFactorKeyword-1)+3)
        nume_comp(1) = v_info_vali(4*(iFactorKeyword-1)+4)
        rela_comp = v_info_valk(16*(iFactorKeyword-1)+1)
        defo_comp = v_info_valk(16*(iFactorKeyword-1)+2)
        type_comp = v_info_valk(16*(iFactorKeyword-1)+3)
        type_cpla = v_info_valk(16*(iFactorKeyword-1)+4)
        kit_comp(1) = v_info_valk(16*(iFactorKeyword-1)+5)
        kit_comp(2) = v_info_valk(16*(iFactorKeyword-1)+6)
        kit_comp(3) = v_info_valk(16*(iFactorKeyword-1)+7)
        kit_comp(4) = v_info_valk(16*(iFactorKeyword-1)+8)
        mult_comp = v_info_valk(16*(iFactorKeyword-1)+14)
        post_iter = v_info_valk(16*(iFactorKeyword-1)+16)

! ----- Set options in COMPOR <CARTE>
        comporValv(RELA_NAME) = rela_comp
        write (comporValv(NVAR), '(I16)') nb_vari
        comporValv(DEFO) = defo_comp
        comporValv(INCRELAS) = type_comp
        comporValv(PLANESTRESS) = type_cpla
        write (comporValv(NUME), '(I16)') nume_comp(1)
        comporValv(MULTCOMP) = mult_comp
        comporValv(KIT1_NAME) = kit_comp(1)
        comporValv(KIT2_NAME) = kit_comp(2)
        comporValv(KIT3_NAME) = kit_comp(3)
        comporValv(KIT4_NAME) = kit_comp(4)
        comporValv(POSTITER) = post_iter
        write (comporValv(KIT1_NUME), '(I16)') nume_comp(2)
        write (comporValv(KIT2_NUME), '(I16)') nume_comp(3)
        write (comporValv(KIT3_NUME), '(I16)') 0
        write (comporValv(KIT4_NUME), '(I16)') 0
        write (comporValv(KIT1_NVAR), '(I16)') nb_vari_comp(1)
        write (comporValv(KIT2_NVAR), '(I16)') nb_vari_comp(2)
        write (comporValv(KIT3_NVAR), '(I16)') nb_vari_comp(3)
        write (comporValv(KIT4_NVAR), '(I16)') nb_vari_comp(4)

! ----- Get list of elements where comportment is defined
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            list_elem_affe, l_affe_all, nb_elem_affe)
!
! ----- Affect in COMPOR <CARTE>
!
        if (l_affe_all) then
            call nocart(compor, 1, nbCmp)
        else
            call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
            call nocart(compor, 3, nbCmp, mode='NUM', nma=nb_elem_affe, &
                        limanu=v_elem_affe)
            call jedetr(list_elem_affe)
        end if
    end do
!
    call jedetr(compor//'.NCMP')
    call jedetr(compor//'.VALV')
!
end subroutine
