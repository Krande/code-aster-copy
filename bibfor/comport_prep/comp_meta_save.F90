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
subroutine comp_meta_save(mesh, comporMeta, nbCmp, metaPrepBehaviour)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nocart.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: comporMeta
    integer(kind=8), intent(in) :: nbCmp
    type(META_PrepBehaviour), intent(in) :: metaPrepBehaviour
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Save informations in COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  comporMeta       : name of map for behaviour in metallurgy
! In  nbCmp            : number of components in <CARTE> COMPOR
! In  metaPrepBehaviour: datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: factorKeyword
    character(len=24), parameter :: list_elem_affe = '&&COMPMETASAVE.LIST'
    aster_logical :: l_affe_all
    integer(kind=8) :: nb_elem_affe
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword, nbPhase
    integer(kind=8) :: nbVari, numeComp
    character(len=16) :: metaType, metaLaw
    character(len=16), pointer :: comporMetaValv(:) => null()
    integer(kind=8), pointer :: v_elem_affe(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    factorKeyword = metaPrepBehaviour%factorKeyword
    nbFactorKeyword = metaPrepBehaviour%nbFactorKeyword

! - Access to COMPOR <CARTE>
    call jeveuo(comporMeta(1:19)//'.VALV', 'E', vk16=comporMetaValv)

! - Read list
    do iFactorKeyword = 1, nbFactorKeyword

! ----- Get options
        metaType = metaPrepBehaviour%paraBehaviour(iFactorKeyword)%metaType
        metaLaw = metaPrepBehaviour%paraBehaviour(iFactorKeyword)%metaLaw
        nbVari = metaPrepBehaviour%paraBehaviour(iFactorKeyword)%nbVari
        numeComp = metaPrepBehaviour%paraBehaviour(iFactorKeyword)%numeComp
        nbPhase = metaPrepBehaviour%paraBehaviour(iFactorKeyword)%nbPhase

! ----- Check define in Metallurgy_type.h
        if (metaLaw .eq. "WAECKEL") then
            ASSERT(nbVari .eq. NBVARIWAECKEL)
        elseif (metaLaw .eq. "JMA") then
            ASSERT(nbVari .eq. NBVARIJMA)
        end if

! ----- Set options
        comporMetaValv(1) = metaType
        write (comporMetaValv(2), '(I16)') nbVari
        comporMetaValv(3) = metaLaw
        write (comporMetaValv(4), '(I16)') numeComp
        write (comporMetaValv(5), '(I16)') nbPhase

! ----- Get list of elements where comportment is defined
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            list_elem_affe, l_affe_all, nb_elem_affe)

! ----- Affect in COMPOR <CARTE>
        if (l_affe_all) then
            call nocart(comporMeta, 1, nbCmp)
        else
            call jeveuo(list_elem_affe, 'L', vi=v_elem_affe)
            call nocart(comporMeta, 3, nbCmp, mode='NUM', nma=nb_elem_affe, &
                        limanu=v_elem_affe)
            call jedetr(list_elem_affe)
        end if
    end do
!
    call jedetr(comporMeta(1:19)//'.NCMP')
    call jedetr(comporMeta(1:19)//'.VALV')
!
end subroutine
