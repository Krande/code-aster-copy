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

subroutine wkvectc(nom, carac, dim, pc)
    use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/assert.h"
#include "asterfort/jgetptc.h"
!
    character(len=*), intent(in) :: nom
    character(len=*), intent(in) :: carac
    integer(kind=8), intent(in) :: dim
!
! --------------------------------------------------------------------
! Creation d'un vecteur jeveux
! --------------------------------------------------------------------
! in  nom   : ch*24 : nom du vecteur jeveux
! in  carac : ch    : descrivion des caracteristiques pour jecreo
! in  dim   : is    : taille du vecteur
! out pc    : c_ptr : pointeur C vers le vecteur de valeurs
! --------------------------------------------------------------------
    integer(kind=8) :: jad
    character(len=8) :: ktyp
    type(c_ptr) :: pc
! --------------------------------------------------------------------
    call jecreo(nom, carac)
    call jeecra(nom, 'LONMAX', ival=dim)
    call jeecra(nom, 'LONUTI', ival=dim)
    call jeveuo(nom, 'E', jad)
!
!
!     -- cas : on veut l'adresse
!
!
!     -- cas : on veut un pointeur
    call jelira(nom, 'TYPELONG', cval=ktyp)
    if (ktyp .eq. 'L') then
        call jgetptc(jad, pc, vl=zl(1))
!
    else if (ktyp .eq. 'I') then
        call jgetptc(jad, pc, vi=zi(1))
!
    else if (ktyp .eq. 'S') then
        call jgetptc(jad, pc, vi4=zi4(1))
!
    else if (ktyp .eq. 'R') then
        call jgetptc(jad, pc, vr=zr(1))
!
    else if (ktyp .eq. 'C') then
        call jgetptc(jad, pc, vc=zc(1))
!
    else if (ktyp .eq. 'K8') then
        call jgetptc(jad, pc, vk8=zk8(1))
!
    else if (ktyp .eq. 'K16') then
        call jgetptc(jad, pc, vk16=zk16(1))
!
    else if (ktyp .eq. 'K24') then
        call jgetptc(jad, pc, vk24=zk24(1))
!
    else if (ktyp .eq. 'K32') then
        call jgetptc(jad, pc, vk32=zk32(1))
!
    else if (ktyp .eq. 'K80') then
        call jgetptc(jad, pc, vk80=zk80(1))
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
