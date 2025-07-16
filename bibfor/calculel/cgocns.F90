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

subroutine cgocns(cnoz, basez, cnsz, meshz, undf0_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/cheksd.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cmpcha.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: jacques.pellet at edf.fr
!
    character(len=*), intent(in) :: cnoz
    character(len=*), intent(in) :: cnsz, meshz
    character(len=*), intent(in) :: basez
    aster_logical, optional, intent(in) :: undf0_
!
! --------------------------------------------------------------------------------------------------
!
! BUT : TRANSFORMER UN CHAM_GEOM (CNOZ) EN CHAM_NO_S (CNSZ)
!
! --------------------------------------------------------------------------------------------------
!
!     ARGUMENTS:
! CNOZ    IN/JXIN  K19 : SD CHAM_GEOM A TRANSFORMER
! BASEZ   IN       K1  : BASE DE CREATION POUR CNSZ : G/V/L
! CNSZ    IN/JXOUT K19 : SD CHAM_NO_S A CREER
! undf0   IN           : flag to init field with zeros
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1) :: base
    character(len=3) :: tsca
    character(len=8) :: mesh, gran_name
    character(len=19) :: cno, cns
    integer(kind=8) :: nb_ec, idx_gd, nb_cmp_mx, nb_node, jvale, ierr, nb_node2
    integer(kind=8) :: i_node, nb_cmp, jcnsl, jcnsv
    integer(kind=8) :: ieq, i_cmp_field
    logical :: sdveri
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: cata_to_field(:) => null()
    integer(kind=8), pointer :: field_to_cata(:) => null()
    character(len=8), pointer :: cmp_name(:) => null()
    aster_logical :: undf0
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    cno = cnoz
    cns = cnsz
    base = basez
    mesh = meshz

    undf0 = ASTER_FALSE
    if (present(undf0_)) then
        undf0 = undf0_
    end if
!
! - Check ?
!
    sdveri = .false.
    if (sdveri) then
        call cheksd(cno, 'sd_cham_geom', ierr)
        ASSERT(ierr .eq. 0)
    end if
!
! - Old is erased
!
    call detrsd('CHAM_NO_S', cns)
!
! - Informations
!
    call dismoi('NOM_GD', cno, 'CHAM_GEOM', repk=gran_name)
    if (gran_name .ne. "GEOM_R") then
        ASSERT(ASTER_FALSE)
    end if
!
    call dismoi('NB_NO_MAILLA', cno, 'CHAM_GEOM', repi=nb_node)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node2)
    ASSERT(nb_node == nb_node2)
!
    call dismoi('NB_EC', gran_name, 'GRANDEUR', repi=nb_ec)
    call dismoi('NUM_GD', gran_name, 'GRANDEUR', repi=idx_gd)
    call dismoi('TYPE_SCA', gran_name, 'GRANDEUR', repk=tsca)

    call jeveuo(cno//'.VALE', 'L', jvale)
    call jeveuo(cno//'.DESC', 'L', vi=desc)
!
! - Create objects for global components (catalog) <=> local components (field)
!
    call cmpcha(cno, cmp_name, cata_to_field, field_to_cata, nb_cmp, &
                nb_cmp_mx)
!
! - Allocate CNS
!
    call cnscre(mesh, gran_name, nb_cmp, cmp_name, base, &
                cns, undf0)
!
! - Access to CNS
!
    call jeveuo(cns//'.CNSL', 'E', jcnsl)
    call jeveuo(cns//'.CNSV', 'E', jcnsv)
!
! - Set values in CNS
!
    do i_node = 1, nb_node
        do i_cmp_field = 1, nb_cmp
            zl(jcnsl-1+(i_node-1)*nb_cmp+i_cmp_field) = .true.
            ieq = (i_node-1)*nb_cmp+i_cmp_field
            if (tsca .eq. 'R') then
                zr(jcnsv-1+ieq) = zr(jvale-1+ieq)
            else if (tsca .eq. 'I') then
                zi(jcnsv-1+ieq) = zi(jvale-1+ieq)
            else if (tsca .eq. 'C') then
                zc(jcnsv-1+ieq) = zc(jvale-1+ieq)
            else if (tsca .eq. 'L') then
                zl(jcnsv-1+ieq) = zl(jvale-1+ieq)
            else if (tsca .eq. 'K8') then
                zk8(jcnsv-1+ieq) = zk8(jvale-1+ieq)
            else
                ASSERT(.false.)
            end if
        end do
    end do
!
! - Check
!
    if (sdveri) then
        call cheksd(cns, 'sd_cham_no_s', ierr)
        ASSERT(ierr .eq. 0)
    end if
!
    AS_DEALLOCATE(vi=cata_to_field)
    AS_DEALLOCATE(vi=field_to_cata)
    AS_DEALLOCATE(vk8=cmp_name)
    call jedema()
end subroutine
