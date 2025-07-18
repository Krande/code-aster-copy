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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine rescmp(cnfnod, cnequi, &
                  r_comp_vale, r_comp_name, r_comp_indx)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    character(len=19), intent(in) :: cnfnod, cnequi
    real(kind=8), intent(out) :: r_comp_vale
    character(len=8), intent(out) :: r_comp_name
    integer(kind=8), intent(out) :: r_comp_indx
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Compute RESI_COMP_RELA
!
! --------------------------------------------------------------------------------------------------
!
! In  cnfnod           : nodal field for internal force
! In  cnequi           : nodal field for out-of-balance force
! Out r_comp_indx      : number of node where RESI_COMP_RELA is maximum
! Out r_comp_vale      : value of RESI_COMP_RELA
! Out r_comp_name      : name of component where RESI_COMP_RELA is maximum
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nddmax = 6
    character(len=8) :: nomddl(nddmax)
    real(kind=8) :: maxddf(nddmax), maxddr(nddmax)
    integer(kind=8) :: numnod(nddmax)
    character(len=3) :: tsca
    integer(kind=8) :: cmpmax
    character(len=19) :: cnfnod_s, cnequi_s
    integer(kind=8) :: i, k
    real(kind=8) :: resim, fonam, res
    integer(kind=8) :: jcnsl
    integer(kind=8) :: licmpu(999)
    integer(kind=8) :: nb_cmp, nb_node, i_cmp, i_node, nbcmpu
    character(len=8) :: nomgd
    real(kind=8), parameter :: epsi = 1.d-50
    real(kind=8), pointer :: v_cnequi_cnsv(:) => null()
    real(kind=8), pointer :: v_cnfnod_cnsv(:) => null()
    character(len=8), pointer :: v_cnfnod_cnsc(:) => null()
    character(len=8), pointer :: v_cnfnod_cnsk(:) => null()
    integer(kind=8), pointer :: v_cnfnod_cnsd(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Convert to "simple" fields
!
    cnfnod_s = '&&RESCMP.CNFNOD_S'
    cnequi_s = '&&RESCMP.CNEQUI_S'
    call cnocns(cnfnod, 'V', cnfnod_s)
    call jeveuo(cnfnod_s(1:19)//'.CNSV', 'L', vr=v_cnfnod_cnsv)
    call jeveuo(cnfnod_s(1:19)//'.CNSD', 'L', vi=v_cnfnod_cnsd)
    call jeveuo(cnfnod_s(1:19)//'.CNSL', 'L', jcnsl)
    call jeveuo(cnfnod_s(1:19)//'.CNSC', 'L', vk8=v_cnfnod_cnsc)
    call jeveuo(cnfnod_s(1:19)//'.CNSK', 'L', vk8=v_cnfnod_cnsk)
    call cnocns(cnequi, 'V', cnequi_s)
    call jeveuo(cnequi_s(1:19)//'.CNSV', 'L', vr=v_cnequi_cnsv)
!
    nb_cmp = v_cnfnod_cnsd(2)
    nb_node = v_cnfnod_cnsd(1)

!
! --- NB DE CMP DANS LE CHAMP
!
    nbcmpu = 0
    do i_cmp = 1, nb_cmp
        do i_node = 1, nb_node
            if (zl(jcnsl-1+(i_node-1)*nb_cmp+i_cmp)) goto 20
        end do
        goto 30
20      continue
        nbcmpu = nbcmpu+1
        ASSERT(nbcmpu .lt. 999)
        licmpu(nbcmpu) = i_cmp
30      continue
    end do
!
! - Some checks
!
    nomgd = v_cnfnod_cnsk(2)
    ASSERT(nbcmpu .le. nddmax)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R')
!
    do i_cmp = 1, nbcmpu
        nomddl(i_cmp) = v_cnfnod_cnsc(licmpu(i_cmp))
        maxddf(i_cmp) = 0.d0
        maxddr(i_cmp) = 0.d0
        numnod(i_cmp) = 0
    end do
!
!
    do i_node = 1, nb_node
        do i_cmp = 1, nbcmpu
            k = licmpu(i_cmp)
            if (zl(jcnsl-1+(i_node-1)*nb_cmp+k)) then
                i = nb_cmp*(i_node-1)+k
                resim = abs(v_cnequi_cnsv(i))
                fonam = abs(v_cnfnod_cnsv(i))
                if (resim .gt. maxddr(i_cmp)) then
                    maxddr(i_cmp) = resim
                    numnod(i_cmp) = i_node
                end if
                maxddf(i_cmp) = max(fonam, maxddf(i_cmp))
            end if
        end do
    end do
    r_comp_vale = 0.d0
!
!
    do i_cmp = 1, nbcmpu
        if (maxddf(i_cmp) .gt. 0.d0) then
            res = maxddr(i_cmp)/maxddf(i_cmp)
        else
            res = -1
        end if
        if (res .gt. r_comp_vale) then
            r_comp_vale = res
            cmpmax = i_cmp
        end if
    end do
!
!  POUR INFO SI BESOIN NUMDDL  : NUMERO DU DDL PENALISANT
!    NUMDDL   = NUMN(CMPMAX)
!
    if (r_comp_vale .lt. epsi) then
        r_comp_vale = -1.d0
        r_comp_indx = 0
        r_comp_name = '   '
    else
        r_comp_indx = numnod(cmpmax)
        r_comp_name = nomddl(cmpmax)
    end if
!
    call detrsd('CHAM_NO_S', cnequi_s)
    call detrsd('CHAM_NO_S', cnfnod_s)
!
    call jedema()
!
end subroutine
