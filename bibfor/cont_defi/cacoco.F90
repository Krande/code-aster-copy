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

subroutine cacoco(sdcont, keywf, mesh)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cfdisi.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: sdcont
    character(len=16), intent(in) :: keywf
    character(len=8), intent(in) :: mesh
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get supplementary gap: shells
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  keywf            : factor keyword to read
! In  mesh             : name of mesh
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, noc
    integer(kind=8) :: nb_para_maxi, nb_cont_elem, nb_cont_zone, nb_slav_elem
    integer(kind=8) :: shell_ep_indx, shell_exc_indx, iad1
    integer(kind=8) :: elem_slav_indx, elem_slav_nume
    integer(kind=8) :: jdecme
    integer(kind=8) :: i_zone, i_slav_elem
    real(kind=8) :: shell_ep, shell_excent
    aster_logical :: l_dist_exist
    character(len=8) :: cara_elem, elem_slav_name
    character(len=19) :: cara_elem_s
    aster_logical :: l_dist_shell
    real(kind=8), pointer :: v_caraelem_cesv(:) => null()
    character(len=8), pointer :: v_caraelem_cesc(:) => null()
    integer(kind=8) :: j_caraelem_cesd, j_caraelem_cesl
    character(len=24) :: sdcont_defi
    character(len=24) :: sdcont_mailco
    integer(kind=8), pointer :: v_sdcont_mailco(:) => null()
    character(len=24) :: sdcont_jeucoq
    real(kind=8), pointer :: v_sdcont_jeucoq(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    cara_elem_s = '&&CAPOCO.CARGEOPO'
!
! - Datastructure for contact definition
!
    sdcont_defi = sdcont(1:8)//'.CONTACT'
    sdcont_mailco = sdcont_defi(1:16)//'.MAILCO'
    sdcont_jeucoq = sdcont_defi(1:16)//'.JEUCOQ'
    call jeveuo(sdcont_mailco, 'L', vi=v_sdcont_mailco)
!
! - Parameters
!
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    nb_cont_elem = cfdisi(sdcont_defi, 'NMACO')
!
! - Create shell gap datastructure
!
    call wkvect(sdcont_jeucoq, 'G V R', nb_cont_elem, vr=v_sdcont_jeucoq)
!
! - Get elementary characteristics datastructure
!
    l_dist_exist = .false.
    do i_zone = 1, nb_cont_zone
        l_dist_shell = mminfl(sdcont_defi, 'DIST_COQUE', i_zone)
        if (l_dist_shell) then
            l_dist_exist = .true.
            call getvid(keywf, 'CARA_ELEM', iocc=i_zone, scal=cara_elem, nbret=noc)
            ASSERT(noc .ne. 0)
        end if
    end do
!
    if (.not. l_dist_exist) then
        goto 999
    end if
!
! - Access to elementary characteristics
!
    call carces(cara_elem//'.CARCOQUE', 'ELEM', ' ', 'V', cara_elem_s, &
                'A', iret)
    call jeveuo(cara_elem_s//'.CESC', 'L', vk8=v_caraelem_cesc)
    call jeveuo(cara_elem_s//'.CESD', 'L', j_caraelem_cesd)
    call jeveuo(cara_elem_s//'.CESL', 'L', j_caraelem_cesl)
    call jeveuo(cara_elem_s//'.CESV', 'L', vr=v_caraelem_cesv)
!
! - Get index for storing shell parameters
!
    nb_para_maxi = zi(j_caraelem_cesd-1+2)
    shell_ep_indx = indik8(v_caraelem_cesc, 'EP      ', 1, nb_para_maxi)
    shell_exc_indx = indik8(v_caraelem_cesc, 'EXCENT  ', 1, nb_para_maxi)
!
! - Loop on contact zones
!
    do i_zone = 1, nb_cont_zone
        l_dist_shell = mminfl(sdcont_defi, 'DIST_COQUE', i_zone)
        if (l_dist_shell) then
            nb_slav_elem = mminfi(sdcont_defi, 'NBMAE', i_zone)
            jdecme = mminfi(sdcont_defi, 'JDECME', i_zone)
            do i_slav_elem = 1, nb_slav_elem
!
! ------------- Current element
!
                elem_slav_indx = jdecme+i_slav_elem
                elem_slav_nume = v_sdcont_mailco(elem_slav_indx)
                elem_slav_name = int_to_char8(elem_slav_nume)
!
! ------------- Get thickness
!
                call cesexi('C', j_caraelem_cesd, j_caraelem_cesl, elem_slav_nume, 1, &
                            1, shell_ep_indx, iad1)
                if (iad1 .gt. 0) then
                    shell_ep = v_caraelem_cesv(iad1)
                else
                    call utmess('F', 'CONTACT3_39', sk=elem_slav_name)
                end if
!
! ------------- Get excentricity
!
                call cesexi('C', j_caraelem_cesd, j_caraelem_cesl, elem_slav_nume, 1, &
                            1, shell_exc_indx, iad1)
                if (iad1 .gt. 0) then
                    shell_excent = v_caraelem_cesv(iad1)
                    if (shell_excent .ge. r8prem()) then
                        call utmess('F', 'CONTACT3_40', sk=elem_slav_name)
                    end if
                else
                    call utmess('F', 'CONTACT3_41', sk=elem_slav_name)
                end if
!
! ------------- Save
!
                v_sdcont_jeucoq(elem_slav_indx) = 0.5d0*shell_ep
            end do
        end if
    end do
!
999 continue
!
    call detrsd('CHAM_ELEM_S', cara_elem_s)
!
end subroutine
