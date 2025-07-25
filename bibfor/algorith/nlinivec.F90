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

subroutine nlinivec(sd_nl_, ip, lonvec, iocc, vi, &
                    vr, vc, vk8, vk16, vk24, address)
    implicit none
! Extract the value of a parameter in the temporary data structure for the
! non linearities of DYNA_VIBRA//TRAN/GENE
!
!  sd_nl_  [Obl]: Name of the nl data structure requested [K24]
!  ip      [Obl]: Index of the parameter relating to the vector to be created [I]
!  lonvec  [Obl]: Length of the vector to be created [I]
!  iocc    [Opt]: Index of the occurence in the case of a non global vector [I]
!  vi      [Opt]: Pointer to the created integer vector     [I]
!  vr      [Opt]: Pointer to the created real vector        [R]
!  vc      [Opt]: Pointer to the created complex vector     [C]
!  vk8     [Opt]: Pointer to the created K8 strings vector  [K8]
!  vk16    [Opt]: Pointer to the created K16 strings vector [K16]
!  vk24    [Opt]: Pointer to the created K24 strings vector [K24]
!
! Examples : call nlinivec(''&&OP29NL'', ROTR_DFK, 5, vi=indarch)
!
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "nldef.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/wkvect.h"
!
!   ====================================================================
!   = 0 =   Variable declarations and initialization
!   ====================================================================
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), intent(in) :: ip
    integer(kind=8), intent(in) :: lonvec
    integer(kind=8), optional, intent(in) :: iocc
    integer(kind=8), pointer, optional :: vi(:)
    real(kind=8), pointer, optional :: vr(:)
    complex(kind=8), pointer, optional :: vc(:)
    character(len=8), pointer, optional :: vk8(:)
    character(len=16), pointer, optional :: vk16(:)
    character(len=24), pointer, optional :: vk24(:)
    integer(kind=8), optional, intent(out) :: address
!
!   -0.2- Local variables
!   --- For strings copying
    character(len=8) :: sd_nl

!   --- For general usage
    integer(kind=8)           :: iret, add, i
    character(len=6)  :: k_iocc
    character(len=24) :: savename
!
#include "nlinc.h"
!
    savename = '                        '
    add = 0

!   Copying the input strings, in order to allow in-command truncated input
    sd_nl = sd_nl_

!   The parameter to be saved was not found in the predefined list
    if (ip .gt. _NL_NBPAR) then
        ASSERT(.false.)
    end if

    savename(1:8) = sd_nl
    if (present(iocc)) then
!       The parameter to be extracted is global but an occurence index was given
        ASSERT(parind(ip) .gt. 0)
        call codent(iocc, 'G', k_iocc)
        savename(9:15) = '.'//k_iocc(1:6)
    else
        ASSERT(parind(ip) .lt. 0)
    end if
    savename(16:24) = '.'//params(ip)

    call jeexin(savename, iret)
    if (iret .ne. 0) then
        ! write (*,*) "Careful, using nlinivec with a pre-existing vector"
        ! write (*,*) "The original <", params(ip) ,"> vector is deleted"
        call jelibe(savename)
        call jedetr(savename)
    end if

    if (partyp(ip) .eq. 'I') then
        if (present(vi)) then
            call wkvect(savename, 'V V I', lonvec, vi=vi)
            do i = 1, lonvec
                vi(i) = 0
            end do
        else
            call wkvect(savename, 'V V I', lonvec, add)
            do i = 1, lonvec
                zi(add+i-1) = 0
            end do
        end if
    else if (partyp(ip) .eq. 'R') then
        if (present(vr)) then
            call wkvect(savename, 'V V R', lonvec, vr=vr)
            do i = 1, lonvec
                vr(i) = 0.d0
            end do
        else
            call wkvect(savename, 'V V R', lonvec, add)
            do i = 1, lonvec
                zr(add+i-1) = 0.d0
            end do
        end if
    else if (partyp(ip) .eq. 'C') then
        if (present(vc)) then
            call wkvect(savename, 'V V C', lonvec, vc=vc)
            do i = 1, lonvec
                vc(i) = dcmplx(0.d0, 0.d0)
            end do
        else
            call wkvect(savename, 'V V C', lonvec, add)
            do i = 1, lonvec
                zc(add+i-1) = dcmplx(0.d0, 0.d0)
            end do
        end if
    else if (partyp(ip) .eq. 'K8') then
        if (present(vk8)) then
            call wkvect(savename, 'V V K8', lonvec, vk8=vk8)
            do i = 1, lonvec
                vk8(i) = ' '
            end do
        else
            call wkvect(savename, 'V V K8', lonvec, add)
            do i = 1, lonvec
                zk8(add+i-1) = ' '
            end do
        end if
    else if (partyp(ip) .eq. 'K16') then
        if (present(vk16)) then
            call wkvect(savename, 'V V K16', lonvec, vk16=vk16)
            do i = 1, lonvec
                vk16(i) = ' '
            end do
        else
            call wkvect(savename, 'V V K16', lonvec, add)
            do i = 1, lonvec
                zk16(add+i-1) = ' '
            end do
        end if
    else if (partyp(ip) .eq. 'K24') then
        if (present(vk24)) then
            call wkvect(savename, 'V V K24', lonvec, vk24=vk24)
            do i = 1, lonvec
                vk24(i) = ' '
            end do
        else
            call wkvect(savename, 'V V K24', lonvec, add)
            do i = 1, lonvec
                zk24(add+i-1) = ' '
            end do
        end if
    end if
    if (present(address)) then
        address = add
    end if

end subroutine
