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
subroutine intget(sd_int_, ip, iocc, lonvec, savejv, &
                  iscal, rscal, cscal, kscal, ivect, &
                  rvect, cvect, kvect, vi, vr, &
                  vc, vk8, vk16, vk24, address, &
                  buffer)
    use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
    implicit none
! Extract the value of a parameter in the temporary data structure for an
! integration scheme in linear dynamics (DYNA_VIBRA)
!
!  sd_int_ [Obl]: Name of the integration data structure requested [K24]
!  ip      [Obl]: Index of the parameter to be extracted [I]
!  iocc    [Opt]: Index of the occurence, by default = 1 [I]
!  lonvec  [Opt]: Get the length of the vector [I]
!  savejv  [Opt]: Get the corresponding jeveux object name [K24]
!  kscal   [Opt]: Value to be extracted in the case of a character parameter [K24]
!  iscal   [Opt]: Value to be extracted in the case of an integer parameter [I]
!  rscal   [Opt]: Value to be extracted in the case of a float parameter [R8]
!  kvect   [Opt]: Vector to be extracted in the case of character parameters [K24]
!  ivect   [Opt]: Vector to be extracted in the case of integer parameters [I]
!  rvect   [Opt]: Vector to be extracted in the case of float parameters [R8]
!
! 1 - First, we verify that the parameter name is valid
! 2 - Second, we extract the value/vector from the correct work vector
!     according to the sd_int_ data structure's map
!
! Examples : call intget('&&OP0074','RESU_IN ',kscal=resu)
!
!            call intget('&&OP0074','NOM_CMP ',iocc=2, lonvec=nbcmp)
!            AS_ALLOCATE(vk8=cmp , size=nbcmp)
!            call intget('&&OP0074','NOM_CMP ',iocc=2, kvect=cmp)
!
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "intdef.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jgetptc.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
!
!
!   ====================================================================
!   = 0 =   Variable declarations and initialization
!   ====================================================================
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_int_
    integer(kind=8), intent(in) :: ip
    integer(kind=8), optional, intent(in) :: iocc
    character(len=24), optional, intent(out) :: savejv
    integer(kind=8), optional, intent(out) :: lonvec
    integer(kind=8), optional, intent(out) :: iscal
    real(kind=8), optional, intent(out) :: rscal
    complex(kind=8), optional, intent(out) :: cscal
    character(len=*), optional, intent(out) :: kscal
    integer(kind=8), optional, intent(out) :: ivect(*)
    real(kind=8), optional, intent(out) :: rvect(*)
    complex(kind=8), optional, intent(out) :: cvect(*)
    character(len=*), optional, intent(out) :: kvect(*)
!
    integer(kind=8), pointer, optional :: vi(:)
    real(kind=8), pointer, optional :: vr(:)
    complex(kind=8), pointer, optional :: vc(:)
    character(len=8), pointer, optional :: vk8(:)
    character(len=16), pointer, optional :: vk16(:)
    character(len=24), pointer, optional :: vk24(:)
!
    integer(kind=8), optional, intent(out) :: address
    integer(kind=8), pointer, optional :: buffer(:)
!
!   -0.2- Local variables
!   --- For strings copying
    character(len=8) :: sd_int
!
!   --- For general usage
    aster_logical :: output_test
    integer(kind=8) :: i, addr, jscal, lvec, level, dec
    character(len=6) :: k_iocc
    character(len=24) :: savename
    type(c_ptr) :: pc
    blas_int :: b_incx, b_incy, b_n
!
#include "intinc.h"
!
!
    savename = '                        '
    addr = 0
!
!   Copying the input strings, in order to allow in-command truncated input
    sd_int = sd_int_
!
!   ====================================================================
!   = 1 = Validation of the input arguments, distinguishing global vars
!   ====================================================================
!
    if ((.not. present(lonvec)) .and. (.not. present(savejv))) then
        output_test = UN_PARMI4(kscal, iscal, rscal, cscal) .or. &
                      UN_PARMI4(kvect, ivect, rvect, cvect) .or. UN_PARMI3(vk8, vk16, vk24) .or. &
                      UN_PARMI4(vi, vr, vc, address)
!
        ASSERT(output_test)
    end if
!
!   The parameter to be saved was not found in the predefined list
    if (ip .gt. _INT_NBPAR) then
        ASSERT(.false.)
    end if
!
    if (present(buffer) .and. (.not. present(savejv))) then
!
        dec = 0
        if (present(iocc)) then
            level = size(buffer)/(2*_INT_NBPAR)
            if (iocc .le. level) then
                dec = (iocc-1)*2*_INT_NBPAR
            else if (present(lonvec)) then
                lonvec = 0
                goto 99
            else
                goto 20
            end if
        end if
!
        addr = buffer(dec+ip)
        lvec = buffer(dec+_INT_NBPAR+ip)
!
        if (present(lonvec)) lonvec = lvec
        if (present(address)) address = addr
!
        if (addr .ne. 0) then
            if (present(iscal)) then
                iscal = zi(addr)
            else if (present(rscal)) then
                rscal = zr(addr)
            else if (present(cscal)) then
                cscal = zc(addr)
            else if (present(kscal)) then
                kscal = zk24(addr)
            else if (present(ivect)) then
                do i = 1, lvec
                    ivect(i) = zi(addr+i-1)
                end do
            else if (present(rvect)) then
                b_n = to_blas_int(lvec)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(addr), b_incx, rvect, b_incy)
            else if (present(cvect)) then
                b_n = to_blas_int(lvec)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, zc(addr), b_incx, cvect, b_incy)
            else if (present(kvect)) then
                do i = 1, lvec
                    kvect(i) = zk24(addr+i-1)
                end do
            else if (present(vi)) then
                call jgetptc(addr, pc, vi=zi(1))
                call c_f_pointer(pc, vi, [lvec])
            else if (present(vr)) then
                call jgetptc(addr, pc, vr=zr(1))
                call c_f_pointer(pc, vr, [lvec])
            else if (present(vc)) then
                call jgetptc(addr, pc, vc=zc(1))
                call c_f_pointer(pc, vc, [lvec])
            else if (present(vk8)) then
                call jgetptc(addr, pc, vk8=zk8(1))
                call c_f_pointer(pc, vk8, [lvec])
            else if (present(vk16)) then
                call jgetptc(addr, pc, vk16=zk16(1))
                call c_f_pointer(pc, vk16, [lvec])
            else if (present(vk24)) then
                call jgetptc(addr, pc, vk24=zk24(1))
                call c_f_pointer(pc, vk24, [lvec])
            end if
            if (present(address)) address = addr
            goto 99
        else if (present(lonvec)) then
            goto 99
        end if
    end if
!
20  continue
    savename(1:8) = sd_int
    if (present(iocc)) then
!       The parameter to be extracted is global but an occurence index was given
        ASSERT(parind(ip) .gt. 0)
        call codent(iocc, 'G', k_iocc)
        savename(9:15) = '.'//k_iocc(1:6)
    else
        ASSERT(parind(ip) .lt. 0)
    end if
    savename(16:24) = '.'//params(ip)
!
!   ====================================================================
!   = 2 = Extracting data
!   ====================================================================
!
!   --- Length of vectors
    if (present(savejv)) savejv = savename
!
    if (present(lonvec) .or. UN_PARMI4(kscal, iscal, rscal, cscal) .or. &
        UN_PARMI4(kvect, ivect, rvect, cvect) .or. UN_PARMI3(vk8, vk16, vk24) .or. &
        UN_PARMI3(vi, vr, vc)) then
        call jeexin(savename, lvec)
        if (lvec .le. 0) then
            if (present(lonvec)) then
                lonvec = 0
                goto 99
            end if
        else
            call jelira(savename, 'LONMAX', lvec)
        end if
    end if
!
    if (present(lonvec)) lonvec = lvec
!
    if (UN_PARMI4(kscal, iscal, rscal, cscal) .or. UN_PARMI4(kvect, ivect, rvect, cvect) &
        .or. UN_PARMI3(vk8, vk16, vk24) .or. UN_PARMI3(vi, vr, vc)) then
!   --- Vectors
        if (abs(parind(ip)) .eq. 2) then
!
            if (UN_PARMI4(kvect, ivect, rvect, cvect) .or. UN_PARMI3(vk8, vk16, vk24) .or. &
                UN_PARMI3(vi, vr, vc)) then
                if (partyp(ip) .eq. 'I') then
                    if (present(ivect)) then
                        call jeveuo(savename, 'L', addr)
                        do i = 1, lvec
                            ivect(i) = zi(addr+i-1)
                        end do
                    else
                        call jeveuo(savename, 'E', vi=vi)
                    end if
                else if (partyp(ip) .eq. 'R') then
                    if (present(rvect)) then
                        call jeveuo(savename, 'L', addr)
                        b_n = to_blas_int(lvec)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        call dcopy(b_n, zr(addr), b_incx, rvect, b_incy)
                    else
                        call jeveuo(savename, 'E', vr=vr)
                    end if
                else if (partyp(ip) .eq. 'C') then
                    if (present(cvect)) then
                        call jeveuo(savename, 'L', addr)
                        b_n = to_blas_int(lvec)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        call zcopy(b_n, zc(addr), b_incx, cvect, b_incy)
                    else
                        call jeveuo(savename, 'E', vc=vc)
                    end if
                else if (partyp(ip) .eq. 'K8') then
                    if (present(kvect)) then
                        call jeveuo(savename, 'L', addr)
                        do i = 1, lvec
                            kvect(i) = zk8(addr+i-1)
                        end do
                    else
                        call jeveuo(savename, 'E', vk8=vk8)
                    end if
                else if (partyp(ip) .eq. 'K16') then
                    if (present(kvect)) then
                        call jeveuo(savename, 'L', addr)
                        do i = 1, lvec
                            kvect(i) = zk16(addr+i-1)
                        end do
                    else
                        call jeveuo(savename, 'E', vk16=vk16)
                    end if
                else if (partyp(ip) .eq. 'K24') then
                    if (present(kvect)) then
                        call jeveuo(savename, 'L', addr)
                        do i = 1, lvec
                            kvect(i) = zk24(addr+i-1)
                        end do
                    else
                        call jeveuo(savename, 'E', vk24=vk24)
                    end if
                end if
            end if
!
!   --- Scalars
        else if (abs(parind(ip)) .eq. 1) then
!
!           The parameter to get is a scalar but no scalar output was found
            output_test = UN_PARMI4(kscal, iscal, rscal, cscal)
            ASSERT(output_test)
!
            call jeveuo(savename, 'L', jscal)
            if (partyp(ip) .eq. 'I') then
                iscal = zi(jscal)
            else if (partyp(ip) .eq. 'R') then
                rscal = zr(jscal)
            else if (partyp(ip) .eq. 'C') then
                cscal = zc(jscal)
            else if (partyp(ip) .eq. 'K8') then
                kscal = zk8(jscal)
            else if (partyp(ip) .eq. 'K16') then
                kscal = zk16(jscal)
            else if (partyp(ip) .eq. 'K24') then
                kscal = zk24(jscal)
            end if
!
        end if
    end if
!
!   jeveux memory address
    if (present(address)) then
        if (addr .ne. 0) then
            address = addr
        else
            call jeveuo(savename, 'E', address)
        end if
    end if
!
99  continue
!
end subroutine
