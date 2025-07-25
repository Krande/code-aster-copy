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
subroutine nlsav(sd_nl_, ip, lonvec, iocc, kscal, &
                 iscal, rscal, cscal, kvect, ivect, &
                 rvect, cvect, buffer)
    implicit none
! Save a parameter in the temporary data structure for the nonlinear
! forces in DYNA_VIBRA//TRAN/GENE
!
!  sd_nl_  [Obl]: Name of the nl data structure to be saved into [K24]
!  ip      [Obl]: Index of the parameter to be saved [I]
!  lonvec  [Obl]: Length of the Vector to be saved, = 1 for scalars [I]
!  iocc    [Opt]: Index of the occurence, by default = 1 [I]kscal_
!  kscal   [Opt]: Value to be saved in the case of a character parameter [K24]
!  iscal   [Opt]: Value to be saved in the case of an integer parameter [I]
!  rscal   [Opt]: Value to be saved in the case of a float parameter [R8]
!  cscal   [Opt]: Value to be saved in the case of a complex parameter [C8]
!  kvect   [Opt]: Vector to be saved in the case of character parameters [K24]
!  ivect   [Opt]: Vector to be saved in the case of integer parameters [I]
!  rvect   [Opt]: Vector to be saved in the case of float parameters [R8]
!  cvect   [Opt]: Vector to be saved in the case of complex parameters [C8]
!
! 1 - First, we verify that the parameter name is valid
! 2 - Second, we save the given value/vector in the correct work vector
!     according to the sd_nl_ data structure's map
!
! Examples : call nlsav('&&OP0074','RESU_IN ',1, kscal=resuin)
!            call nlsav('&&OP0074','NOM_CMP ',1, iocc=2, kvect=('DX','DY'))
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "nldef.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
!   ====================================================================
!   = 0 =   Variable declarations and initialization
!   ====================================================================
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_nl_
    integer(kind=8), intent(in) :: ip
    integer(kind=8), intent(in) :: lonvec
    integer(kind=8), optional, intent(in) :: iocc
    character(len=*), optional, intent(in) :: kscal
    integer(kind=8), optional, intent(in) :: iscal
    real(kind=8), optional, intent(in) :: rscal
    complex(kind=8), optional, intent(in) :: cscal
    character(len=*), optional, intent(in) :: kvect(lonvec)
    integer(kind=8), optional, intent(in) :: ivect(lonvec)
    real(kind=8), optional, intent(in) :: rvect(lonvec)
    complex(kind=8), optional, intent(in) :: cvect(lonvec)
    integer(kind=8), pointer, optional :: buffer(:)
!
!   -0.2- Local variables
!   --- For strings copying
    character(len=8) :: sd_nl
    character(len=24) :: kscal_
    character(len=24) :: savejv
    character(len=24), pointer :: kvect_(:) => null()
!
!   --- For general usage
    aster_logical :: input_test
    integer(kind=8) :: i, jvect, jscal, iret
    integer(kind=8) :: dec, level, lvec, addr
    character(len=6) :: k_iocc
    blas_int :: b_incx, b_incy, b_n
!
#include "nlinc.h"
!
    savejv = '                        '
!
!   Copying the input strings, in order to allow in-command truncated input
    sd_nl = sd_nl_
!
    if (present(kscal)) kscal_ = kscal
    if (present(kvect)) then
        AS_ALLOCATE(vk24=kvect_, size=lonvec)
        do i = 1, lonvec
            kvect_(i) = kvect(i)
        end do
    end if
!
!   ====================================================================
!   = 1 = Validation of the input arguments, distinguishing global vars
!   ====================================================================
!
    input_test = UN_PARMI4(kscal, iscal, rscal, cscal) .or. UN_PARMI4(kvect, ivect, rvect, cvect)
!
    ASSERT(input_test)
!
    if (lonvec .gt. 1) then
        input_test = UN_PARMI4(kvect, ivect, rvect, cvect)
        ASSERT(input_test)
    end if
!
!   The parameter to be saved was not found in the predefined list
    if (ip .gt. _NL_NBPAR) then
        ASSERT(.false.)
    end if
!   Some verifications on parameter
    ASSERT(params(ip) .ne. 'XXXXXXXX')
    ASSERT(partyp(ip) .ne. 'XXX')
    ASSERT(parind(ip) .ne. 0)
!
    if (present(buffer)) then
!
        dec = 0
        if (present(iocc)) then
            level = size(buffer)/(2*_NL_NBPAR)
            if (iocc .le. level) then
                dec = (iocc-1)*2*_NL_NBPAR
            else
                goto 20
            end if
        end if
!
        addr = buffer(dec+ip)
        lvec = buffer(dec+_NL_NBPAR+ip)
!
        if ((addr .ne. 0) .and. (lonvec .eq. lvec)) then
            if (present(iscal)) then
                ASSERT(partyp(ip) .eq. 'I  ')
                zi(addr) = iscal
            else if (present(rscal)) then
                ASSERT(partyp(ip) .eq. 'R  ')
                zr(addr) = rscal
            else if (present(cscal)) then
                ASSERT(partyp(ip) .eq. 'C  ')
                zc(addr) = cscal
            else if (present(kscal)) then
                if (partyp(ip) .eq. 'K8 ') then
                    zk8(addr) = kscal
                else if (partyp(ip) .eq. 'K16') then
                    zk16(addr) = kscal
                else if (partyp(ip) .eq. 'K24') then
                    zk24(addr) = kscal
                else
                    ASSERT(.false.)
                end if
            else if (present(ivect)) then
                ASSERT(partyp(ip) .eq. 'I  ')
                do i = 1, lvec
                    zi(addr+i-1) = ivect(i)
                end do
            else if (present(rvect)) then
                ASSERT(partyp(ip) .eq. 'R  ')
                b_n = to_blas_int(lvec)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, rvect, b_incx, zr(addr), b_incy)
            else if (present(cvect)) then
                ASSERT(partyp(ip) .eq. 'C  ')
                b_n = to_blas_int(lvec)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, cvect, b_incx, zc(addr), b_incy)
            else if (present(kvect)) then
                if (partyp(ip) .eq. 'K8 ') then
                    do i = 1, lvec
                        zk8(addr+i-1) = kvect(i)
                    end do
                else if (partyp(ip) .eq. 'K16') then
                    do i = 1, lvec
                        zk16(addr+i-1) = kvect(i)
                    end do
                else if (partyp(ip) .eq. 'K24') then
                    do i = 1, lvec
                        zk24(addr+i-1) = kvect(i)
                    end do
                else
                    ASSERT(.false.)
                end if
            end if
            goto 99
        end if
    end if
!
20  continue
!
    savejv(1:8) = sd_nl
    if (present(iocc)) then
!       The parameter to be saved is global but an occurence index was given
        ASSERT(parind(ip) .gt. 0)
        call codent(iocc, 'G', k_iocc)
        savejv(9:15) = '.'//k_iocc(1:6)
    end if
    savejv(16:24) = '.'//params(ip)
!
!   ====================================================================
!   = 2 = Saving data
!   ====================================================================
!
!   --- Vectors
    if (abs(parind(ip)) .eq. 2) then
!
!       The parameter to be saved is a vector but no vector input was found
        input_test = UN_PARMI4(kvect, ivect, rvect, cvect)
        ASSERT(input_test)
        ASSERT(lonvec .ge. 1)
!
        call jeexin(savejv, iret)
        if (iret .gt. 0) then
            call jeveuo(savejv, 'E', jvect)
        else
            call wkvect(savejv, 'V V '//partyp(ip), lonvec, jvect)
        end if
!
        if (partyp(ip) .eq. 'K8') then
            do i = 1, lonvec
                zk8(jvect+i-1) = kvect_(i) (1:8)
            end do
        else if (partyp(ip) .eq. 'K16') then
            do i = 1, lonvec
                zk16(jvect+i-1) = kvect_(i) (1:16)
            end do
        else if (partyp(ip) .eq. 'K24') then
            do i = 1, lonvec
                zk24(jvect+i-1) = kvect_(i)
            end do
        else if (partyp(ip) .eq. 'R') then
            b_n = to_blas_int(lonvec)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rvect, b_incx, zr(jvect), b_incy)
        else if (partyp(ip) .eq. 'C') then
            b_n = to_blas_int(lonvec)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, cvect, b_incx, zc(jvect), b_incy)
        else if (partyp(ip) .eq. 'I') then
            do i = 1, lonvec
                zi(jvect+i-1) = ivect(i)
            end do
        end if
!
!   --- Scalars
    else if (abs(parind(ip)) .eq. 1) then
!
!       The parameter to be saved is a scalar but no scalar input was found
        input_test = UN_PARMI4(kscal, iscal, rscal, cscal)
        ASSERT(input_test)
!
        call jeexin(savejv, iret)
        if (iret .gt. 0) then
            call jeveuo(savejv, 'E', jscal)
        else
            call wkvect(savejv, 'V V '//partyp(ip), 1, jscal)
        end if
!
        if (partyp(ip) .eq. 'K8 ') then
            zk8(jscal) = kscal_(1:8)
        else if (partyp(ip) .eq. 'K16') then
            zk16(jscal) = kscal_(1:16)
        else if (partyp(ip) .eq. 'K24') then
            zk24(jscal) = kscal_
        else if (partyp(ip) .eq. 'R') then
            zr(jscal) = rscal
        else if (partyp(ip) .eq. 'C') then
            zc(jscal) = cscal
        else if (partyp(ip) .eq. 'I') then
            zi(jscal) = iscal
        end if
!
    end if
!
    if (present(kvect)) then
        AS_DEALLOCATE(vk24=kvect_)
    end if
!
99  continue
!
end subroutine
