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
subroutine pgpsav(sd_pgp, param, lonvec, iobs, kscal, &
                  iscal, rscal, cscal, kvect, ivect, &
                  rvect, cvect, savejv)
    implicit none
! Save a parameter in the temporary data structure for the command
! POST_GENE_PHYS
!
!  sd_pgp [Obl]: Name of the pgp data structure to be saved into [K24]
!  param  [Obl]: Name of the parameter to be saved [K24]
!  lonvec [Obl]: Length of the Vector to be saved, = 1 for scalars [I]
!  iobs   [Opt]: Index of the observation, by default = 1 [I]kscal_
!  kscal  [Opt]: Value to be saved in the case of a character parameter [K24]
!  iscal  [Opt]: Value to be saved in the case of an integer parameter [I]
!  rscal  [Opt]: Value to be saved in the case of a float parameter [R8]
!  cscal  [Opt]: Value to be saved in the case of a complex parameter [C8]
!  kvect  [Opt]: Vector to be saved in the case of character parameters [K24]
!  ivect  [Opt]: Vector to be saved in the case of integer parameters [I]
!  rvect  [Opt]: Vector to be saved in the case of float parameters [R8]
!  cvect  [Opt]: Vector to be saved in the case of complex parameters [C8]
!
! 1 - First, we verify that the parameter name is valid
! 2 - Second, we save the given value/vector in the correct work vector
!     according to the sd_pgp data structure's map
!
! Examples : call pgpsav('&&OP0058','RESU_IN ',1, kscal=resuin)
!            call pgpsav('&&OP0058','NOM_CMP ',1, iobs=2, kvect=('DX','DY'))
!
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
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
    character(len=*), intent(in) :: sd_pgp
    character(len=*), intent(in) :: param
    integer(kind=8), intent(in) :: lonvec
    integer(kind=8), optional, intent(in) :: iobs
    character(len=*), optional, intent(in) :: kscal
    integer(kind=8), optional, intent(in) :: iscal
    real(kind=8), optional, intent(in) :: rscal
    complex(kind=8), optional, intent(in) :: cscal
    character(len=*), optional, intent(in) :: kvect(lonvec)
    integer(kind=8), optional, intent(in) :: ivect(lonvec)
    real(kind=8), optional, intent(in) :: rvect(lonvec)
    complex(kind=8), optional, intent(in) :: cvect(lonvec)
    character(len=24), optional, intent(out) :: savejv
!
!   -0.2- Local variables
!   --- For strings copying
    character(len=8) :: sd_pgp_
    character(len=8) :: param_
    character(len=24) :: kscal_
    character(len=24) :: savejv_
    character(len=24), pointer :: kvect_(:) => null()
!
!
!   --- For general usage
    aster_logical :: input_test
    integer(kind=8) :: nbparams
    parameter(nbparams=24)
    integer(kind=8) :: parind(nbparams), ip, i, jvect, jscal
    character(len=3) :: partyp(nbparams)
    character(len=6) :: k_iobs
    character(len=8) :: params(nbparams)
    blas_int :: b_incx, b_incy, b_n
!
!   -0.3- Initialization
    data params/'RESU_OUT', 'RESU_IN ', 'TYP_RESU', 'BASE    ', 'MODELE  ', &
        'MAILLAGE', 'NB_OBSER', 'NOM_CHAM', 'TYP_CHAM', 'NOM_CMP ', &
        'TYP_SCAL', 'NUM_NOEU', 'NUM_MAIL', 'NUM_ORDR', 'DISC    ', &
        'ADD_CORR', 'ACC_MO_A', 'ACC_DIR ', 'VEC_PR_R', &
        'VEC_PR_C', 'REF_SUP1', 'REF_SUP2', 'REF_COMP', 'REF_INDI'/
!
    data partyp/'K24', 'K24', 'K24', 'K24', 'K24', &
        'K24', 'I', 'K24', 'K24', 'K24', &
        'K24', 'I', 'I', 'I', 'R8', &
        'I', 'K24', 'R8', 'R8', &
        'C8', 'K24', 'K24', 'K24', 'I'/
!
!   parind = -2 : vector global          ; = -1 : scalar global ;
!          =  2 : vector per observation ; =  1 : scalar per observation
    data parind/-1, -1, -1, -1, -1, &
        -1, -1, 1, 1, 2, &
        1, 2, 2, 2, 2, &
        1, 2, 2, 2, &
        2, 2, 2, 2, 2/
!
    call jemarq()
    savejv_ = '                        '
!
!   Copying the input strings, in order to allow in-command truncated input
    sd_pgp_ = sd_pgp
    param_ = param
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
        ASSERT(UN_PARMI4(kvect, ivect, rvect, cvect))
    end if
!
    do ip = 1, nbparams
        if (params(ip) .eq. param_) goto 10
    end do
10  continue
!
!   The parameter to be saved was not found in the predefined list
    if (ip .eq. nbparams+1) then
        ASSERT(.false.)
    end if
!
    savejv_(1:8) = sd_pgp_
    if (present(iobs)) then
!       The parameter to be saved is global but an observation index was given
        ASSERT(parind(ip) .gt. 0)
        call codent(iobs, 'G', k_iobs)
        savejv_(9:15) = '.'//k_iobs(1:6)
    end if
    savejv_(16:24) = '.'//param_
!
!   ====================================================================
!   = 2 = Saving data
!   ====================================================================
!
!   --- Vectors
    if (abs(parind(ip)) .eq. 2) then
!
!       The parameter to be saved is a vector but no vector input was found
        ASSERT(UN_PARMI4(kvect, ivect, rvect, cvect))
        ASSERT(lonvec .ge. 1)
!
        call wkvect(savejv_, 'V V '//partyp(ip), lonvec, jvect)
        if (partyp(ip) .eq. 'K24') then
            do i = 1, lonvec
                zk24(jvect+i-1) = kvect_(i)
            end do
        else if (partyp(ip) .eq. 'R8') then
            b_n = to_blas_int(lonvec)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rvect, b_incx, zr(jvect), b_incy)
        else if (partyp(ip) .eq. 'C8') then
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
        ASSERT(UN_PARMI4(kscal, iscal, rscal, cscal))
!
        call wkvect(savejv_, 'V V '//partyp(ip), 1, jscal)
        if (partyp(ip) .eq. 'K24') then
            zk24(jscal) = kscal_
        else if (partyp(ip) .eq. 'R8') then
            zr(jscal) = rscal
        else if (partyp(ip) .eq. 'C8') then
            zc(jscal) = cscal
        else if (partyp(ip) .eq. 'I') then
            zi(jscal) = iscal
        end if
!
    end if
!
    if (present(savejv)) savejv = savejv_
    if (present(kvect)) AS_DEALLOCATE(vk24=kvect_)
!
    call jedema()
!
end subroutine
