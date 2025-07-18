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
subroutine char_read_tran(keywordfact, iocc, ndim, &
                          l_tran_, tran_, &
                          l_cent_, cent_, &
                          l_angl_naut_, angl_naut_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8dgrd.h"
#include "asterfort/getvr8.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    integer(kind=8), intent(in) :: ndim
    aster_logical, optional, intent(out) :: l_tran_
    real(kind=8), optional, intent(out) :: tran_(3)
    aster_logical, optional, intent(out) :: l_cent_
    real(kind=8), optional, intent(out) :: cent_(3)
    aster_logical, optional, intent(out) :: l_angl_naut_
    real(kind=8), optional, intent(out) :: angl_naut_(3)
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_CHAR_MECA
!
! Read transformation
!
! --------------------------------------------------------------------------------------------------
!
! In  keywordfact      : factor keyword to read
! In  iocc             : factor keyword index in AFFE_CHAR_MECA
! In  ndim             : space dimension
! Out l_tran           : .true. if TRAN defined (translation)
! Out tran             : vector defining translation
! Out l_cent           : .true. if center defined (rotation)
! Out cent             : vector defining center
! Out l_angl_naut      : .true. if angl defined (rotation)
! Out angl_naut        : angle defining rotation
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nangmx, i
    integer(kind=8) :: ntran, ncent, nangl, vali(2)
    aster_logical :: l_tran, l_cent, l_angl_naut
    real(kind=8) :: tran(3), cent(3), angl_naut(3)
!
! --------------------------------------------------------------------------------------------------
!
    l_tran = ASTER_FALSE
    l_angl_naut = ASTER_FALSE
    l_cent = ASTER_FALSE
    tran(1:3) = 0.d0
    cent(1:3) = 0.d0
    angl_naut(1:3) = 0.d0
    if (ndim .eq. 3) then
        nangmx = 3
    else
        nangmx = 1
    end if

! - Translation
    call getvr8(keywordfact, 'TRAN', iocc=iocc, nbval=0, nbret=ntran)
    ntran = -ntran
    if (ntran .ne. 0) then
        l_tran = ASTER_TRUE
        if (ntran .ne. ndim) then
            vali(1) = ndim
            vali(2) = ndim
            call utmess('F', 'CHARGES2_42', sk='TRAN', ni=2, vali=vali)
        end if
        call getvr8(keywordfact, 'TRAN', iocc=iocc, nbval=ntran, vect=tran)
    end if

! - Rotation
    call getvr8(keywordfact, 'CENTRE', iocc=iocc, nbval=0, nbret=ncent)
    ncent = -ncent
    if (ncent .ne. 0) then
        l_cent = ASTER_TRUE
        if (ncent .ne. ndim) then
            vali(1) = ndim
            vali(2) = ndim
            call utmess('F', 'CHARGES2_42', sk='CENTRE', ni=2, vali=vali)
        end if
        call getvr8(keywordfact, 'CENTRE', iocc=iocc, nbval=ncent, vect=cent)
    end if
!
    call getvr8(keywordfact, 'ANGL_NAUT', iocc=iocc, nbval=0, nbret=nangl)
    nangl = -nangl
    if (nangl .ne. 0) then
        l_angl_naut = .true.
        if (nangl .ne. nangmx) then
            vali(1) = nangmx
            vali(2) = ndim
            call utmess('F', 'CHARGES2_42', sk='ANGL_NAUT', ni=2, vali=vali)
        end if
        call getvr8(keywordfact, 'ANGL_NAUT', iocc=iocc, nbval=nangmx, vect=angl_naut)
        do i = 1, 3
            angl_naut(i) = angl_naut(i)*r8dgrd()
        end do
    end if

! - Copy output
    if (present(l_tran_)) then
        l_tran_ = l_tran
        tran_(1:3) = tran(1:3)
    end if
    if (present(l_cent_)) then
        l_cent_ = l_cent
        cent_(1:3) = cent(1:3)
    end if
    if (present(l_angl_naut_)) then
        l_angl_naut_ = l_angl_naut
        angl_naut_(1:3) = angl_naut(1:3)
    end if
!
end subroutine
