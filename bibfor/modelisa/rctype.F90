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

subroutine rctype(jmat, nb_para_list, para_list_name, para_list_vale, para_vale, &
                  para_type, keyw_factz, keywz, materi)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
!
    integer(kind=8), intent(in) :: jmat
    integer(kind=8), intent(in) :: nb_para_list
    character(len=*), intent(in) :: para_list_name(*)
    real(kind=8), intent(in) :: para_list_vale(*)
    real(kind=8), intent(out) :: para_vale
    character(len=*), intent(out) :: para_type
    character(len=*), optional, intent(in) :: keyw_factz
    character(len=*), optional, intent(in) :: keywz
    character(len=*), optional, intent(in) :: materi
!
! --------------------------------------------------------------------------------------------------
!
! Comportment - Utility
!
! Get informations (value and type) about parameters for traction curve
!
! --------------------------------------------------------------------------------------------------
!
! In  ktrac          : simple keyword for traction curve (SIGM by default)
!                       1  - 'TRACTION'
!                       2  - 'META_TRACTION'
! In  jmat           : JEVEUX adress to coded material
! In  nb_para_list   : number of possible parameters (except EPSI)
! In  para_list_name : list of possible parameters name (except EPSI)
! In  para_list_vale : list of parameters value (except EPSI)
! Out para_vale      : parameter value (except EPSI)
! Out para_type      : parameter type (except EPSI)
! In  keyw_factz     : factor keyword for traction curve (TRACTION by default)
! In  keywz          : simple keyword for traction curve (SIGM by default)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lfct, lmat, lsup
    parameter(lmat=9, lfct=10, lsup=2)
!
    integer(kind=8) :: icomp, ipi, idf, nbf, ivalk, ik, ipif, jpro, kmat, inom
    integer(kind=8) :: imate, nbmat
    character(len=24) :: para_name(2)
    character(len=16) :: keyw_fact
    character(len=8) :: keyw, nomi
    integer(kind=8) :: ipara, nb_para, i_para_list
!
! --------------------------------------------------------------------------------------------------
!
    para_type = ' '
    para_vale = 0.d0
!
    ipi = 0; ipif = 0; imate = 0
!
!   Number of materials data for current element
    nbmat = zi(jmat)
    if (nbmat .ne. 1) then
        ASSERT(present(materi))
        do kmat = 1, nbmat
            inom = zi(jmat+kmat)
            nomi = zk8(inom)
            if (nomi .eq. materi) then
!               Coded material
                imate = jmat+zi(jmat+nbmat+kmat)
                goto 5
            end if
        end do
        call utmess('F', 'CALCUL_45', sk=materi)
    else
!       Coded material
        imate = jmat+zi(jmat+nbmat+1)
    end if
5   continue
!
!   Simple keyword for traction curve
    if (present(keywz)) then
        keyw = keywz
    else
        keyw = 'SIGM'
    end if
!
!   Factor keyword for traction curve
    if (present(keyw_factz)) then
        keyw_fact = keyw_factz
    else
        keyw_fact = 'TRACTION'
    end if
!
!   Get index for factor keyword
    do icomp = 1, zi(imate+1)
        if (keyw_fact .eq. zk32(zi(imate)+icomp-1)) then
            ipi = zi(imate+2+icomp-1)
            goto 11
        end if
    end do
    call utmess('F', 'COMPOR5_1', sk=keyw_fact)
11  continue
!
! - Get index for simple keyword
!
    idf = zi(ipi)+zi(ipi+1)
    nbf = zi(ipi+2)
    ivalk = zi(ipi+3)
    do ik = 1, nbf
        if (keyw .eq. zk16(ivalk+idf+ik-1)) then
            if (keyw .eq. 'SIGM') then
                ipif = ipi+lmat-1+lfct*(ik-1)
            else if ((keyw .eq. 'SIGM_F1') .or. &
                     (keyw .eq. 'SIGM_F2') .or. &
                     (keyw .eq. 'SIGM_F3') .or. &
                     (keyw .eq. 'SIGM_F4') .or. &
                     (keyw .eq. 'SIGM_C')) then
                ipif = ipi+lmat-1+lfct*(ik-1)+lsup*(ik-1)
            else
                ASSERT(.false.)
            end if
            goto 21
        end if
    end do
    ASSERT(.false.)
21  continue
!
! - Get name of parameter(s)
!
    jpro = zi(ipif+1)
!
    if (zk24(jpro) .eq. 'NAPPE') then
        nb_para = 2
        para_name(1) = zk24(jpro+2)
        para_name(2) = zk24(jpro+5)
    else
        nb_para = 1
        para_name(1) = zk24(jpro+2)
        if (para_name(1) .eq. 'EPSI') then
            para_vale = para_list_vale(1)
            goto 999
        else
            call utmess('F', 'COMPOR5_2', sk=para_name(1))
        end if
    end if
!
! - <NAPPE> case
!
    do ipara = 1, nb_para
        if (para_name(ipara) (1:4) .ne. 'EPSI') then
            do i_para_list = 1, nb_para_list
                if (para_list_name(i_para_list) .eq. para_name(ipara)) then
                    para_vale = para_list_vale(i_para_list)
                    para_type = para_list_name(i_para_list)
                    goto 999
                end if
            end do
        end if
    end do
!
    call utmess('F', 'COMPOR5_4')
!
999 continue
!
end subroutine
