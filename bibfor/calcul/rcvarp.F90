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

subroutine rcvarp(arret, varc_name_, poum, varc_vale, iret)
!
    use calcul_module, only: ca_iredec_, ca_jvcnom_, ca_jvcval_, ca_nbcvrc_, &
                             ca_td1_, ca_tf1_, ca_timed1_, ca_timef1_
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=1), intent(in) :: arret
    character(len=*), intent(in) :: varc_name_
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(out) :: iret
    real(kind=8), intent(out) :: varc_vale
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Get value of external state variable when used in SIMU_POINT_MAT
!
! --------------------------------------------------------------------------------------------------
!
! In  arret            : in case of problem
!                     ' '   no message and output return code (iret)
!                     'F'   fatal error
! In  varc_name        : name of external state variable
! In  poum             : when get external state variables
!                     '-'   at beginning of step time
!                     '+'   at end of step time
!                     'REF' for reference value
! Out varc_vale        : value of external state variable
! Out iret             : code if error
!                      0    if OK
!                      1    if not found
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: varc_name
    integer(kind=8) :: varc_indx
    integer(kind=8), save :: iprem = 0
    real(kind=8) :: valvrm, valvrp, tdef
    real(kind=8), save :: rundf
!
! --------------------------------------------------------------------------------------------------
!
    if (iprem .eq. 0) then
        rundf = r8nnem()
        iprem = 1
    end if
    tdef = rundf
    iret = 0
!
! - Get index of external state variable
!
    varc_name = varc_name_
    varc_indx = indik8(zk8(ca_jvcnom_), varc_name, 1, ca_nbcvrc_)
!
! - Not found: NaN
!
    if (varc_indx .eq. 0) then
        iret = 1
        if (arret .eq. ' ') then
            varc_vale = rundf
            goto 999
        else
            call utmess('F', 'CALCUL_50', sk=varc_name)
        end if
    end if
!
! - Get value
!
    if (poum .eq. 'REF') then
        varc_vale = zr(ca_jvcval_-1+3*(varc_indx-1)+3)

    else if (poum .eq. '+' .and. ca_iredec_ .eq. 0) then
        varc_vale = zr(ca_jvcval_-1+3*(varc_indx-1)+2)

    else if (poum .eq. '-' .and. ca_iredec_ .eq. 0) then
        varc_vale = zr(ca_jvcval_-1+3*(varc_indx-1)+1)

    else if (ca_iredec_ .eq. 1) then
        valvrm = zr(ca_jvcval_-1+3*(varc_indx-1)+1)
        valvrp = zr(ca_jvcval_-1+3*(varc_indx-1)+2)
        if ((.not. isnan(valvrm)) .and. (.not. isnan(valvrp))) then
            if (poum .eq. '-') then
                varc_vale = valvrm+(ca_td1_-ca_timed1_)*(valvrp-valvrm)/(ca_timef1_-ca_timed1_)
            else if (poum .eq. '+') then
                varc_vale = valvrm+(ca_tf1_-ca_timed1_)*(valvrp-valvrm)/(ca_timef1_-ca_timed1_)
            else
                ASSERT(.false.)
            end if
        else
            varc_vale = rundf
        end if

    else
        ASSERT(.false.)
    end if
!
    iret = 0
    if (isnan(varc_vale)) then
        iret = 1
    end if
!
! - Manage error
!
    if (iret .eq. 1) then
        if (varc_name .eq. 'TEMP') then
            varc_vale = tdef
            iret = 1
            goto 999
        end if
        if (arret .eq. ' ') then
            varc_vale = rundf
        else
            call utmess('F', 'CALCUL_50', sk=varc_name)
        end if
    end if
!
999 continue
!
end subroutine
