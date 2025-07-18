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

subroutine ccchcr(crit, name_gd, nb_val_in, val_in, cmp_in, &
                  nb_cmp_out, val_out, ichk)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/lciv2e.h"
#include "asterfort/lciv2s.h"
#include "asterfort/utmess.h"
!
!
    character(len=16), intent(in) :: crit
    character(len=8), intent(in) :: name_gd
    integer(kind=8), intent(in) :: nb_val_in
    real(kind=8), intent(in) :: val_in(nb_val_in)
    character(len=8), intent(in) :: cmp_in(nb_val_in)
    integer(kind=8), intent(in) :: nb_cmp_out
    real(kind=8), intent(out) :: val_out(nb_cmp_out)
    integer(kind=8), intent(out) :: ichk
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP - CHAM_UTIL
!
! Compute CRITERE
!
! --------------------------------------------------------------------------------------------------
!
! In  crit       : type of criterion
! In  name_gd    : name of <GRANDEUR> of input field
! In  nb_val_in  : number of input values
! In  val_in     : input values
! In  cmp_in     : name of input components
! In  nb_cmp_out : number of output values
! Out val_out    : output values
! Out ichk       : 0 if OK
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ncsig
    parameter(ncsig=6)
    character(len=4) :: cmpsig(ncsig), cmpeps(ncsig), nomcmp(ncsig)
!
    integer(kind=8) :: i
    real(kind=8) :: rac2, vale(6)
    character(len=16) :: valk(2)
!   COMMON LCIV2S/LCIV2E
    integer(kind=8) :: nt, nd
    common/tdim/nt, nd
!
    data cmpsig/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
    data cmpeps/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
!
! --------------------------------------------------------------------------------------------------
!
    ichk = 1
    rac2 = sqrt(2.d0)
! - COMMON LCIV2S/LCIV2E
    nt = nb_val_in
    nd = 3
!
! - Checking
!
    if (crit .eq. 'VMIS' .or. crit .eq. 'INVA_2' .or. crit .eq. 'TRACE') then
        ASSERT(nb_cmp_out .eq. 1)
        if (name_gd .eq. 'SIEF_R') then
            nomcmp = cmpsig
        else if (name_gd .eq. 'EPSI_R') then
            nomcmp = cmpeps
        else
            valk(1) = crit
            valk(2) = name_gd
            call utmess('F', 'CHAMPS_7', nk=2, valk=valk)
        end if
!       CELCES ASSURE QUE LES COMPOSANTES SONT DANS L'ORDRE DU CATALOGUE
!       LE MEME QUE CELUI DE CMPSIG/CMPEPS
!       IL SUFFIT DONC DE VERIFIER QUE :
!       - SI NBCMP=6, CE SONT CELLES DE CMPSIG/CMPEPS
!       - SI NBCMP=4, CE SONT LES 4 PREMIERES
        if (nb_val_in .ne. 4 .and. nb_val_in .ne. 6) then
            goto 999
        end if
        do i = 1, nb_val_in
            if (nomcmp(i) .ne. cmp_in(i)) then
                goto 999
            end if
        end do
    end if
!
! - VMIS
    if (crit .eq. 'VMIS') then
        if (name_gd .ne. 'SIEF_R') then
            valk(1) = crit
            valk(2) = name_gd
            call utmess('F', 'CHAMPS_7', nk=2, valk=valk)
        end if
        do i = 1, nd
            vale(i) = val_in(i)
        end do
        do i = nd+1, nt
            vale(i) = rac2*val_in(i)
        end do
        val_out(1) = lciv2s(vale)
        ichk = 0
!
! - INVA_2
    else if (crit .eq. 'INVA_2') then
        if (name_gd .ne. 'EPSI_R') then
            valk(1) = crit
            valk(2) = name_gd
            call utmess('F', 'CHAMPS_7', nk=2, valk=valk)
        end if
        do i = 1, nd
            vale(i) = val_in(i)
        end do
        do i = nd+1, nt
            vale(i) = rac2*val_in(i)
        end do
        val_out(1) = lciv2e(vale)
        ichk = 0
!
!-- TRACE
    else if (crit .eq. 'TRACE') then
        val_out(1) = val_in(1)+val_in(2)+val_in(3)
        ichk = 0
!
    else
        ASSERT(.false.)
    end if
!
999 continue
!
end subroutine
