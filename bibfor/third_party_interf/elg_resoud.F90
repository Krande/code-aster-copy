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

subroutine elg_resoud(matas1, matpre, nsecm, chsecm, chsolu, &
                      base, rsolu, csolu, criter, prepos, &
                      istop, iret)
    use elg_data_module
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elg_calc_rhs_red.h"
#include "asterfort/elg_calc_solu.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/resou1.h"
#include "asterfort/vtcrem.h"
#include "asterfort/detrsd.h"
!-----------------------------------------------------------------------
! But : faire "resoud" si ELIM_LAGR='OUI'
!-----------------------------------------------------------------------
!
    character(len=19), intent(in) :: matas1
    character(len=*), intent(in) :: matpre
    integer(kind=8), intent(in) :: nsecm
    character(len=*), intent(in) :: chsecm
    character(len=*), intent(in) :: chsolu
    character(len=*), intent(in) :: base
    real(kind=8), intent(inout) :: rsolu(*)
    complex(kind=8), intent(inout) :: csolu(*)
    character(len=*), intent(in) :: criter
    aster_logical, intent(in) :: prepos
    integer(kind=8), intent(in) :: istop
    integer(kind=8), intent(out) :: iret
    integer(kind=8) :: jsolu1, jsolu2, nsecmb
    character(len=19) :: matas2, secm19, solve2, solu19, chcin2
    character(len=24) :: solu2
    real(kind=8), pointer :: secm(:) => null()
    character(len=24), pointer :: refa(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
!
    call jeveuo(matas1//'.REFA', 'L', vk24=refa)
    matas2 = refa(19) (1:19)
    ASSERT(matas2 .ne. ' ')
    call dismoi('SOLVEUR', matas2, 'MATR_ASSE', repk=solve2)
!
!   -- mise a  jour du COMMON ELIMLG :
!   ---------------------------------------------
    call elg_gest_data('CHERCHE', matas1, matas2, ' ')
!
!
!   -- Second-membre(s) : passage complet -> reduit :
!   --------------------------------------------------
    solu2 = '&&ELG_RESOUD.SOLU2'
    if (nsecm .eq. 0) then
        secm19 = chsecm
        call jeveuo(secm19//'.VALE', 'L', vr=secm)
        call elg_calc_rhs_red(matas1, 1_8, secm, solu2)
    else
        ASSERT(nsecm .gt. 0)
        call elg_calc_rhs_red(matas1, nsecm, rsolu, solu2)
    end if
!
!
!   -- On cree chcin2 (qui est nul) car ELIM_LAGR n'admet pas AFFE_CHAR_CINE :
!   ----------------------------------------------------------------------------
    chcin2 = '&&elg_resoud.cine'
    call vtcrem(chcin2, matas2, 'V', 'R')

!
!   -- on appelle resou1 avec le(s) second-membre(s) reduit(s) :
!   ------------------------------------------------------------
    call jeveuo(solu2, 'E', jsolu2)
    nsecmb = max(nsecm, 1)
    call resou1(matas2, matpre, solve2, chcin2, nsecmb, &
                ' ', ' ', 'V', zr(jsolu2), csolu, &
                criter, prepos, istop, iret)
    call detrsd('CHAMP', chcin2)
!
!
!   -- Solution(s) : passage reduit -> complet :
!   ---------------------------------------------
    if (nsecm .eq. 0) then
        solu19 = chsolu
        call copisd('CHAMP', base, chsecm, solu19)
        call jeveuo(solu19//'.VALE', 'E', jsolu1)
        call elg_calc_solu(matas1, 1_8, zr(jsolu2), zr(jsolu1))
    else
        call elg_calc_solu(matas1, nsecm, zr(jsolu2), rsolu)
    end if
    call jedetr(solu2)
!
!
!
    call jedema()
end subroutine
