! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine digric(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
use te0047_type
implicit none
!
#include "jeveux.h"
#include "asterfort/dicrgr.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
!
type(te0047_dscr), intent(in) :: for_discret
integer, intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iadzi, iazk24, imat, ivarim, jtp, ifono, icontp, ivarip, neq, icontm
    real(kind=8) :: klv(78)
    character(len=24) :: messak(5)
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    call jevech('PCONTMR', 'L', icontm)
!
    if (for_discret%nomte .ne. 'MECA_DIS_TR_L') then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_11', nk=5, valk=messak)
    endif
!   paramètres en entrée
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTPR', 'L', jtp)
!
    ifono = 1
    icontp = 1
    ivarip = 1
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
    endif
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    endif
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    endif
    neq = for_discret%nno*for_discret%nc
!
    call dicrgr('RIGI', for_discret%option, neq, for_discret%nc, zi(imat), &
                for_discret%ulm, for_discret%dul, zr(icontm), zr(ivarim), for_discret%pgl, &
                klv, zr(ivarip), zr(ifono), zr(icontp))
!
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
    endif
end subroutine
