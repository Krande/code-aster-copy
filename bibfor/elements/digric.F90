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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine digric(DD, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    DD      : voir l'appel
! OUT   iret    : code retour
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
    type(te0047_dscr), intent(in) :: DD
    integer(kind=8), intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iadzi, iazk24, imat, ivarim, ifono, icontp, ivarip
    real(kind=8) :: klv(78)
    character(len=24) :: messak(5)
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!
    if (DD%nomte .ne. 'MECA_DIS_TR_L') then
        messak(1) = DD%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = DD%type_comp
        messak(4) = DD%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_11', nk=5, valk=messak)
    end if
    ! Paramètres en entrée
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
    !
    ifono = 1
    icontp = 1
    ivarip = 1
    if (DD%lVect) then
        call jevech('PVECTUR', 'E', ifono)
    end if
    if (DD%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (DD%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    !
    call dicrgr(DD, zi(imat), zr(ivarim), klv, zr(ivarip), zr(ifono), zr(icontp))
    !
    if (DD%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        call utpslg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
    end if
end subroutine
