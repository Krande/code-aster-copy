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
!
interface
    subroutine maillagefibre(nogfma, ulnbnoeuds, maxmailgrp, nbgf, vcoord, nbnoeuds, &
                             vigroup, vngroup, vmailgrp, vimailles, ulnbmailles, ncarma )
        character(len=8), intent(in) :: nogfma
        integer(kind=8), intent(in) :: ulnbnoeuds
        integer(kind=8), intent(in) :: maxmailgrp
        integer(kind=8), intent(in) :: nbgf
        integer(kind=8), intent(in) :: ulnbmailles
        real(kind=8), intent(in) :: vcoord(2*ulnbnoeuds)
        integer(kind=8), intent(in) :: nbnoeuds
        integer(kind=8), intent(in) :: vmailgrp(nbgf)
        integer(kind=8), intent(in) :: vigroup(nbgf*maxmailgrp)
        character(len=24), intent(in) :: vngroup(nbgf)
        integer(kind=8), intent(in) :: ncarma
        integer(kind=8), intent(in) :: vimailles(ulnbmailles*ncarma)
    end subroutine maillagefibre
end interface
