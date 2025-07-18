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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmcrdn(sd_suiv, keyw_fact, nb_dof_monitor, nb_keyw_fact)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/wkvect.h"
!
    character(len=24), intent(in) :: sd_suiv
    character(len=16), intent(in) :: keyw_fact
    integer(kind=8), intent(in) :: nb_dof_monitor
    integer(kind=8), intent(in) :: nb_keyw_fact
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear operators - DOF monitor
!
! Create datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  sd_suiv          : datastructure for dof monitor parameters
! In  keyw_fact        : factor keyword to read extraction parameters
! In  nb_keyw_fact     : number of factor keyword to read extraction parameters
! In  nb_dof_monitor   : number of factor dofs monitored
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_keyw_fact, i_dof_monitor, nb_line_title
    character(len=16) :: title(3)
    character(len=1) :: chaine
    character(len=24) :: dofm_titl
    character(len=16), pointer :: v_dofm_titl(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

!
! - Create title vector
!
    dofm_titl = sd_suiv(1:14)//'     .TITR'
    call wkvect(dofm_titl, 'V V K16', 3*nb_dof_monitor, vk16=v_dofm_titl)
!
! - Title from user
!
    do i_keyw_fact = 1, nb_keyw_fact
        write (chaine, '(I1)') i_keyw_fact
        title(1) = '    SUIVI '
        title(2) = '     DDL  '
        title(3) = '     '//chaine
        call getvtx(keyw_fact, 'TITRE', iocc=i_keyw_fact, nbval=0, nbret=nb_line_title)
        nb_line_title = -nb_line_title
        ASSERT(nb_line_title .le. 3)
        if (nb_line_title .ne. 0) then
            call getvtx(keyw_fact, 'TITRE', iocc=i_keyw_fact, nbval=nb_line_title, vect=title)
        end if
        v_dofm_titl(3*(i_keyw_fact-1)+1) = title(1)
        v_dofm_titl(3*(i_keyw_fact-1)+2) = title(2)
        v_dofm_titl(3*(i_keyw_fact-1)+3) = title(3)
    end do
!
! - Automatic
!
    if (nb_dof_monitor .gt. nb_keyw_fact) then
        do i_dof_monitor = nb_keyw_fact+1, nb_dof_monitor
            write (chaine, '(I1)') i_dof_monitor
            title(1) = '    SUIVI '
            title(2) = '     DDL  '
            title(3) = '     '//chaine
            v_dofm_titl(3*(i_dof_monitor-1)+1) = title(1)
            v_dofm_titl(3*(i_dof_monitor-1)+2) = title(2)
            v_dofm_titl(3*(i_dof_monitor-1)+3) = title(3)
        end do
    end if
!
end subroutine
