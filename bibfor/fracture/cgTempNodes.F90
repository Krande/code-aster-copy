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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgTempNodes(cgStudy, cgTable)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/vrcins.h"
#include "asterfort/exisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/chpchd.h"
#include "asterfort/detrsd.h"
#include "asterfort/imprsd.h"
#include "asterfort/jemarq.h"
#include "asterfort/celces.h"
#include "asterfort/cescns.h"
#include "asterfort/celfpg.h"
#include "asterfort/cesred.h"
!
    type(CalcG_study), intent(in) :: cgStudy
    type(CalcG_table), intent(inout) :: cgTable
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Get tempature on the crack
!
! IO compor    : name of <CARTE> COMPORTEMENT
! --------------------------------------------------------------------------------------------------
!
    !   character(len=19), parameter :: celvrc = '&&TABLEG.CELVRC'
    ! character(len=19), parameter :: cesvrc = '&&TABLEG.CESVRC'
    ! character(len=19), parameter :: redvrc = '&&TABLEG.REDVRC'
!    character(len=19), parameter :: cnovrc = '&&TABLEG.CNOVRC'
    ! character(len=19), parameter :: cnsvrc = '&&TABLEG.CNSVRC'
    ! character(len=19), parameter :: fpgvrc = '&&TABLEG.CELPFG'
!    character(len=8) :: tych
!    character(len=2) :: codret
    character(len=3) :: repk
    ! integer :: iret, nbma, i
    ! integer, pointer :: listma(:) => null()
!
    call jemarq()
!
    call dismoi('EXI_VARC', cgStudy%material, 'CHAM_MATER', repk=repk)

!    print*, cgStudy%model, cgStudy%material, cgStudy%carael, cgStudy%time
    ! codret = 'XX'
    ! call vrcins(cgStudy%model, cgStudy%material, "       ", cgStudy%time, celvrc, codret)
    ! !, nompaz='PVARCNO')
    ! if( codret .ne. 'OK' ) then
    !     go to 999
    ! end if

! pour le moment on oublie la suite
    print *, cgTable%table_g
    go to 999

! On n'a un problème avec les éléments tardifs - passez par un TE ?
!     call celces(celvrc, 'V', cesvrc)
!     call celfpg(celvrc, fpgvrc, iret)
!     ASSERT( iret == 0 )

!     call dismoi('NB_MA_MAILLA', cgStudy%mesh, 'MAILLAGE', repi=nbma)
!     allocate(listma(nbma))
!     do i = 1, 10
!         listma(i) = i
!     end do

!     call cesred(cesvrc, 10, listma, 0, ["XXX"],'V', redvrc)
!     call imprsd("CHAMP", redvrc, 6, "VARINS")

!     call cescns(redvrc, fpgvrc, 'V', cnsvrc, 'F', iret)
! !    call imprsd("CHAMP", celvrc, 6, "VARINS")
!     deallocate(listma)
!
!     call exisd('CHAM_ELEM', celvrc, iret)
!     ASSERT( iret .eq. 1 )
!     call dismoi('TYPE_CHAMP', celvrc, 'CHAMP', repk=tych)
!     ASSERT( tych .eq. 'ELGA' )
! !   transformation cham_elem (ELGA) -> cham_no : celvrc -> cnovrc
!     call chpchd(celvrc, 'NOEU', ' ', 'NON', 'V', cnovrc)
!     call detrsd('CHAM_ELEM', celvrc)
!     call imprsd("CHAM_NO", cnovrc, 6, "VARINS")
!
!    do i_node = 1, nb_point
!        node_nume = fondNoeudNume(i_node)
!    end do
!
999 continue
!
    call jedema()
!
end subroutine
