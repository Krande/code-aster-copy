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
subroutine nmiret(codret, tabret)
    !
    implicit none
    !
#include "asterf_types.h"
#include "asterfort/asmpi_comm_logical.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdmpic.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
    !
    aster_logical :: tabret(0:10)
    character(len=19) :: codret
    !
    ! ----------------------------------------------------------------------
    !
    ! ROUTINE MECA_NON_LINE (CALCUL)
    !
    ! RESUME LES CODES RETOURS DES TE
    !
    ! ----------------------------------------------------------------------
    !
    !
    ! IN  CODRET  : CHAM_ELEM ISSU DES TE
    ! OUT TABRET  : TABRET(0) = .TRUE. UN CODE RETOUR NON NUL EXISTE
    !                TABRET(I) = .TRUE. CODE RETOUR I RENCONTRE
    !                             SINON .FALSE.
    !                I VALANT DE 1 A 10
    !
    !
    !
    !
    integer(kind=8) :: iret, jcesd, jcesl, nbmail, icmp
    integer(kind=8) :: ima, iad, vali
    character(len=8) :: nomgd, mesh
    character(len=19) :: chamns
    integer(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    !
    ! ----------------------------------------------------------------------
    !
    call jemarq()
    !
    do iret = 0, 10
        tabret(iret) = .false.
    end do
    !
    ! --- ON TRANSFORME LE "CHAM_ELEM" EN UN "CHAM_ELEM_S"
    !
    chamns = '&&NMIRET.CHAMNS'
    call exisd('CHAM_ELEM', codret, iret)
    if (iret .eq. 0) then
        goto 99
    end if

    !
    !     -- EN ATTENDANT DE FAIRE MIEUX, POUR PERMETTRE MUMPS/DISTRIBUE :
    call sdmpic('CHAM_ELEM', codret)
    !
    call celces(codret, 'V', chamns)
    !
    ! --- ACCES AU CHAM_ELEM_S
    !
    call jeveuo(chamns//'.CESK', 'L', vk8=cesk)
    call jeveuo(chamns//'.CESD', 'L', jcesd)
    call jeveuo(chamns//'.CESV', 'L', vi=cesv)
    call jeveuo(chamns//'.CESL', 'L', jcesl)
    !
    !     CHAM_ELEM/ELGA MAIS EN FAIT : 1 POINT ET 1 SOUS_POINT PAR ELEMENT
    if ((zi(jcesd-1+3) .ne. 1) .or. (zi(jcesd-1+4) .ne. 1)) then
        ASSERT(.false.)
    end if
    !
    nomgd = cesk(2)
    if (nomgd .ne. 'CODE_I') then
        ASSERT(.false.)
    end if
    !
    nbmail = zi(jcesd-1+1)
    icmp = zi(jcesd-1+2)
    if (icmp .ne. 1) then
        ASSERT(.false.)
    end if
    !
    do ima = 1, nbmail
        call cesexi('C', jcesd, jcesl, ima, 1, 1, icmp, iad)
        if (iad .le. 0) cycle
        iret = cesv(iad)
        if (iret .eq. 0) then
        else if (iret .lt. 11 .and. iret .gt. 0) then
            tabret(iret) = .true.
        else
            vali = iret
            call utmess('A', 'MECANONLINE2_67', si=vali)
        end if
    end do
    !
    ! --- Il faut faire une synthèse en HPC
    !
    call dismoi('NOM_MAILLA', codret, 'CHAM_ELEM', repk=mesh)
    if (isParallelMesh(mesh)) then
        call asmpi_comm_logical("MPI_LOR", nbval=10, vl=tabret)
    end if
    !
    do iret = 1, 10
        if (tabret(iret)) tabret(0) = .true.
    end do
    !
    call detrsd('CHAM_ELEM_S', chamns)
    !
99  continue
    !
    call jedema()
end subroutine
