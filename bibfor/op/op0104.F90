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

subroutine op0104()
    implicit none
!     OPERATEUR: DEFI_GROUP
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/addGroupElem.h"
#include "asterfort/addGroupNode.h"
#include "asterfort/cgrcbp.h"
#include "asterfort/detgnm.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/sscgma.h"
#include "asterfort/sscgno.h"
#include "asterfort/utmess.h"
!
!
    integer(kind=8) :: n1, n2, nbgrma, nbgmin, iret, nbgma
    integer(kind=8) :: nbocc, nbgrno, iocc, nbgnin, n3
    character(len=8) :: k8b, ma, ma2
    character(len=16) :: nomcmd, typcon, option
    character(len=24) :: grpmai, grpnoe, grpmav, grpnov, gpptnm, gpptnn
    aster_logical :: l_write
!     ------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(ma2, typcon, nomcmd)
    call getvid(' ', 'MAILLAGE', scal=ma, nbret=n1)
    if (n1 .eq. 0) then
        call getvid(' ', 'GRILLE', scal=ma, nbret=n1)
    end if
    if (ma .ne. ma2) then
        call utmess('F', 'SOUSTRUC_15')
    end if
    grpmai = ma//'.GROUPEMA'
    grpnoe = ma//'.GROUPENO'
    gpptnm = ma//'.PTRNOMMAI'
    gpptnn = ma//'.PTRNOMNOE'
    grpmav = '&&OP0104'//'.GROUPEMA'
    grpnov = '&&OP0104'//'.GROUPENO'
!
!
!     1. ON DETRUIT SI DEMANDE (DETR_GROUP_MA ET DETR_GROUP_NO):
!     ------------------------------------------------------------
    call getfac('DETR_GROUP_MA', n1)
    call getfac('DETR_GROUP_NO', n2)
    if (n1 .ne. 0 .or. n2 .ne. 0) then
        call detgnm(ma)
    end if
!
!
!
!     2. CREA_GROUP_MA :
!     ------------------------------------------------------------
!     --- ON COMPTE LE NOMBRE DE NOUVEAUX GROUP_MA :
    call getfac('CREA_GROUP_MA', nbgrma)
    if (nbgrma .eq. 0) goto 107
!
!     --- ON AGRANDIT LA COLLECTION SI NECESSAIRE :
    call addGroupElem(ma, nbgrma)
    nbgmin = 0
    call jeexin(grpmai, iret)
    if (iret > 0) then
        call jelira(grpmai, 'NOMUTI', nbgmin)
    end if
107 continue
!
!
!     3. CREA_GROUP_NO :
!     ------------------------------------------------------------
!     --- ON COMPTE LE NOMBRE DE NOUVEAUX GROUP_NO :
    call getfac('CREA_GROUP_NO', nbocc)
    nbgrno = 0
    do iocc = 1, nbocc
        call getvtx('CREA_GROUP_NO', 'TOUT_GROUP_MA', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call jelira(grpmai, 'NMAXOC', nbgma)
            nbgrno = nbgrno+nbgma
            goto 10
        end if
        call getvem(ma, 'GROUP_MA', 'CREA_GROUP_NO', 'GROUP_MA', iocc, &
                    0, k8b, n2)
        if (n2 .ne. 0) then
            nbgrno = nbgrno-n2
            goto 10
        end if
        call getvtx('CREA_GROUP_NO', 'OPTION', iocc=iocc, nbval=0, nbret=n3)
        if (n3 .ne. 0) then
            call getvtx('CREA_GROUP_NO', 'OPTION', iocc=iocc, scal=option, nbret=n3)
            if (option .eq. 'RELA_CINE_BP') then
                l_write = .false.
                call cgrcbp('CREA_GROUP_NO', iocc, ma, l_write, nbgma)
                nbgrno = nbgrno+nbgma
                goto 10
            end if
        end if
!        -- ON CREE UN GROUP_NO PAR MOT CLE FACTEUR --
        nbgrno = nbgrno+1
10      continue
    end do
    if (nbgrno .eq. 0) goto 207
!
!
!     --- ON AGRANDIT LA COLLECTION SI NECESSAIRE :
    call addGroupNode(ma, nbgrno)
    nbgnin = 0
    call jeexin(grpnoe, iret)
    if (iret > 0) then
        call jelira(grpnoe, 'NOMUTI', nbgnin)
    end if
207 continue
!
!     --- TRAITEMENT DU MOT CLEF CREA_GROUP_MA :
    if (nbgrma .gt. 0) call sscgma(ma, nbgrma, nbgmin)
!
!     --- TRAITEMENT DU MOT CLEF CREA_GROUP_NO :
    if (nbgrno .gt. 0) call sscgno(ma, nbgnin)
!
! --- Pour un ParallelMesh, il faut partager les groupes entre les domaines.
!     C'est fait dans le post_exec
!
!
    call jedema()
end subroutine
