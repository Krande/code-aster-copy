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

subroutine op0090()
    implicit none
! person_in_charge: mathieu.courtois at edf.fr
!     OPERATEUR "RECU_FONCTION"
!     ------------------------------------------------------------------
#include "asterfort/chpve2.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/rfbefl.h"
#include "asterfort/rfinte.h"
#include "asterfort/rfnapp.h"
#include "asterfort/rfnoch.h"
#include "asterfort/rfrcha.h"
#include "asterfort/rfresu.h"
#include "asterfort/rfrgen.h"
#include "asterfort/rftabl.h"
    integer(kind=8) :: nreg, nrb, nch, ng, ier
    integer(kind=8) :: nta, nres, nc, nna
    character(len=8) :: k8b
    character(len=19) :: cham19, resu, tabres, tabtyp(8), nappe
    data tabtyp/'NOEU#DEPL_R', 'NOEU#TEMP_R', 'NOEU#PRES_R',&
     &            'ELXX#SIEF_R', 'ELXX#VARI_R', 'ELXX#EPSI_R',&
     &            'ELXX#FLUX_R', 'ELXX#PRES_R'/
!     ------------------------------------------------------------------
!
!     -----------------------------------------------------------------
!                      --- CAS D'UN CHAM_GD ---
!     -----------------------------------------------------------------
    call getvid(' ', 'CHAM_GD', scal=cham19, nbret=nch)
    if (nch .ne. 0) then
        call chpve2(cham19, 8, tabtyp, ier)
        call rfrcha()
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                       --- CAS D'UN RESULTAT ---
!     -----------------------------------------------------------------
    call getvid(' ', 'RESULTAT ', scal=resu, nbret=nres)
    if (nres .ne. 0) then
        call rfresu()
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                   --- CAS D'UN NOEUD DE CHOC ---
!     -----------------------------------------------------------------
    call getvtx(' ', 'NOEUD_CHOC', scal=k8b, nbret=nc)
    call getvtx(' ', 'GROUP_NO_CHOC', scal=k8b, nbret=ng)
    if (nc+ng .ne. 0) then
        call rfnoch()
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                    --- CAS D'UN RESU_GENE ---
!     -----------------------------------------------------------------
    call getvid(' ', 'RESU_GENE', scal=resu, nbret=nreg)
    if (nreg .ne. 0) then
        call rfrgen(resu)
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                       --- CAS D'UNE TABLE ---
!     -----------------------------------------------------------------
    call getvid(' ', 'TABLE', scal=tabres, nbret=nta)
    if (nta .ne. 0) then
        call rftabl(tabres)
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                 --- CAS D'UNE BASE_ELAS_FLUI ---
!     -----------------------------------------------------------------
    call getvid(' ', 'BASE_ELAS_FLUI', scal=resu, nbret=nrb)
    if (nrb .ne. 0) then
        call rfbefl(resu)
        goto 10
    end if
!
!     -----------------------------------------------------------------
!                 --- CAS D'UNE SD_INTERSPECTRE ---
!     -----------------------------------------------------------------
    call getvid(' ', 'INTE_SPEC', scal=resu, nbret=nrb)
    if (nrb .ne. 0) then
        call rfinte(resu)
        goto 10
    end if
!
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!                     --- CAS D'UNE NAPPE ---
!     -----------------------------------------------------------------
    call getvid(' ', 'NAPPE', scal=nappe, nbret=nna)
    if (nna .ne. 0) then
        call rfnapp(nappe)
        goto 10
    end if
!
!     -----------------------------------------------------------------
10  continue
!
end subroutine
