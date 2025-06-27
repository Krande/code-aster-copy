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

subroutine op0001()
!
!-----------------------------------------------------------------------
!   LIRE_MAILLAGE with FORMAT='ASTER'
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/chckma.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/infoma.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lrmast.h"
#include "asterfort/mavegr.h"
#include "asterfort/wkvect.h"

!
! ----- DECLARATIONS
!
    integer(kind=8) :: iaux, niv, ifl, ifm
    integer(kind=8) :: nbnoeu, nbmail, nbcoor
    integer(kind=8) :: iret, infmed
    character(len=8) :: nomu, fmt, veri
    character(len=16) :: concep, cmd
    character(len=64) :: nomamd
    real(kind=8) :: dtol
    integer(kind=8), pointer :: dime(:) => null()

!---------------------------------------------------------------------------------------
    call jemarq()
!
! --- RECUPERATION DES ARGUMENTS  DE LA COMMANDE
!
    ifl = 0
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomu, concep, cmd)
!
    call getvis(' ', 'UNITE', scal=ifl, nbret=iaux)
!
    call getvtx(' ', 'FORMAT', scal=fmt, nbret=iaux)
    ASSERT(fmt .eq. 'ASTER')
!
! --- LECTURE DU MAILLAGE AU FORMAT ASTER :
!     -----------------------------------
    call lrmast(nomu, ifm, ifl, nbnoeu, nbmail, &
                nbcoor)

! --- SUPPRESSION DES GROUPES DE NOEUDS OU MAILLES DE NOM ' ' :
!     -------------------------------------------------------
    call mavegr(nomu)

! --- CREATION DE L'OBJET .DIME :
!     -------------------------
    call wkvect(nomu//'.DIME', 'G V I', 6, vi=dime)
    dime(1) = nbnoeu
    dime(3) = nbmail
    dime(6) = nbcoor

! --- CARACTERISTIQUES GEOMETRIQUES :
!     -----------------------------
    call cargeo(nomu)

! --- PHASE DE VERIFICATION DU MAILLAGE :
!     ---------------------------------
    call getvtx('VERI_MAIL', 'VERIF', iocc=1, scal=veri, nbret=iret)
    if (veri .eq. 'OUI') then
        call getvr8('VERI_MAIL', 'APLAT', iocc=1, scal=dtol, nbret=iret)
        call chckma(nomu, dtol)
    end if

!     IMPRESSIONS DU MOT CLE INFO :
!     ---------------------------
    call infoma(nomu)
!
    call jedema()
end subroutine
