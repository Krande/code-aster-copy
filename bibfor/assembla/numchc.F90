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

subroutine numchc(nu, ccid, nbchc, lchci, base)
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: lchci(*), nu, ccid
    character(len=1) :: base
    integer(kind=8) :: nbchc
!-----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
!  BUT : ON NOTE LES DDLS ELIMINES PAR LES CHARGES CINEMATIQUES
!
!  REMARQUE : LE RESTE DU TRAITEMENT DES CHARGES CINEMATIQUES EST FAIT
!             AU DERNIER MOMENT (ASMCHC+CSMBGG)
!
!-----------------------------------------------------------------------
! IN   NOMNU   K*14    : NOM DE LA NUMEROTATION
! IN   NBCHC   I       : NOMBRE DE CHARGE CINEMATIQUES
! IN   LCHCI   K*      : LISTE DES NOMS DES CHARGES CINEMATIQUES
!                        L'EFFET DE CES CHARGES EST CUMULE DANS MATAS
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!     VARIABLES LOCALES
!----------------------------------------------------------------------
    character(len=8) :: gd
    character(len=14) :: num_ddl
    character(len=19) :: nomch
    integer(kind=8) :: neq, numgd, iddes, nec, jccid, idprno
    integer(kind=8) :: nelim, jafci, nimp, imp, ino, iddl, ieq, ich
    integer(kind=8), pointer :: nequ(:) => null()
!----------------------------------------------------------------------
!
    call jemarq()
    num_ddl = nu
!
    call jeveuo(num_ddl//'.NUME.NEQU', 'L', vi=nequ)
    neq = nequ(1)
    call dismoi('NOM_GD', num_ddl, 'NUME_DDL', repk=gd)
    call jenonu(jexnom('&CATA.GD.NOMGD', gd), numgd)
    call jeveuo(jexnum('&CATA.GD.DESCRIGD', numgd), 'L', iddes)
    nec = zi(iddes+2)
!
    call wkvect(ccid, base//' V I ', neq, jccid)
!
    if (nbchc .eq. 0) goto 40

!
!     -- IL N'Y A PEUT-ETRE AUCUN DDL A ELIMINER (CHAR_CINE VIDES) :
    nimp = 0
    do ich = 1, nbchc
        nomch = lchci(ich)
        call jeveuo(nomch//'.AFCI', 'L', jafci)
        nimp = nimp+zi(jafci)
    end do
!
    if (nimp .eq. 0) goto 40
!
!
    call jeveuo(jexnum(num_ddl//'.NUME.PRNO', 1), 'L', idprno)
    nelim = 0
    do ich = 1, nbchc
        nomch = lchci(ich)
        call jeveuo(nomch//'.AFCI', 'L', jafci)
        nimp = zi(jafci)
        do imp = 1, nimp
            ino = zi(jafci+3*(imp-1)+1)
            iddl = zi(jafci+3*(imp-1)+2)
            ieq = zi(idprno-1+(nec+2)*(ino-1)+1)+iddl-1
            zi(jccid-1+ieq) = 1
        end do
    end do
!
!
40  continue
    call jedema()
end subroutine
