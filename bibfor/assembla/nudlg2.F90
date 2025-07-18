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

subroutine nudlg2(nu)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: nu
! ----------------------------------------------------------------------
!     BUT : CREER L'OBJET NU.DLG2 PERMETTANT D'ASSOCIER LES
!           LAGRANGES 1 ET 2 D'UN NUME_DDL (SUR LA BASE VOLATILE)
!
!     IN:
!        NU     (K14) : NOM D'UN NUME_DDL
!
!     INOUT : NU EST ENRICHI DE L'OBJET .DLG2
!
!        .DLG2 = V I LONG=NEQ
!        POUR IEQ1=1,NEQ :
!         .DLG2(IEQ1) = 0 => IEQ CORRESPOND A UN DDL PHYSIQUE
!         .DLG2(IEQ1) = IEQ2 => IEQ CORRESPOND A UN DDL DE LAGRANGE
!                               ET IEQ2 EST L'AUTRE LAGRANGE QUI LUI
!                               EST ASSOCIE.
!                               DANS CE CAS .DLG2(IEQ2)=IEQ1
! ----------------------------------------------------------------------
!
    integer(kind=8) :: neq, ili, nbligr, iexi, jprno, jdlg2
    integer(kind=8) :: ima, nn, n1, n2, n3, n4
    integer(kind=8) ::  j2nema, nbma
    integer(kind=8) :: ieq2, ieq3, nueq2, nueq3, nec
    character(len=19) :: ligr19
    character(len=14) :: nu14
    character(len=8) :: nogd
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: nema(:) => null()
!
!     -- ZZNSUP : NOMBRE DE NOEUDS DE LA MAILLE TARDIVE
#define zznsup(ili,ima) zi(j2nema+ima) - zi(j2nema+ima-1) - 1
!
!     -- ZZNEMA : NUMERO DES NOEUDS DE LA MAILLE TARDIVE
#define zznema(ili,ima,j) nema(zi(j2nema+ima-1)+j-1)
!
!
!
    call jemarq()
    nu14 = nu
!
    call jeexin(nu14//'.NUME.DLG2', iexi)
    if (iexi .gt. 0) goto 999
!
    call dismoi('NB_EQUA', nu14, 'NUME_DDL', repi=neq)
    call dismoi('NOM_GD', nu14, 'NUME_DDL', repk=nogd)
    call dismoi('NB_EC', nogd, 'GRANDEUR', repi=nec)
    call jelira(nu14//'.NUME.PRNO', 'NMAXOC', nbligr)
    call jeveuo(nu14//'.NUME.NUEQ', 'L', vi=nueq)
!
!
!     - ALLOCATION ET CALCUL DE .DLG2:
    call wkvect(nu14//'.NUME.DLG2', 'V V I', neq, jdlg2)
!
!
    do ili = 2, nbligr
        call jenuno(jexnum(nu14//'.NUME.LILI', ili), ligr19)
        call jeexin(ligr19//'.LGNS', iexi)
        if (iexi .eq. 0) goto 31
!
        call jeveuo(ligr19//'.NEMA', 'L', vi=nema)
        call jelira(ligr19//'.NEMA', 'NUTIOC', nbma)
        call jeveuo(jexatr(ligr19//'.NEMA', 'LONCUM'), 'L', j2nema)
!
        call jelira(jexnum(nu14//'.NUME.PRNO', ili), 'LONMAX', n4)
        if (n4 .eq. 0) goto 31
        call jeveuo(jexnum(nu14//'.NUME.PRNO', ili), 'L', jprno)
        do ima = 1, nbma
            nn = zznsup(ili, ima)
            if (nn .eq. 3) then
                n1 = zznema(ili, ima, 1)
                n2 = zznema(ili, ima, 2)
                n3 = zznema(ili, ima, 3)
                if (((n1 .gt. 0) .and. (n2 .lt. 0)) .and. (n3 .lt. 0)) then
!             L'ELEMENT IMA EST UN SEG3 DE DUALISATION (PHYS,LAG1,LAG2)
!             N2 ET N3 SONT 2 LAGRANGES APPARIES
                    ieq2 = zi(jprno-1+(-n2-1)*(nec+2)+1)
                    ieq3 = zi(jprno-1+(-n3-1)*(nec+2)+1)
                    nueq2 = nueq(ieq2)
                    nueq3 = nueq(ieq3)
                    zi(jdlg2-1+nueq2) = nueq3
                    zi(jdlg2-1+nueq3) = nueq2
                end if
            end if
        end do
31      continue
    end do
!
!
!
!
999 continue
    call jedema()
end subroutine
