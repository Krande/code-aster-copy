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

subroutine rvchgr(mailla, nlsnac, repere, sdnewr, &
                  iret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rvegal.h"
#include "asterfort/rvrepn.h"
#include "asterfort/utmess.h"
    character(len=24) :: nlsnac
    character(len=19) :: sdnewr
    character(len=8) :: mailla, repere
    integer(kind=8) :: iret
!
!
!  OPERATION REALISEE
!  ------------------
!
!     CALCUL DU REPERE DE TRAVAIL (LOCALE, DIRECTION OU POLAIRE)
!
!  ARGUMENTS EN ENTREE
!  -------------------
!
!     MAILLA : NOM DU CONCEPT MAILLAGE
!     NLSNAC : NOM DU VECTEUR DES NOEUDS LIEU DE POST-TRAITEMENT
!     REPERE : VAUT 'LOCAL' OU 'POLAIRE'
!
!  ARGUMENTS EN SORTIE
!  -------------------
!
!     IRET   : CODE RETOUR : 1 RAS, 0 ERREUR (EMISSION D' UN MESSAGE)
!     SDNEWR : NOM DE LA SD CONSERVANT LE NOUVEAU REPERE
!
!              .VEC1 : XD V R8 -->   COORD. DU VECTEUR 1 DANS (X,Y)
!              .VEC2 : XD V R8 -->   COORD. DU VECTEUR 2 DANS (X,Y)
!
!              VECJ(2I-1) <--  X(VECT_J,POINT_I)
!              VECJ(2I  ) <--  Y(VECT_J,POINT_I)
!
!              NB_OC = NB_PARTIE DU LIEU
!
!***********************************************************************
!
!  -----------------------------------------
!
!
!
!  ---------------------------------
!
!  VARIABLES LOCALES
!  -----------------
!
    integer(kind=8) :: i, nd, nbnac, ind, alsnac, ierd
    aster_logical :: egal
    real(kind=8) :: znd, zref, aux
    character(len=8) :: k8b
    real(kind=8), pointer :: vale(:) => null()
!
!====================== CORPS DE LA ROUTINE ===========================
!
    call jemarq()
    i = 0
    iret = 1
!
    if (repere(1:7) .eq. 'POLAIRE') then
        call dismoi('Z_CST', mailla, 'MAILLAGE', repk=k8b, arret='C', &
                    ier=ierd)
        if (k8b(1:3) .eq. 'NON') then
            iret = 0
            call utmess('A', 'POSTRELE_28')
            goto 999
        end if
    end if
!
    ind = 1
!
    call jelira(nlsnac, 'LONMAX', nbnac)
    call jeveuo(nlsnac, 'L', alsnac)
    call jeveuo(mailla//'.COORDO    .VALE', 'L', vr=vale)
!
    nd = zi(alsnac+1-1)
!
    zref = vale(3)
!
10  continue
    if ((iret .ne. 0) .and. (ind .le. nbnac)) then
!
        nd = zi(alsnac+ind-1)
        znd = vale(1+3*nd-1)
!
        call rvegal(1.0d-3, 'R', zref, znd, egal, &
                    aux)
!
        if (.not. egal) then
!
            iret = 0
!
        end if
!
        ind = ind+1
!
        goto 10
!
    end if
!
    if (iret .ne. 0) then
!
        call rvrepn(mailla, nlsnac, repere, sdnewr)
!
    else
!
        call utmess('A', 'POSTRELE_28')
!
    end if
!
999 continue
!
    call jedema()
end subroutine
