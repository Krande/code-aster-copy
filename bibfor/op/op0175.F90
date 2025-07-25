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

subroutine op0175()
!     COMMANDE :  CALC_FERRAILLAGE
! ----------------------------------------------------------------------
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterc/getres.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/w175af.h"
#include "asterfort/w175ca.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: ifm, niv, n0, nuord
    integer(kind=8) :: iret, jpara, ie, nbordr, i, nuordr
    character(len=8) :: resu, model, caraElem
    character(len=16) :: crit, concep, nomcmd
    character(len=19) :: chfer1, chfer2, chefge, resu19, resuc1
    real(kind=8) :: prec
    integer(kind=8), pointer :: nume_ordre(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(resuc1, concep, nomcmd)
    call getvid(' ', 'RESULTAT', scal=resu, nbret=n0)
    resu19 = resu
!
!
!     -- CHOIX DES INSTANTS DE CALCUL :
!     ---------------------------------
    call getvr8(' ', 'PRECISION', scal=prec, nbret=ie)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=ie)
    call rsutnu(resu19, ' ', 0, '&&OP0175.NUME_ORDRE', nbordr, prec, crit, iret)
    ASSERT(iret .eq. 0)
    ASSERT(nbordr .gt. 0)
    call jeveuo('&&OP0175.NUME_ORDRE', 'L', vi=nume_ordre)
!
!
!     -- ON PREND LE MODELE POUR LE 1ER INSTANT :
!     --------------------------------------------
    nuord = nume_ordre(1)
!
    call rsadpa(resu, 'L', 1, 'MODELE', nuord, 0, sjv=jpara)
    model = zk8(jpara)
    ASSERT(model .ne. ' ')
    call getvtx(' ', 'CARA_ELEM', scal=caraElem, nbret=ie)
    ASSERT(caraElem .ne. ' ')
!
!     -- 1. ON CREE LE CHAMP DE DONNEES (CHFER1) :
!     ---------------------------------------------
    chfer1 = '&&OP0175.CHFER1'
    call w175af(model, chfer1)
    if (niv .gt. 1) then
        call imprsd('CARTE', chfer1, 6, 'CHFER1=')
    end if
!
!     -- 2. ON APPELLE L'OPTION FERRAILLAGE :
!     -------------------------------------------
    do i = 1, nbordr
        nuordr = nume_ordre(i)
        call rsexch('F', resu19, 'EFGE_ELNO', nuordr, chefge, iret)
        call rsexch(' ', resu19, 'FERR_ELEM', nuordr, chfer2, iret)
        if (resu19 .eq. resuc1) then
            if (iret .eq. 0) then
                call utmess('A', 'CALCULEL_88', si=nuordr, sk=resu19)
            end if
        end if
        call w175ca(model, caraElem, chfer1, chefge, chfer2)
        if (niv .gt. 1) then
            call imprsd('CHAMP', chfer2, 6, 'CHFER2=')
        end if
        call rsnoch(resu19, 'FERR_ELEM', nuordr)
    end do
!
    call jedema()
end subroutine
