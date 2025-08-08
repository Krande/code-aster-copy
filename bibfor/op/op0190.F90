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

subroutine op0190()
!     COMMANDE :  VERI_FERRAILLAGE
! ----------------------------------------------------------------------
    implicit none
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim3.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/imprsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
#include "asterfort/w190af.h"
#include "asterfort/w190ca.h"
!
    integer(kind=8) :: ifm, niv, n0, nuord
    integer(kind=8) :: iret, jpara, ie, nbordr, i, nuordr, iret0, iret99
    character(len=8) :: resu, model, caraElem, noma, noma2, noma3, tych, nogd
    character(len=16) :: crit, concep, nomcmd
    character(len=19) :: chmar1, chmar2, chefge, resu19, resuc1, chamfer, ligrel
    character(len=24) ::chefge0
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
    call rsutnu(resu19, ' ', 0, '&&OP0190.NUME_ORDRE', nbordr, prec, crit, iret)
    ASSERT(iret .eq. 0)
    ASSERT(nbordr .gt. 0)
    call jeveuo('&&OP0190.NUME_ORDRE', 'L', vi=nume_ordre)

!     -- ON PREND LE MODELE POUR LE 1ER INSTANT :
!     --------------------------------------------
    nuord = nume_ordre(1)
!
    call rsadpa(resu, 'L', 1, 'MODELE', nuord, 0, sjv=jpara)
    model = zk8(jpara)
    call exlim3('AFFE', 'G', model, ligrel)
    ASSERT(model .ne. ' ')
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=ie)
    ASSERT(caraElem .ne. ' ')
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=noma)

! LECTURE DU CHAMP DE FERRAILLAGE EN ENTREE

    call getvid('', 'CHAM_FERR', scal=chamfer, nbret=iret)
    call dismoi('TYPE_CHAMP', chamfer, 'CHAMP', repk=tych)
    if (tych(1:4) .ne. 'ELEM') then
        call utmess('F', 'VERIFERRAILLAGE_1')
    end if
    call dismoi('NOM_GD', chamfer, 'CHAMP', repk=nogd)
    if (nogd(1:6) .ne. 'FER2_R') then
        call utmess('F', 'VERIFERRAILLAGE_2')
    end if
    call dismoi('NOM_MAILLA', chamfer, 'CHAMP', repk=noma2)
    if (noma .ne. noma2) then
        call utmess('F', 'VERIFERRAILLAGE_3')
    end if

! LECTURE DU CHAMP DES EFFORTS DE REFERENCE EN ENTREE
    call getvid('', "CHAM_REFE", scal=chefge0, nbret=iret0)
    ! call jecreo(chefge0,'V V R')
    if (iret0 .eq. 1) then
        call dismoi('TYPE_CHAMP', chefge0, 'CHAMP', repk=tych)
        ! if (tych(1:4) .ne. 'ELNO') then
        !     call utmess('F', 'VERIFERRAILLAGE_4')
        ! end if
        call dismoi('NOM_GD', chefge0, 'CHAMP', repk=nogd)
        if (nogd(1:6) .ne. 'SIEF_R') then
            call utmess('F', 'VERIFERRAILLAGE_5')
        end if
        call dismoi('NOM_MAILLA', chefge0, 'CHAMP', repk=noma3)
        if (noma .ne. noma3) then
            call utmess('F', 'VERIFERRAILLAGE_6')
        end if
    end if

! ! INITIALISATION DE CHEFGE0 A 0 S IL EXISTE PAS
    if (iret0 .eq. 0) then
        chefge0 = '&&OP0190.EFFORTS0'
        call alchml(ligrel, 'MARG_ELEM', 'PEFFOR0', 'V', chefge0, iret99, '')
        ASSERT(iret99 .eq. 0)
    end if
    !
!     -- 1. ON CREE LE CHAMP DE DONNEES (CHMAR1) :
!     ---------------------------------------------
    chmar1 = '&&OP0190.CHMAR1'
    call w190af(model, chmar1)

    if (niv .gt. 1) then
        call imprsd('CARTE', chmar1, 6, 'CHMAR1=')
    end if
!
!     -- 2. ON APPELLE L'OPTION FERRAILLAGE :
!     -------------------------------------------
    do i = 1, nbordr
        nuordr = nume_ordre(i)
        call rsexch('F', resu19, 'EFGE_ELNO', nuordr, chefge, iret)
        call rsexch(' ', resu19, 'MARG_ELEM', nuordr, chmar2, iret)
        ! if (resu19 .eq. resuc1) then
        !     if (iret .eq. 0) then
        !         call utmess('A', 'VERIFERRAILLAGE_7', si=nuordr, sk=resu19)
        !     end if
        ! end if
        call w190ca(model, caraElem, chmar1, chefge, chamfer, chefge0, chmar2)
        if (niv .gt. 1) then
            call imprsd('CHAMP', chmar2, 6, 'chmar2=')
        end if
        call rsnoch(resu19, 'MARG_ELEM', nuordr)
    end do
!
    call jedema()
end subroutine
