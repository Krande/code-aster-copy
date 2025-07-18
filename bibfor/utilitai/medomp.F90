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

subroutine medomp(result, modele, mater, mateco, carele, nh)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medome_once.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
!
!
    character(len=8), intent(in) :: result
    character(len=8), intent(out) :: modele
    character(len=24), optional, intent(out) :: mateco, mater
    character(len=8), optional, intent(out) :: carele
    integer(kind=8), optional, intent(out) ::  nh
!
! ----------------------------------------------------------------------
!
!  OPERATEUR POST_ELEM
!
!  SAISIE ET VERIFICATION DE LA COHERENCE DES DONNEES MECANIQUES
!  DU PROBLEME
!
! ----------------------------------------------------------------------
!
! IN  RESULT : NOM DE LA SD RESULTAT
! OUT MODELE : NOM DU MODELE
! OUT MATECO : MATERIAU CODE
! OUT CARELE : CARACTERISTIQUES ELEMENTAIRES
! OUT NH     : MODE DE FOURIER
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iret
    integer(kind=8) :: nbordr, numord, inuord, numlu
    integer(kind=8) :: n1, n2, n3
    real(kind=8) :: prec
    character(len=8) :: materi, carel
    character(len=16) :: repons
    character(len=19) :: knum
    character(len=8) :: crit
    character(len=5) :: phen
    aster_logical :: lrdm, lmater, l_ther
    integer(kind=8) :: lfour
    integer(kind=8), pointer :: v_list_store(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    materi = ' '
    carel = ' '
    modele = ' '
    nbordr = 0
    numord = 0
    l_ther = ASTER_FALSE
!
    if (result(1:1) .eq. ' ') then
!
! ----- RECUPERATION DU MODELE DANS LA COMMANDE
!
        call getvid(' ', 'MODELE', scal=modele, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'POSTELEM_20')
        end if
        call dismoi('PHENOMENE', modele, 'MODELE', repk=phen)
        if (phen .eq. 'THERM') then
            l_ther = ASTER_TRUE
        end if
        call dismoi('EXI_RDM', modele, 'MODELE', repk=repons)
        lrdm = repons .eq. 'OUI'
        call dismoi('BESOIN_MATER', modele, 'MODELE', repk=repons)
        lmater = repons .eq. 'OUI'
!
! ----- RECUPERATION DU CARA_ELEM DANS LA COMMANDE
!
        if (present(carele)) then
            call getvid(' ', 'CARA_ELEM', scal=carel, nbret=n2)
            if ((n2 .eq. 0) .and. lrdm) then
                call utmess('A', 'CALCULEL3_39')
            end if
        end if
!
! ----- RECUPERATION DU CHAM_MATER DANS LA COMMANDE
!
        if (present(mateco)) then
            call getvid(' ', 'CHAM_MATER', scal=materi, nbret=n3)
            if ((n3 .eq. 0) .and. lmater) then
                call utmess('A', 'CALCULEL3_40')
            end if
        end if
!
    else
!
        call getvis(' ', 'NUME_ORDRE', scal=numlu, nbret=inuord)
!
! ----- L'UTILISATEUR N'A PAS FOURNI DE NUMERO D'ORDRE :
! ----- RECUPERATION DU PREMIER NUMERO D'ORDRE DANS LA SD RESULTAT
!
        knum = '&&MEDOMP.NUME_ORDRE'
        if (inuord .eq. 0) then
            call getvr8(' ', 'PRECISION', scal=prec, nbret=n1)
            call getvtx(' ', 'CRITERE', scal=crit, nbret=n2)
            call rsutnu(result, ' ', 0, knum, nbordr, &
                        prec, crit, iret)
            call jeveuo(knum, 'L', vi=v_list_store)
            numlu = v_list_store(1)
        end if
!
! ----- VERIFICATION DE L'UNICITE DU MODELE DANS LE RESULTAT
!
        if (inuord .eq. 0) then
            call medome_once(result, v_list_store, nbordr, &
                             model_=modele)
            call jedetr(knum)
        else
            call medome_once(result, v_list_store, nbordr, numlu, &
                             model_=modele)
        end if
!
! ----- RECUPERATION MODELE, MATERIAU ET CARA_ELEM DANS LA SD RESULTAT
!
        call rslesd(result, numlu, modele, materi, carel)
!
    end if
!
! --- CODAGE DU MATERIAU
!
    if (present(mateco)) then
        mateco = ' '
        if (materi .ne. ' ') call rcmfmc(materi, mateco, l_ther_=l_ther)
    end if
!
    if (present(mater)) then
        mater = materi
    end if
!
! --- CARA_ELEM SI NECESSAIRE
!
    if (present(carele)) then
        carele = carel
    end if
!
! --- MODE FOURIER SI NECESSAIRE
!
    if (present(nh)) then
        nh = 0
        lfour = getexm(' ', 'MODE_FOURIER')
        if (lfour .eq. 1) then
            call getvis(' ', 'MODE_FOURIER', scal=nh, nbret=n1)
        end if
    end if
!
    call jedema()
end subroutine
