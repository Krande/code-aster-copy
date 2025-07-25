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

subroutine ctresu(nomtb)
    implicit none
#include "asterf_types.h"
#include "asterfort/ctacce.h"
#include "asterfort/ctcrtb.h"
#include "asterfort/ctdata.h"
#include "asterfort/cteltb.h"
#include "asterfort/ctnotb.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
    character(len=8) :: nomtb
!     ----- OPERATEUR CREA_TABLE , MOT-CLE FACTEUR RESU   --------------
!
!        BUT : CREER UNE TABLE A PARTIR D'UN RESULTAT OU D'UN CHAMP
!
!        IN/OUT : NOMTB (K8) : NOM DE LA TABLE
!
! ----------------------------------------------------------------------
    integer(kind=8) :: nbcmp, nbno, nbma, nbval
    aster_logical :: toucmp
    character(len=1) :: tygd
    character(len=4) :: tych
    character(len=8) :: typac, sdres, noma
    character(len=16) :: nsymb
    character(len=19) :: chpgs, chpsu
    character(len=24) :: nival, nrval, niord, nkcha, nkcmp, nkvari, mesnoe, mesmai
    character(len=24) :: label
!     ------------------------------------------------------------------
    call jemarq()
!
!  -- 1.INITIALISATION
!     -----------------
    chpgs = '&&CTRESU.PT_GAUSS_S'
    chpsu = '&&CTRESU.PT_GAUSS_P'
    nival = '&&CTRESU.ACCES_IS'
    nrval = '&&CTRESU.ACCES_R8'
    nkcha = '&&CTRESU.SD_CHAM'
    niord = '&&CTRESU.ORDRE'
    nkcmp = '&&CTRESU.CMP_USER'
    nkvari = '&&CTRESU.VARI_USER'
    mesmai = '&&CTRESU.MES_MAILLES'
    mesnoe = '&&CTRESU.MES_NOEUDS'
!
!   RECUPERATIONS DES CHAMPS
    call ctacce(nsymb, typac, nbval, nival, nrval, &
                niord, nkcha, sdres)
!
!   RECUPERATION DES NOEUDS,MAILLES,COMPOSANTES
    call ctdata(mesnoe, mesmai, nkcha, tych, toucmp, &
                nkcmp, nkvari, nbcmp, chpgs, &
                chpsu, noma, nbno, nbma, nbval, tygd)
!   field or result user name
    call getvtx('RESU', 'INTITULE', iocc=1, scal=label)

!   CREATION DE LA TABLE
    call ctcrtb(nomtb, tych, sdres, nkcha, typac, &
                toucmp, nbcmp, nbval, nkcmp, nkvari)
!
!   REMPLISSAGE DE LA TABLE
    if (tych .eq. 'NOEU') then
        call ctnotb(nbno, mesnoe, noma, nbval, nkcha, &
                    nkcmp, toucmp, nbcmp, typac, &
                    nrval, sdres, nomtb, nsymb, nival, &
                    niord, label)
    else if (tych(1:2) .eq. 'EL' .or. tych .eq. 'CART') then
        call cteltb(nbma, mesmai, noma, nbval, nkcha, &
                    nkcmp, nkvari, toucmp, nbcmp, typac, &
                    nrval, sdres, nomtb, nsymb, &
                    chpgs, chpsu, tych, nival, niord, &
                    label)
    end if
!
!   NETTOYAGE
    call jedetr(chpgs)
    call jedetr(chpsu)
    call jedetr(nival)
    call jedetr(nrval)
    call jedetr(nkcha)
    call jedetr(niord)
    call jedetr(nkcmp)
    call jedetr(nkvari)
    call jedetr(mesmai)
    call jedetr(mesnoe)
!
    call jedema()
!
end subroutine
