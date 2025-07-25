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

subroutine w039c3(carele, modele, ifi, form, titre, aunoeud)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/carelo.h"
#include "asterfort/chpchd.h"
#include "asterfort/dismoi.h"
#include "asterfort/detrsd.h"
#include "asterfort/imprsd.h"
#include "asterfort/irceme.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/lxlgut.h"
!
    integer(kind=8) :: ifi
    character(len=8) :: carele, modele
    character(len=80) :: titre
    character(len=*) :: form
    logical :: aunoeud
!
!     BUT:
!       IMPRIMER LES REPERES LOCAUX DES ELEMENTS
! ----------------------------------------------------------------------
!     IN MODELE  : MODELE
!     IN CARELE  : CARA_ELEM
! ----------------------------------------------------------------------
!     VARIABLES LOCALES
!
    integer(kind=8) :: iret, jaux, nbCmpDyna
    character(len=1), parameter :: nomcmp(3) = ['X', 'Y', 'Z']
    character(len=8) :: typech, sdcarm, carele8
    character(len=19) :: chrel1, chrel2, chrel3, chrelno1, chrelno2, chrelno3, modelLigrel, celmod
    character(len=19) :: chrmed(3)
    character(len=64) :: nommed(3)
    character(len=85) :: titrz, messk(3)
    character(len=16) :: field_type
    aster_logical :: l3d
! ----------------------------------------------------------------------
    call jemarq()
!
    chrel1 = carele//'.REPLO_1'
    chrel2 = carele//'.REPLO_2'
    chrel3 = carele//'.REPLO_3'
!
    call carelo(modele, carele, 'V', chrel1, chrel2, chrel3)
!
    call jeexin(chrel3//'.CELD', iret)
!   il n'y a que deux vecteurs dans le cas 2d
    if (iret .ne. 0) then
        l3d = .true.
    else
        l3d = .false.
    end if
!   IMPRESSION DES CHAMPS DE VECTEURS
!
    if (aunoeud .and. (form .eq. 'MED')) then
!       Passage des champs ELEM en ELNO
!       on complete le carele par des '_'
        jaux = lxlgut(carele)
        carele8 = '________'
        carele8(1:jaux) = carele(1:jaux)
!
        chrelno1 = carele8//'.REPLC_1'
        chrelno2 = carele8//'.REPLC_2'
        chrelno3 = carele8//'.REPLC_3'
!
!       Récupération du LIGREL
        call dismoi('NOM_LIGREL', modele, 'MODELE', repk=modelLigrel)
!
        celmod = '&&W039C3.CELMOD'
        call alchml(modelLigrel, 'TOU_INI_ELNO', 'PGEOM_R', 'V', celmod, iret, ' ')
        if (iret .ne. 0) then
            messk(1) = modelLigrel
            messk(2) = 'PGEOM_R'
            messk(3) = 'TOU_INI_ELNO'
            call utmess('F', 'UTILITAI3_23', nk=3, valk=messk)
        end if
!
        call chpchd(chrel1, 'ELNO', celmod, 'OUI', 'V', chrelno1, modele)
        call chpchd(chrel2, 'ELNO', celmod, 'OUI', 'V', chrelno2, modele)
        if (l3d) then
            call chpchd(chrel3, 'ELNO', celmod, 'OUI', 'V', chrelno3, modele)
        end if
        call detrsd('CHAMP', celmod)
!
        nommed(1) = 'RepLocal_X'
        nommed(2) = 'RepLocal_Y'
        nommed(3) = 'RepLocal_Z'
        chrmed(1) = chrelno1
        chrmed(2) = chrelno2
        chrmed(3) = chrelno3
        sdcarm = ' '
        typech = 'ELNO'
        modele = ' '
    else
        nommed(1) = chrel1
        nommed(2) = chrel2
        nommed(3) = chrel3
        chrmed(1) = chrel1
        chrmed(2) = chrel2
        chrmed(3) = chrel3
        sdcarm = ' '
        typech = 'ELEM'
    end if
!
    field_type = 'Unknown'
    if (form .eq. 'MED') then
!     -------------------------
        call irceme(ifi, nommed(1), chrmed(1), typech, modele, 0, nomcmp, ' ', ' ', 0, &
                    0.d0, 0, 0, [0], sdcarm, sdcarm, field_type, nbCmpDyna, .false._1, iret)
        ASSERT(iret .eq. 0)
!
        call irceme(ifi, nommed(2), chrmed(2), typech, modele, 0, nomcmp, ' ', ' ', 0, &
                    0.d0, 0, 0, [0], sdcarm, sdcarm, field_type, nbCmpDyna, .false._1, iret)
        ASSERT(iret .eq. 0)
!
        if (l3d) then
            call irceme(ifi, nommed(3), chrmed(3), typech, modele, 0, nomcmp, ' ', ' ', 0, &
                        0.d0, 0, 0, [0], sdcarm, sdcarm, field_type, nbCmpDyna, .false._1, iret)
            ASSERT(iret .eq. 0)
        end if
!
    else if (form .eq. 'RESULTAT') then
!     ---------------------------
        titrz = '1er '//titre
        call imprsd('CHAMP', chrel1, ifi, titrz)
        titrz = '2eme '//titre
        call imprsd('CHAMP', chrel2, ifi, titrz)
        if (l3d) then
            titrz = '3eme '//titre
            call imprsd('CHAMP', chrel3, ifi, titrz)
        end if
    else
        ASSERT(.false.)
    end if
!
    call detrsd('CHAM_ELEM', chrel1)
    call detrsd('CHAM_ELEM', chrel2)
    if (l3d) call detrsd('CHAM_ELEM', chrel3)
!
    if (aunoeud .and. (form .eq. 'MED')) then
        call detrsd('CHAM_ELEM', chrelno1)
        call detrsd('CHAM_ELEM', chrelno2)
        if (l3d) call detrsd('CHAM_ELEM', chrelno3)
    end if
!
    call jedema()
end subroutine
