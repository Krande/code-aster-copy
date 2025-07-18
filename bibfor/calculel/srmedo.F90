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

subroutine srmedo(modele, mate, mateco, cara, kcha, ncha, &
                  result, nuord, nbordr, base, &
                  npass, ligrel)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/exlim1.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/medom1.h"
#include "asterfort/srlima.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ncha, nuord, nbordr, npass
    character(len=1) :: base
    character(len=8) :: modele, cara, result
    character(len=19) :: kcha
    character(len=24) :: mate, ligrel, mateco
!
!     BUT: APPEL DE MEDOM1 AVEC CONSTRUCTION DU BON LIGREL
!          POUR LE CALCUL DE L'OPTION SIRO_ELEM
!
!
! OUT    : MODELE : NOM DU MODELE
! OUT    : MATE   : CHAMP MATERIAU
! OUT    : MATECO : MATERIAU CODE
! OUT    : CARA   : NOM DU CHAMP DE CARACTERISTIQUES
! IN     : KCHA   : NOM JEVEUX POUR STOCKER LES CHARGES
! OUT    : NCHA   : NOMBRE DE CHARGES
! IN     : RESULT : NOM DE LA SD RESULTAT
! IN     : NUORD  : NUMERO D'ORDRE
! IN     : BASE   : 'G' OU 'V' POUR LA CREATION DU LIGREL
! IN/OUT : NPASS  : NOMBRE DE PASSAGE DANS LA ROUTINE
! OUT    : LIGREL : NOM DU LIGREL
!
! ----------------------------------------------------------------------
!
#include "jeveux.h"
!
    integer(kind=8) :: nbmxba
    parameter(nbmxba=2)
!
    integer(kind=8) :: nbligr, i, kmod, nbmato, nbma2d
    integer(kind=8) :: iligrs, imodls, ibases, jlisma
!
    character(len=1) :: baslig
    character(len=24) :: ligr1
    character(len=24) :: mail2d, mail3d, mailto, noobj
!
    save nbligr, iligrs, imodls, ibases
! ----------------------------------------------------------------------
!
    call jemarq()
!
!   -- recuperation du modele, cara, charges a partir du resultat et du
!      numero ordre
    call medom1(modele, mate, mateco, cara, kcha, ncha, &
                result, nuord)
!
!
!   -- pour le premier passage on initialise les tableaux sauves
    if (npass .eq. 0) then
        npass = npass+1
        nbligr = 0
        call jedetr('&&SRMEDO.LIGRS    ')
        call jedetr('&&SRMEDO.MODELS   ')
        call jedetr('&&SRMEDO.BASES    ')
        call wkvect('&&SRMEDO.LIGRS    ', 'V V K24', nbordr*nbmxba, iligrs)
        call wkvect('&&SRMEDO.MODELS   ', 'V V K8', nbordr*nbmxba, imodls)
        call wkvect('&&SRMEDO.BASES    ', 'V V K8', nbordr*nbmxba, ibases)
        call jeveut('&&SRMEDO.LIGRS    ', 'L', iligrs)
        call jeveut('&&SRMEDO.MODELS   ', 'L', imodls)
        call jeveut('&&SRMEDO.BASES    ', 'L', ibases)
    end if
!
!   -- on regarde si le modele a deja ete rencontre
    kmod = indik8(zk8(imodls-1), modele, 1, nbligr+1)
    baslig = ' '
    do i = 1, nbligr
        if (zk8(imodls-1+i) .eq. modele) then
            kmod = 1
            baslig = zk8(ibases-1+i) (1:1)
        end if
    end do
!
!   --  on regarde si le ligrel a ete cree sur la meme base
!       que la base demandee
    if ((kmod .gt. 0) .and. (baslig .eq. base)) then
!
!       -- si oui, on le reprend
        ligrel = zk24(iligrs-1+nbligr)
!
    else
!       -- sinon, on cree un nouveau ligrel
        mail2d = '&&SRMEDO.MAILLE_FACE'
        mail3d = '&&SRMEDO.MAILLE_3D_SUPP'
        mailto = '&&SRMEDO.MAILLE_2D_3D'
!
!       recuperation des mailles de face et des mailles 3d support
        call srlima(modele, mail2d, mail3d, mailto, nbma2d)
        nbmato = 2*nbma2d
        call jeveuo(mailto, 'L', jlisma)
!
        noobj = '12345678.LIGR000000.LIEL'
        call gnomsd(' ', noobj, 14, 19)
        ligr1 = noobj(1:19)
        ASSERT(ligr1 .ne. ' ')
!
        call exlim1(zi(jlisma), nbmato, modele, base, ligr1)
!
        call jedetr(mail2d)
        call jedetr(mail3d)
        call jedetr(mailto)
!
        nbligr = nbligr+1
        zk24(iligrs-1+nbligr) = ligr1
        zk8(imodls-1+nbligr) = modele
        zk8(ibases-1+nbligr) = base
        ligrel = ligr1
    end if
!
    call jedema()
!
end subroutine
