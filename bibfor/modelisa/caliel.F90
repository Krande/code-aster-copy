! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine caliel(fonrez, chargz)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/aflrch.h"
#include "asterfort/caarle.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/nueffe.h"
#include "asterfort/rapo2d.h"
#include "asterfort/rapo3d.h"
#include "asterfort/rapoco.h"
#include "asterfort/wkvect.h"
    character(len=*) :: chargz, fonrez
! -------------------------------------------------------
!     MODELISATION DU RACCORD ENTRE DES ELEMENTS
!     AYANT DES MODELISATIONS DIFFERENTES PAR DES RELATIONS
!     LINEAIRES ENTRE DDLS.
!     CES RELATIONS SONT AFFECTEES A LA CHARGE CHARGZ.
!     TYPES DES RACCORDS TRAITES :
!       1) RACCORD POUTRE-3D PAR DES RELATIONS LINEAIRES
!          ENTRE LES NOEUDS DES MAILLES DE SURFACE MODELISANT
!          LA TRACE DE LA SECTION DE LA POUTRE SUR LE MASSIF 3D
!          ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
!
!       2) RACCORD POUTRE-COQUE PAR DES RELATIONS LINEAIRES
!          ENTRE LES NOEUDS DES MAILLES DE BORD DE COQUE MODELISANT
!          LA TRACE DE LA SECTION DE LA POUTRE SUR A COQUE
!          ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
! -------------------------------------------------------
!  FONREZ        - IN    - K4   - : 'REEL' OU 'FONC'
!  CHARGZ        - IN    - K8   - : NOM DE LA SD CHARGE
!                - JXVAR -      -
! -------------------------------------------------------
!
! -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       19/03/04 (OB): PAR ADHERENCE A NUEFFE
!--------------------------------------------------------------------
!
!
! --------- VARIABLES LOCALES ---------------------------
!
    character(len=8) :: mod, charge
    character(len=14) :: numddl
    character(len=16) :: motfac, option
    character(len=19) :: ligrmo, lisrel
    integer :: iocc, nliai, iop, nb_ligr
    character(len=24), pointer :: list_ligr(:) => null()
!
! --------- FIN  DECLARATIONS  VARIABLES LOCALES --------
!
    call jemarq()
    charge = chargz
    motfac = 'LIAISON_ELEM'
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 999
!
! --- NOM DE LA LISTE DE RELATIONS
!
    lisrel = '&&CALIEL.RLLISTE'
!
! --- MODELE ASSOCIE AU LIGREL DE CHARGE
!     ----------------------------------
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', repk=mod)
!
! ---  LIGREL DU MODELE
!
    ligrmo = mod(1:8)//'.MODELE'
!
! --- CREATION SUR LA VOLATILE DU NUMEDDL ASSOCIE AU LIGREL
! --- DU MODELE
!     -----------------------------------------------------
    nb_ligr = 1
    call wkvect('&&CALIEL.LIGRMO', 'V V K24', 1, vk24=list_ligr)
    list_ligr(1) = ligrmo

    numddl = '&&CALIEL.NUMED'
    call nueffe(nb_ligr, list_ligr, 'VV', numddl, 'SANS', mod)
!
    do iocc = 1, nliai
        call getvtx(motfac, 'OPTION', iocc=iocc, scal=option, nbret=iop)
        if (option .eq. '3D_POU') then
            call rapo3d(numddl, iocc, fonrez, lisrel, chargz)
        else if (option .eq. '3D_POU_ARLEQUIN') then
            call caarle(numddl, iocc, lisrel, chargz)
        else if (option .eq. '2D_POU') then
            call rapo2d(numddl, iocc, fonrez, lisrel, chargz)
        else if (option .eq. '3D_TUYAU') then
            call rapo3d(numddl, iocc, fonrez, lisrel, chargz)
        else if (option .eq. 'PLAQ_POUT_ORTH') then
            call rapo3d(numddl, iocc, fonrez, lisrel, chargz)
        else if (option .eq. 'COQ_POU') then
            call rapoco(numddl, iocc, fonrez, lisrel, chargz)
        else if (option .eq. 'COQ_TUYAU') then
            call rapoco(numddl, iocc, fonrez, lisrel, chargz)
        end if
    end do
!
!     -- AFFECTATION DE LA LISTE_RELA A LA CHARGE :
!     ---------------------------------------------
    call aflrch(lisrel, charge, 'NLIN')
!
!
! --- MENAGE
!
    call jedetc('V', '&&CALIEL.RLLISTE', 1)
    call jedetr('&&CALIEL.LIGRMO')
    call jedetr('&&CALIEL.NUMED')
!
999 continue
    call jedema()
end subroutine
