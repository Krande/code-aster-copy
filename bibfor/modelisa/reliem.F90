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

subroutine reliem(mo, ma, typem, motfaz, iocc, &
                  nbmocl, limocl, tymocl, litroz, nbtrou, l_keep_propz, l_allz)
    implicit none
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/reliem1.h"
!
    integer(kind=8) :: iocc, nbmocl, nbtrou
    character(len=8) :: ma
    character(len=*) :: limocl(nbmocl), tymocl(nbmocl), mo
    character(len=*) :: litroz, typem, motfaz
    aster_logical, optional, intent(in) :: l_keep_propz
    aster_logical, optional, intent(in) :: l_allz
! ----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
!
!     CE MODULE PERMET DE CREER UN OBJET JEVEUX CONTENANT UNE LISTE
!     DE NOMS OU NUMEROS DE MAILLES OU DE NOEUDS CORRESPONDANT AUX
!     MOTS-CLES TRANSMIS EN ARGUMENTS.
!
! IN  : MO     : NOM DU MODELE (FACULTATIF SINON : ' ')
!           SI LE NOM DU MODELE EST DONNE, ON VERIFIERA QUE LES MAILLES
!           (OU LES NOEUDS) RECUPERES FONT PARTIE DU MODELE.
!           S'ILS NE FONT PAS PARTIE DU LGREL => ALARME
! IN  : MA     : NOM DU MAILLAGE
! IN  : TYPEM  : PRECISE LE TYPE DE LISTE QUE L'ON VEUT RECUPERER
!              : 'NU_MAILLE'  : NUMEROS DE MAILLES
!              : 'NO_MAILLE'  : NOMS    DE MAILLES
!              : 'NU_NOEUD'   : NUMEROS DE NOEUDS
!              : 'NO_NOEUD'   : NOMS    DE NOEUDS
! IN  : MOTFAZ : NOM DU MOT CLE FACTEUR (OU ' ')
! IN  : IOCC   : NUMERO DE L'OCCURENCE DU MOT CLE FACTEUR
! IN  : NBMOCL : NOMBRE DE MOTS CLES A SCRUTER
!                (DIMENSION DE LIMOCL)
! IN  : LIMOCL : LISTE DES MOTS CLE A SCRUTER
! IN  : TYMOCL : LISTE DES TYPES DE MOTS CLE A SCRUTER :
!                / 'GROUP_MA'
!                / 'GROUP_NO'
!                / 'MAILLE'
!                / 'NOEUD'
!                / 'TOUT'   % TOUT:'OUI'
! IN/JXOUT : LITROZ : NOM DE L'OBJET JEVEUX QUI CONTIENDRA LA LISTE DES
!                     ENTITES (MAILLE OU NOEUD) TROUVEES
! OUT : NBTROU : NOMBRE D'ENTITES TROUVEES
! IN, OPTIONAL : L_KEEP_PROP : PRIS EN COMPTE UNIQUEMENT UN MAILLAGE PARALLELE
!    (CELA NE CHANGE RIEN DANS LES AUTRES CAS)
!    POUR UN PARALLEL_MESH, SI TRUE ON NE GARDE QUE LES MAILLES/NOEUDS DONT LE SOUS-DOMAINE
!    EST PROPRIETAIRE SI FALSE ON GARDE TOUT
!    SI L'ARGUMENT N'EST PAS PRESENT ON GARDE TOUT (=FALSE).
! IN, OPTIONAL     : L_ALL : TRUE  : forcer comme TOUT='OUI'
!                            FALSE : par défaut, cas normal
! ----------------------------------------------------------------------
!
    character(len=19) :: ligrel
    aster_logical :: l_keep_prop, l_all
!
    ligrel = ' '
    if (mo .ne. ' ') then
        call dismoi("NOM_LIGREL", mo, "MODELE", repk=ligrel)
    end if
    if (present(l_keep_propz)) then
        l_keep_prop = l_keep_propz
    else
        l_keep_prop = ASTER_FALSE
    end if
    if (present(l_allz)) then
        l_all = l_allz
    else
        l_all = ASTER_FALSE
    end if
    call reliem1(ligrel, ma, typem, motfaz, iocc, &
                 nbmocl, limocl, tymocl, litroz, nbtrou, l_keep_prop, l_all)

!
end subroutine
