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

subroutine cglect(resu, modele, ndim, option, &
                  typfis, nomfis, fonoeu, chfond, basfon, &
                  taillr, conf, lnoff, liss, ndeg, typdis)
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/cgleff.h"
#include "asterfort/cgtyfi.h"
#include "asterfort/cgvefo.h"
#include "asterfort/cgveli.h"
#include "asterfort/cgvemf.h"
#include "asterfort/cgverc.h"
#include "asterfort/cgveth.h"
#include "asterfort/cgvcmo.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndim, lnoff, ndeg
    character(len=8) :: resu, modele, typfis, nomfis, conf
    character(len=16) :: option, typdis
    character(len=24) :: fonoeu, chfond, basfon, taillr, liss
!
! person_in_charge: samuel.geniaut at edf.fr
!
!     SOUS-ROUTINE DE L'OPERATEUR CALC_G
!
!     BUT : LECTURE ET VERIFICATION DES OPERANDES
!
!  IN :
!  OUT :
!     RESU   : MOT-CLE RESULTAT
!     MODELE : MODELE ASSOCIE A RESU
!     NDIM   : DIMENSION DU MODELE
!     OPTION : MOT-CLE OPTION
!     TYPFIS : TYPE D'OBJET POUR DECRIRE LE FOND DE FISSURE
!              'FONDFISS' OU 'FISSURE' OU 'THETA'
!     NOMFIS : NOM DE L'OBJET POUR DECRIRE LE FOND DE FISSURE
!     FONOEU : NOMS DES NOEUDS DU FOND DE FISSURE
!     CHFOND : COORDONNES DES POINTS/NOEUDS DU FOND DE FISSURE
!     BASFON : BASE LOCALE AU FOND DE FISSURE
!     TAILLR : TAILLES DE MAILLES CONNECTEES AUX NOEUDS
!     CONF  : CONFIGURATION DE LA FISSURE EN FEM
!
!     LNOFF  : NOMBRE DE NOEUDS (OU POINTS) DU FOND DE FISSURE
!     LISS   : TYPE DE LISSAGE (NOM UNIQUE CONTRACTE)
!     TYPDIS : TYPE DE DISCONTINUITE SI FISSURE XFEM
!              'FISSURE' OU 'COHESIF'
! ======================================================================
!
    integer(kind=8) :: ier, nexci
!
    call jemarq()
!
!     RECUPERATION DE LA SD RESULTAT : RESU
    call getvid(' ', 'RESULTAT', scal=resu, nbret=ier)
!
!     RECUPERATION DE L'OPTION
    call getvtx(' ', 'OPTION', scal=option, nbret=ier)
!
!     DETERMINATION DU TYPFIS = 'FONDFISS' OU 'FISSURE' OU 'THETA'
!     ET RECUPERATION DE LA SD POUR DECRIRE LE FOND DE FISSURE : NOMFIS
!     TYPE DE DISCONTINUITE SI FISSURE XFEM: 'FISSURE' OU 'COHESIF'
    call cgtyfi(typfis, nomfis, typdis)
!
!     LECTURE DES CHARGES ET VERIFICATION DE LA COMPATIBILITE AVEC RESU
    if (typdis .ne. 'COHESIF') then
        call getfac('EXCIT', nexci)
        call cgverc(resu, nexci)
    end if
!
!     RECUPERATION DU MODELE PUIS DE LA DIMENSION DU MODELE
    call dismoi('MODELE', resu, 'RESULTAT', repk=modele)
    call dismoi('DIM_GEOM', modele, 'MODELE', repi=ndim)
!
!   CALCUL COHESIF OUVERT EN 3D UNIQUEMENT POUR L INSTANT
    if (ndim .eq. 2 .and. typdis .eq. 'COHESIF') then
        call utmess('F', 'RUPTURE2_5')
    end if
!
!     VERIFICATION DE LA COMPATIBILITE ENTRE LA SD ASSOCIEE AU FOND
!     DE FISSURE ET LE MODELE
    call cgvemf(modele, typfis, nomfis, typdis)

!     VERIFICATION DE LA COMPATIBILITE ENTRE LA DIMENSION DU MODELE
!     ET DE CELLES DES ELEMENTS EN FOND DE FISSURE
    call cgvcmo(modele, nomfis, typfis, ndim)
!
!     VERIFICATION DE LA COMPATIBILITE ENTRE OPTION ET TYPE DE FISSURE
    call cgvefo(option, typfis, nomfis, typdis)
!
!     VERIFICATION DES DONNEES RELATIVES AU(X) CHAMP(S) THETA
    call cgveth(typfis, ndim)
!
!     LECTURE DE LA DESCRIPTION DU FOND DE FISSURE
!     ET RECUPERATION DES OBJETS FONOEU, CHFOND, BASFON + LNOFF
    call cgleff(typfis, nomfis, fonoeu, chfond, basfon, &
                taillr, conf, lnoff)
!
!     VERIFICATION DES DONNEES RELATIVES AU LISSAGE
!     ET DETERMINATION DU LISSAGE (NOM UNIQUE CONTRACTE) : LISS ET NDEG
    call cgveli(typfis, typdis, ndim, lnoff, liss, &
                ndeg)
!
    call jedema()
!
end subroutine
