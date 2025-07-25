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

subroutine resth2(modele, ligrel, lchar, nchar, ma, &
                  cartef, nomgdf, carteh, nomgdh, cartet, &
                  nomgdt, cartes, nomgds, chgeom, chsour, &
                  psourc)
! person_in_charge: olivier.boiteau at edf.fr
!-----------------------------------------------------------------------
!    - FONCTION REALISEE:  PREPARATION DU CALCUL DE L'ESTIMATEUR
!                          D'ERREUR EN RESIDU SUR LE PROBLEME THERMIQUE.
!
! IN MODELE  : NOM DU MODELE
! IN LIGREL  : NOM DU LIGREL
! IN LCHAR   : LISTE DES CHARGES
! IN NCHAR   : NOMBRE DE CHARGES
! OUT MA     : NOM DU MAILLAGE
! OUT CARTEF/NOMGDF: INFO SUR LE FLUX RETENU
! OUT CARTEH/NOMGDH: INFO SUR L'ECHANGE RETENU
! OUT CARTET/NOMGDT: INFO SUR LA TEMP_EXT RETENUE
! OUT CARTES/NOMGDS: INFO SUR LA SOURCE RETENUE
! OUT CHGEOM : CHAMP GEOMETRIE
! OUT CHSOUR : CHAMP SOURCE
! OUT PSOURC : NOM DU PARAMETRE ASSOCIE A SOURCE
!
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES:
!       MESSAGE:UTMESS,UTMESG.
!       JEVEUX:JEMARQ,JEDEMA,JEDETR,MEGEOM,DISMOI,EXISD,ETENCA.
!       ELEMENTS FINIS:CALCUL,RESVOI.
!
!     FONCTIONS INTRINSEQUES:
!       AUCUNE.
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       22/08/02 (OB): CREATION DU A LA SEPARATION DE RESTHE, EN UNE
!          PARTIE PRELIMINAIRE (RESTH2) HORS DE LA BOUCLE EN TEMPS ET
!          UNE PARTIE (RESTHE) CALCUL DEPENDANT DU TEMPS.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterfort/assert.h"
#include "asterfort/alchml.h"
#include "asterfort/dismoi.h"
#include "asterfort/etenca.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/megeom.h"
#include "asterfort/resvoi.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nchar
    character(len=8) :: modele, lchar(1), ma, psourc
    character(len=19) :: cartef, carteh, cartet, cartes, nomgdf, nomgdh, nomgdt
    character(len=19) :: nomgds, cartesc
    character(len=24) :: ligrel, chgeom, chsour
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: i, ier, iret, iretf, ireth, irett, irets, iretep, iretsc
    character(len=1) :: base
    character(len=19) :: cartf, carth, cartt, carts, cartep, cartsc
!
! DEBUT DE LA SUBROUTINE
    call jemarq()
!
    ASSERT(ligrel(1:8) .eq. modele)
!
! RECHERCHE DU NOM DU CHAMP GEOMETRIE DANS LA SD MODELE OU CHARGE -----
! SURCOUCHE DE LA REQUETE D'EXISTENCE JEEXIN/JEVEUO DU DESCRIPTEUR DE
! MODELE//'MODELE.NOMA.COORDO': RESULTAT DANS CHGEOM. ON S'ASSURE DE
! LA COHERENCE AVEC CHARGE//'CHTH.MODEL.NOMO'.
    base = 'V'
    call megeom(modele, chgeom)
!
!   RECHERCHE DES ELTS FINIS CONTIGUS ET REMPLISSAGE DU CHAMELEM DE TYPE
!   VOISIN VIA L'OPTION DE CALCUL 'INIT_MAIL_VOIS'.
    call alchml(ligrel, 'INIT_MAIL_VOIS', 'PVOISIN', base, '&&RESTHER.VOISIN', iret, ' ')
!
! SURCOUCHE DE DISMMO RENVOYANT LE NOM DU MAILLAGE (MA) VIA UN JEVEUO
! SUR MODELE//'.MODELE.NOMA'
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=ma)
! REMPLISSAGE DU CHAM_ELEM '&&RESTHER.VOISIN' PAR LES NUMEROS ET LES
! TYPES DE MAILLES VOISINES.
    call resvoi(modele, ma, '&&RESTHER.VOISIN')
!
! BOUCLE SUR LES CHARGEMENTS -----------------------------------------
!
! ATTENTION: POUR UN MEME TYPE DE CHARGEMENT, SEULE LA DERNIERE
! OCCURENCE EST CONSERVEE.
!
    nomgdf = ' '
    nomgdh = ' '
    nomgdt = ' '
    nomgds = ' '
    cartef = ' '
    carteh = ' '
    cartet = ' '
    cartes = ' '
    cartesc = ' '
    chsour = ' '
    iretf = 0
    ireth = 0
    irett = 0
    irets = 0
    iretsc = 0
    iretep = 0
! OPTION DE CALCUL PAR DEFAUT
    psourc = 'PSOURCR'
!
! BOUCLE SUR LES AFFE_CHAR_THER
    do i = 1, nchar
! INIT.
        cartf = lchar(i)//'.CHTH.FLURE'
        carth = lchar(i)//'.CHTH.COEFH'
        cartt = lchar(i)//'.CHTH.T_EXT'
        carts = lchar(i)//'.CHTH.SOURE'
        cartsc = lchar(i)//'.CHTH.SOURC'
        cartep = lchar(i)//'.CHTH.HECHP'
! DETERMINE L'EXISTENCE DES SD DE TYPE CHAMP_GD ET DE NOM CARTE.
        call exisd('CHAMP_GD', cartf, iretf)
        call exisd('CHAMP_GD', carth, ireth)
        call exisd('CHAMP_GD', cartt, irett)
        call exisd('CHAMP_GD', carts, irets)
        call exisd('CHAMP_GD', cartsc, iretsc)
        call exisd('CHAMP_GD', cartep, iretep)
        if (iretep .ne. 0) then
            call utmess('A', 'CALCULEL6_42')
        end if
        if (((ireth .eq. 0) .and. (irett .ne. 0)) .or. ((irett .eq. 0) .and. (ireth .ne. 0))) then
            ASSERT(.false.)
        end if
!
! TRAITEMENT DES CHARGEMENTS DE TYPE FLUX_REP/FLUN
        if (iretf .ne. 0) then
! SURCOUCHE DE DISMCA RENVOYANT LE NOM DE LA SD (NOMGD...) VIA UN
! JEVEUO/JENUNO SUR CARTF//'.DESC'.
            call dismoi('NOM_GD', cartf, 'CARTE', repk=nomgdf)
!
! EXTENSION DE LA CARTE CARTEF VIA CARTEF//'.PTMA' ET '.PTMS' SUR 'V'
            call etenca(cartf, ligrel, ier)
            ASSERT(ier .eq. 0)
!
! SEULE CARTE FLUN CONSERVEE (REGLE SURCHARGE USUELLE DE LA DERNIERE)
            if (cartef .ne. ' ') then
                call utmess('I', 'CALCULEL6_43', sk='FLUX LINEAIRE')
                call jedetr(cartef//'.PTMA')
                call jedetr(cartef//'.PTMS')
            end if
            cartef = cartf
        end if
!
! TRAITEMENT DES CHARGEMENTS DE TYPE ECHANGE/COEF_H
        if (ireth .ne. 0) then
            call dismoi('NOM_GD', carth, 'CARTE', repk=nomgdh)
            call etenca(carth, ligrel, ier)
            ASSERT(ier .eq. 0)
!
! TRAITEMENT DES CHARGEMENTS DE TYPE ECHANGE/TEMP_EXT
            call dismoi('NOM_GD', cartt, 'CARTE', repk=nomgdt)
            call etenca(cartt, ligrel, ier)
            ASSERT(ier .eq. 0)
!
! SEULE CARTE FLUN CONSERVEE (REGLE SURCHARGE USUELLE DE LA DERNIERE)
            if (carteh .ne. ' ') then
                call utmess('I', 'CALCULEL6_43', sk='ECHANGE')
                call jedetr(carteh//'.PTMA')
                call jedetr(carteh//'.PTMS')
                call jedetr(cartet//'.PTMA')
                call jedetr(cartet//'.PTMS')
            end if
            carteh = carth
            cartet = cartt
        end if
!
! TRAITEMENT DES SOURCES VOLUMIQUES
        if (irets .ne. 0) then
            chsour = carts//'.DESC'
            call dismoi('NOM_GD', carts, 'CARTE', repk=nomgds)
! SEULE CARTE FLUN CONSERVEE (REGLE SURCHARGE USUELLE DE LA DERNIERE)
            if (cartes .ne. ' ') then
                call utmess('A', 'CALCULEL6_43', sk='SOURCE')
            end if
! OPTION DE CALCUL POUR SOURCE VARIABLE
            if (nomgds(1:6) .eq. 'SOUR_F') psourc = 'PSOURCF'
            cartes = carts
        end if
!
! TRAITEMENT DES SOURCES VOLUMIQUES CALCULEES
        if (iretsc .ne. 0) then
            chsour = cartsc//'.DESC'
            call dismoi('NOM_GD', cartsc, 'CARTE', repk=nomgds)
! SEULE CARTE FLUN CONSERVEE (REGLE SURCHARGE USUELLE DE LA DERNIERE)
            if (cartesc .ne. ' ') then
                call utmess('A', 'CALCULEL6_43', sk='SOURCE')
            end if
! OPTION DE CALCUL POUR SOURCE VARIABLE INTERDITE
            if (nomgds(1:6) .eq. 'SOUR_F') then
                ASSERT(.false.)
            end if
            cartesc = carts
        end if
!
! FIN BOUCLE AFFE_CHAR_THER
    end do
!
    call jedema()
!
end subroutine
