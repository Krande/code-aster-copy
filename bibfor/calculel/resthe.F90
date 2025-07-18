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

subroutine resthe(ligrel, evol, chtemm, chtemp, chflum, &
                  chflup, mateco, valthe, insold, inst, &
                  resu, niveau, ifm, niv, ma, &
                  cartef, nomgdf, carteh, nomgdh, cartet, &
                  nomgdt, cartes, nomgds, chgeom, chsour, &
                  psourc, iaux)
! person_in_charge: olivier.boiteau at edf.fr
!-----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DE L'ESTIMATEUR D'ERREUR EN RESIDU
!                          SUR LE PROBLEME THERMIQUE
!
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
! IN LIGREL : NOM DU LIGREL
! IN EVOL   : LOGICAL=TRUE SI CALCUL TRANSITOIRE, FALSE SINON
! IN CHTEMM : NOM DU CHAMP DE TEMPERATURE A T-
! IN CHTEMP : NOM DU CHAMP DE TEMPERATURE A T+
! IN CHFLUM : NOM DU CHAMP DU FLUX DE TEMPERATURE A T-
! IN CHFLUP : NOM DU CHAMP DU FLUX DE TEMPERATURE A T+
! IN MATECO : NOM DU CONCEPT MATERIAU CODE
! IN VALTHE : VALEUR DU PARAMETRE THETA
! IN INSOLD : INSTANT CORRESPONDANT AU CALCUL PRECEDENT
! IN INST   : INSTANT CORRESPONDANT AU CALCUL ACTUEL
! IN NIVEAU : NIVEAU DE CALCUL DE L'ESTIMATEUR
! IN IFM/NIV: NIVEAU IMPRESSION
! IN MA     : NOM DU MAILLAGE
! IN CARTEF/NOMGDF: INFO SUR LE FLUX RETENU
! IN CARTEH/NOMGDH: INFO SUR L'ECHANGE RETENU
! IN CARTET/NOMGDT: INFO SUR LA TEMP_EXT RETENUE
! IN CARTES/NOMGDS: INFO SUR LA SOURCE RETENUE
! IN CHGEOM : CHAMP GEOMETRIE
! IN CHSOUR : CHAMP SOURCE
! IN PSOURC : NOM DU PARAMETRE ASSOCIE A SOURCE
! IN IAUX   : NUME_ORDRE
! OUT RESU   : NOM DU CHAM_ELEM_ERREUR PRODUIT
!              SI RESU EXISTE DEJA,ON LE DETRUIT.
!
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES:
!       JEVEUX:JEMARQ,JEDEMA,JEDETR,MEGEOM,DISMOI,EXISD,ETENCA,MECACT,
!              JEVEUO,JELIRA,JEXNUM,WKVECT.
!       ELEMENTS FINIS:CALCUL,RESVOI.
!
!     FONCTIONS INTRINSEQUES:
!       AUCUNE.
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       22/06/01 (OB): CREATION EN S'INSPIRANT DE RESLOC.F.
!       22/08/02 (OB): MODIFICATION DU A SEPARATION DE RESTHE, EN UNE
!          PARTIE PRELIMINAIRE (RESTH2) HORS DE LA BOUCLE EN TEMPS ET
!          UNE PARTIE (RESTHE) CALCUL DEPENDANT DU TEMPS.
!          RAJOUT DE L'OBJET '&&RESTHE.JEVEUO' POUR AMELIORER LES
!          PERFORMANCES DE RESTHE/CALCUL/TE0003.
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1504
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mecact.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: niveau, ifm, niv, iaux
    real(kind=8) :: valthe, insold, inst
    aster_logical :: evol
    character(len=8) :: ma, psourc
    character(len=19) :: cartef, carteh, cartet, cartes, nomgdf, nomgdh, nomgdt
    character(len=19) :: nomgds
    character(len=24) :: chtemm, chtemp, chflum, chflup, resu, ligrel, mateco
    character(len=24) :: chgeom, chsour
!
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: nbcmp, nbin, nbout, iarepe, mceld, mcelv, pceld, pcelv, igd
    integer(kind=8) :: iadef, iavaf, ncmpf, iadeh, iavah, ncmph, ncmpt, iadet, iavat
    integer(kind=8) :: nbjeve, ijeveo
    character(len=1) :: base
    character(len=8) :: lpain(9), lpaout(1), licmp(19)
    character(len=24) :: lchin(9), lchout(1), kcmp(19), charev
!
    call jemarq()
!
! RECHERCHE DES ADRESSES POUR OBTENIR FLUXM/P SUR LES VOISINS
    call jeveuo(ligrel(1:19)//'.REPE', 'L', iarepe)
!
    if (iaux .gt. 1) then
        call jeveuo(chflum(1:19)//'.CELD', 'L', mceld)
        call jeveuo(chflum(1:19)//'.CELV', 'L', mcelv)
    else
        mceld = 1
        mcelv = 1
    end if
    call jeveuo(chflup(1:19)//'.CELD', 'L', pceld)
    call jeveuo(chflup(1:19)//'.CELV', 'L', pcelv)
!
! RECHERCHE DES ADRESSES POUR LES CHARGES SUR LES SEGMENTS ET DE
! LEURS NOMBRES DE COMPOSANTES AU NIVEAU INFO (GRANDEUR PREMIERE)
    iadef = 1
    iavaf = 1
    ncmpf = 0
    iadeh = 1
    iavah = 1
    ncmph = 0
    iadet = 1
    iavat = 1
    ncmpt = 0
    if (cartef .ne. ' ') then
        call jeveuo(cartef//'.DESC', 'L', iadef)
        call jeveuo(cartef//'.VALE', 'L', iavaf)
        igd = zi(iadef)
        call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmpf)
    end if
    if (carteh .ne. ' ') then
        call jeveuo(carteh//'.DESC', 'L', iadeh)
        call jeveuo(carteh//'.VALE', 'L', iavah)
        igd = zi(iadeh)
        call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmph)
        call jeveuo(cartet//'.DESC', 'L', iadet)
        call jeveuo(cartet//'.VALE', 'L', iavat)
        igd = zi(iadet)
        call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmpt)
    end if
!
! STOCKAGE DE CES DONNEES DANS UN VECTEUR D'ENTIER DONT LE NOM SERA
! TRANSMIS VIA LA CARTE
    nbjeve = 14
    call wkvect('&&RESTHE.JEVEUO', ' V V I', nbjeve, ijeveo)
    zi(ijeveo) = niveau
    zi(ijeveo+1) = ifm
    zi(ijeveo+2) = niv
    zi(ijeveo+3) = iarepe
    zi(ijeveo+4) = mceld
    zi(ijeveo+5) = mcelv
    zi(ijeveo+6) = pceld
    zi(ijeveo+7) = pcelv
    zi(ijeveo+8) = iavaf
    zi(ijeveo+9) = ncmpf
    zi(ijeveo+10) = iavah
    zi(ijeveo+11) = ncmph
    zi(ijeveo+12) = iavat
    zi(ijeveo+13) = ncmpt
!
! CONSTITUTION DES CARTES REPRESENTANT LES CHARGEMENTS SURFACIQUES---
! LA REPRESENTATION DE DONNEES PAR CARTE PERMETTANT DE TRANSFERER DES
! DONNEES GLOBALES AU NIVEAU DES CALCULS ELEMENTAIRES.
! NBRE ET VECTEURS DES COMPOSANTES DE LA CARTE
! CALCULS PRELIMINAIRES
    nbcmp = 19
    if (evol) then
        charev = 'EVOL'
    else
        charev = ' '
    end if
    licmp(1) = 'Z1'
    licmp(2) = 'Z2'
    licmp(3) = 'Z3'
    licmp(4) = 'Z4'
    licmp(5) = 'Z5'
    licmp(6) = 'Z6'
    licmp(7) = 'Z7'
    licmp(8) = 'Z8'
    licmp(9) = 'Z9'
    licmp(10) = 'Z10'
    licmp(11) = 'Z11'
    licmp(12) = 'Z12'
    licmp(13) = 'Z13'
    licmp(14) = 'Z14'
    licmp(15) = 'Z15'
    licmp(16) = 'Z16'
    licmp(17) = 'Z17'
    licmp(18) = 'Z18'
    licmp(19) = 'Z19'
    kcmp(1) = ma
    kcmp(2) = ligrel
    kcmp(3) = chflum
    kcmp(4) = chflup
    kcmp(5) = cartef
    kcmp(6) = nomgdf
    kcmp(7) = carteh
    kcmp(8) = nomgdh
    kcmp(9) = cartet
    kcmp(10) = nomgdt
    kcmp(11) = cartes
    kcmp(12) = nomgds
    write (kcmp(13), '(F19.8)') valthe
    write (kcmp(14), '(F19.8)') insold
    write (kcmp(15), '(F19.8)') inst
    kcmp(16) = charev
    write (kcmp(17), '(I24)') ijeveo
! DEUX DERNIERES COMPOSANTES REDONDANTES POUR PLUS TARD
    kcmp(18) = ma
    kcmp(19) = ma
!
!
! CARTE NOMMEE &&RESTHER.CHARGE SUR LE MODELE LIGREL.
    base = 'V'
    call mecact(base, '&&RESTHER.CHARGE', 'MODELE', ligrel, 'NEUT_K24', &
                ncmp=nbcmp, lnomcmp=licmp, vk=kcmp)
!
! LANCEMENT DES CALCULS ELEMENTAIRES-----------------------------------
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PFLUX_M'
    lchin(2) = chflum
    lpain(3) = 'PFLUX_P'
    lchin(3) = chflup
    lpain(4) = psourc
    lchin(4) = chsour
    lpain(5) = 'PMATERC'
    lchin(5) = mateco
    lpain(6) = 'PCHARG'
    lchin(6) = '&&RESTHER.CHARGE'
    lpain(7) = 'PVOISIN'
    lchin(7) = '&&RESTHER.VOISIN'
    lpain(8) = 'PTEMP_M'
    lchin(8) = chtemm
    lpain(9) = 'PTEMP_P'
    lchin(9) = chtemp
    nbin = 9
!
    lpaout(1) = 'PERREUR'
    lchout(1) = resu
    nbout = 1
!
    call calcul('S', 'ERTH_ELEM', ligrel, nbin, lchin, &
                lpain, nbout, lchout, lpaout, 'G', &
                'OUI')
!
    call jedetr('&&RESTHE.JEVEUO')
    call jedema()
!
end subroutine
