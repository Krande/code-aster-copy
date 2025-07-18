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

subroutine reslo2(modele, ligrel, chvois, cvoisx, tabido)
    implicit none
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/alchml.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/megeom.h"
#include "asterfort/resvoi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: tabido(5)
    character(len=8) :: modele
    character(len=24) :: chvois
    character(len=*) :: ligrel
!
!     BUT:
!         PREPARER LE CALCUL DE L'ESTIMATEUR D'ERREUR EN RESIDU.
!         RECHERCHER LES VOISINS DE CHAQUE ELEMENT ET RECUPERER
!         DES ADRESSES.
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   MODELE : NOM DU MODELE
! IN   LIGREL : NOM DU LIGREL
!
!      SORTIE :
!-------------
! OUT  CHVOIS : NOM DU CHAMP DES VOISINS
! OUT  CVOISX : POUR XFEM, NOM DU CHAMP DES VOISINS DES SOUS-ELEMENTS
! OUT  TABIDO : TABLEAU D'ENTIERS CONTENANT DES ADRESSES
!           1 : IATYMA : ADRESSE DU VECTEUR TYPE MAILLE (NUMERO <-> NOM)
!           2 : IAGD   : ADRESSE DU VECTEUR GRANDEUR (NUMERO <-> NOM)
!           3 : IACMP  : ADRESSE DU VECTEUR NOMBRE DE COMPOSANTES
!                  (NUMERO DE GRANDEUR <-> NOMBRE DE COMPOSANTES)
!           4 : ICONX1 : ADRESSE DE LA COLLECTION CONNECTIVITE
!           5 : ICONX2 : ADRESSE DU POINTEUR DE LONGUEUR DE LA
!                        CONNECTIVITE
! ......................................................................
!
!
!
!
    character(len=6) :: nompro
    parameter(nompro='RESLO2')
!
    integer(kind=8) :: nbpin
    parameter(nbpin=2)
    integer(kind=8) :: nbpout
    parameter(nbpout=1)
!
    integer(kind=8) :: iret, nbtm, ity, nbgd, igd, ncmp
    integer(kind=8) :: iacmp, iagd, iatyma, iconx1, iconx2
    character(len=1) :: base
    character(len=8) :: lpain(nbpin), lpaout(nbpout), ma, typema, gd
    character(len=16) :: opt
    character(len=19) :: cnseto, loncha
    character(len=24) :: lchin(nbpin), lchout(nbpout), chgeom
    character(len=24) :: blan24, cvoisx
!
!
! ----------------------------------------------------------------------
!
!               123456789012345678901234
    blan24 = '                        '
!
    base = 'V'
    call megeom(modele, chgeom)
!
! ------- RECHERCHE DES VOISINS ----------------------------------------
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
!
    lpaout(1) = 'PVOISIN'
    chvois = '&&'//nompro//'.CH_VOISIN'
    lchout(1) = chvois
    opt = 'INIT_MAIL_VOIS'
    call alchml(ligrel, opt, 'PVOISIN', base, lchout(1), iret, ' ')
    call exisd('CHAMP_GD', lchout(1), iret)
    if (iret .eq. 0) then
        call utmess('F', 'CALCULEL2_88', sk=opt)
    end if
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=ma)
    call resvoi(modele, ma, chvois)
!
! ------- SI LE MODELE EST XFEM ON CALCULE LA SD VOISIN ----------------
! ---------- POUR LES SOUS-ELEMENTS DES ELEMENTS XFEM ------------------
!
    cvoisx = blan24
    call jeexin(modele//'.FISS', iret)
!
    if (iret .ne. 0) then
        opt = 'CHVOIS_XFEM'
        cnseto = modele//'.TOPOSE.CNS'
        loncha = modele//'.TOPOSE.LON'
        cvoisx = '&&'//nompro//'.CH_CNINVX'
!
        lpain(1) = 'PCNSETO'
        lchin(1) = cnseto
        lpain(2) = 'PLONCHA'
        lchin(2) = loncha
        lpaout(1) = 'PCVOISX'
        lchout(1) = cvoisx
!
        call calcul('C', opt, ligrel, 2, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')
!
        call exisd('CHAMP_GD', lchout(1), iret)
!
        if (iret .eq. 0) then
            call utmess('F', 'CALCULEL2_88', sk=opt)
        end if
!
    end if
!
! ----- CALCUL DE 5 ADRESSES : -----------------------------------------
!      IATYMA : ADRESSE DU VECTEUR TYPE MAILLE (NUMERO <-> NOM)
!      IAGD   : ADRESSE DU VECTEUR GRANDEUR (NUMERO <-> NOM)
!      IACMP : ADRESSE DU VECTEUR NOMBRE DE COMPOSANTES
!                     (NUMERO DE GRANDEUR <-> NOMBRE DE COMPOSANTES)
!      ICONX1 : ADRESSE DE LA COLLECTION CONNECTIVITE
!      ICONX2 : ADRESSE DU POINTEUR DE LONGUEUR DE LA CONNECTIVITE
!
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtm)
    call wkvect('&&'//nompro//'.TYPEMA', 'V V K8', nbtm, iatyma)
!
    do ity = 1, nbtm
        call jenuno(jexnum('&CATA.TM.NOMTM', ity), typema)
        zk8(iatyma-1+ity) = typema
    end do
!
    call jelira('&CATA.GD.NOMGD', 'NOMMAX', nbgd)
    call wkvect('&&'//nompro//'.GD', 'V V K8', nbgd, iagd)
!
    do igd = 1, nbgd
        call jenuno(jexnum('&CATA.GD.NOMGD', igd), gd)
        zk8(iagd-1+igd) = gd
    end do
!
    call wkvect('&&'//nompro//'.NBCMP', 'V V I', nbgd, iacmp)
!
    do igd = 1, nbgd
        call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', ncmp)
        zi(iacmp-1+igd) = ncmp
    end do
!
    call jeveuo(ma//'.CONNEX', 'L', iconx1)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', iconx2)
!
! ----- REMPLISSAGE DU TABLEAU -----------------------------------------
!
    tabido(1) = iatyma
    tabido(2) = iagd
    tabido(3) = iacmp
    tabido(4) = iconx1
    tabido(5) = iconx2
!
end subroutine
