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
!
subroutine singue(cherrs, chenes, nomail, ndim, nnoem, &
                  nelem, xy, prec, ligrmo, chelem, &
                  types)
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterfort/cesexi.h"
#include "asterfort/dsingu.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsingu.h"
#include "asterfort/ssingu.h"
#include "asterfort/tsingu.h"
!
    integer(kind=8) :: ndim, nnoem, nelem
    real(kind=8) :: xy(3, nnoem), prec
    character(len=8) :: nomail
    character(len=16) :: types
    character(len=19) :: cherrs, chenes
    character(len=24) :: ligrmo, chelem
!
!     BUT:
!         1) RECUPERATION DE L ERREUR ET DE L ENERGIE EN CHAQUE EF
!         2) CALCUL DU DEGRE DE LA SINGULARITE
!         3) CALCUL DU RAPPORT ENTRE L ANCIENNE ET LA NOUVELLE TAILLE
!         4) CALCUL DE LA NOUVELLE TAILLE DES EF
!            RQ : CES TROIS QUANTITES SONT CALCULEES DANS CHAQUE ELEMENT
!                 ET SONT CONSTANTES PAR ELEMENT
!         5) STOCKAGE DE CES DEUX COMPOSANTES DANS CHELEM
!         OPTION : 'SING_ELEM'
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   CHERRS      : NOM SD_S OU EST STOCKE L ERREUR
! IN   CHENES      : NOM SD_S OU EST STOCKE L ENERGIE
! IN   NOMAIL      : NOM DU MAILLAGE
! IN   NDIM        : DIMENSION DU PROBLEME
! IN   NNOEM       : NOMBRE DE NOEUDS DU MAILLAGE
! IN   NELEM       : NOMBRE D ELEMENTS FINIS DU MAILLAGE
! IN   XY(3,NNOEM) : COORDONNEES DES NOEUDS
! IN   PREC        : % DE L ERREUR TOTALE SOUHAITE
!                   POUR CALCULER LA NOUVELLE CARTE DE TAILLE
!                   DES EF H*
!                   ERREUR_TOTALE(H*)=PREC*ERREUR_TOTALE
! IN   LIGRMO      : NOM DU LIGREL DU MODELE
! IN   CHELEM      : CHAM_ELEM QUI VA CONTENIR LE DEGRE ET LA TAILLE
! IN   TYPES       : TYPE DE L ESTIMATEUR D ERREUR (NOM DE L OPTION)
!
!      SORTIE :
!-------------
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: jcesc, jcesd, jcesl, jcesv, iad
    integer(kind=8) :: nsommx, nelcom, degre
    integer(kind=8) :: nbcmp, ncmp
    integer(kind=8) :: icmp, inel, nbr(nelem), nalpha
    real(kind=8) :: erreur(nelem), taille(nelem), energi(nelem)
    real(kind=8) :: alpha(nelem), re(nelem), he(nelem)
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: conn(:) => null()
    real(kind=8), pointer :: mesu(:) => null()
    integer(kind=8), pointer :: cinv(:) => null()
!
    call jemarq()
!
! 1 - RECUPERATION DES ADRESSES DES OBJETS CREES DANS SINGUM
!
    call jeveuo('&&SINGUM.DIME           ', 'L', vi=dime)
    call jeveuo('&&SINGUM.MESU           ', 'L', vr=mesu)
    call jeveuo('&&SINGUM.CONN           ', 'L', vi=conn)
    call jeveuo('&&SINGUM.CINV           ', 'L', vi=cinv)
!
! 2 - NSOMMX = NBRE MAX DE NOEUDS SOMMETS CONNECTES AUX EF
!     NELCOM = NBRE MAX D EFS SURF EN 2D OU VOL EN 3D
!              CONNECTES AUX NOEUDS
!     DEGRE  = 1 EF LINEAIRE - 2 EF QUADRATIQUE
!
    nsommx = dime(1)
    nelcom = dime(2)
    degre = dime(3)
!
! 3 - RECUPERATION DE L'ERREUR EN CHAQUE EF ERREUR(EF)
!       ET DE LA TAILLE EN CHAQUE EF TAILLE(EF)
!     NOMBRE DE COMPOSANTES A STOCKER PAR EF NBR(NELEM)
!       2 SI EF SURFACIQUES EN 2D OU VOLUMIQUES EN 3D
!       0 SINON
!
    call jeveuo(cherrs//'.CESC', 'L', jcesc)
    call jelira(cherrs//'.CESC', 'LONMAX', nbcmp)
    call jeveuo(cherrs//'.CESD', 'L', jcesd)
    call jeveuo(cherrs//'.CESL', 'L', jcesl)
    call jeveuo(cherrs//'.CESV', 'L', jcesv)
!
! LECTURE DE LA SD .CESC
    do icmp = 1, nbcmp
        if (zk8(jcesc+icmp-1) (1:6) .eq. 'ERREST') ncmp = icmp
    end do
!
    do inel = 1, nelem
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp, iad)
        if (iad .gt. 0) then
            erreur(inel) = zr(jcesv+iad-1)
            nbr(inel) = 3
        else
            erreur(inel) = 0.d0
            nbr(inel) = 0
        end if
    end do
!
! LECTURE DE LA SD .CESC
    do icmp = 1, nbcmp
        if (zk8(jcesc+icmp-1) (1:6) .eq. 'TAILLE') ncmp = icmp
    end do
!
    do inel = 1, nelem
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp, iad)
        if (iad .gt. 0) then
            taille(inel) = zr(jcesv+iad-1)
            nbr(inel) = 3
        else
            taille(inel) = 0.d0
            nbr(inel) = 0
        end if
    end do
!
! 4 - RECUPERATION DE L'ENERGIE EN CHAQUE EF ENERGI(EF)
!
    call jeveuo(chenes//'.CESC', 'L', jcesc)
    call jelira(chenes//'.CESC', 'LONMAX', nbcmp)
    call jeveuo(chenes//'.CESD', 'L', jcesd)
    call jeveuo(chenes//'.CESL', 'L', jcesl)
    call jeveuo(chenes//'.CESV', 'L', jcesv)
!
! LECTURE DE LA SD .CESC
    do icmp = 1, nbcmp
        if (zk8(jcesc+icmp-1) (1:6) .eq. 'TOTALE') ncmp = icmp
    end do
!
    do inel = 1, nelem
        call cesexi('C', jcesd, jcesl, inel, 1, &
                    1, ncmp, iad)
        if (iad .gt. 0) then
            energi(inel) = zr(jcesv+iad-1)
        else
            energi(inel) = 0.d0
        end if
    end do
!
! 5 - CALCUL DU DEGRE DE LA SINGULARITE ALPHA(NELEM) PAR EF
!
    call dsingu(ndim, nelem, nnoem, nsommx, nelcom, &
                degre, conn, cinv, xy, erreur, &
                energi, mesu, alpha, nalpha)
!
! 6 - CALCUL DU RAPPORT DE TAILLE DES EF RE=HE*/HE
!     HE TAILLE DE L EF ACTUEL - HE* TAILLE DU NOUVEL EF
!
    call rsingu(ndim, nelem, nbr, nalpha, degre, &
                prec, erreur, alpha, types, re)
!
! 7 - CALCUL DE LA NOUVELLE TAILLE DES EF HE*=RE*HE
!     HE TAILLE DE L EF ACTUEL - HE* TAILLE DU NOUVEL EF
!
    call tsingu(nelem, nbr, re, taille, he)
!
! 8 - STOCKAGE DE ALPHA ET RE DANS CHELEM
!
    call ssingu(nomail, nelem, nbr, ligrmo, alpha, &
                re, he, chelem)
!
    call jedema()
!
end subroutine
