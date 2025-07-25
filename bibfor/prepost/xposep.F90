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

subroutine xposep(mo, malini, mailc, mailx, nsetot, &
                  nnntot, ncotot, logrma, listgr)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ismali.h"
#include "asterfort/ligrma.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xpogma.h"
!
    character(len=8) :: mo, malini
    character(len=24) :: mailc, mailx, logrma, listgr
!
!   SEPARATION DES MAILLES DE MALINI EN 2 GROUPES
!              - MAILC : MAILLES MAILLES NON AFFECTEES D'UN MODELE
!                        OU MAILLES NON SOUS-DECOUPEES (CLASSIQUE)
!              - MAILX : MAILLES SOUS-DECOUPEES (X-FEM)
!
!   ET CALCUL DU
!        - NOMBRE TOTAL DE SOUS-ELEMENT (NSETOT)
!        - NOMBRE TOTAL DE NOUVEAUX NOEUDS (NNNTOT)
!        - LONGUEUR DE LA CONNECTIVITE DES NOUVEAUX NOEUDS (NCOTOT)
!
!   IN
!       MO     : MODELE FISSURE
!       MASSMO : TRAITEMENT DES MAILLES SANS MODELE ('OUI' OU 'NON')
!       MALINI : MAILLAGE SAIN
!
!   OUT
!       MAILC  : LISTE DES NUMEROS DES MAILLES NON AFFECTEES D'UN MODELE
!                OU NON SOUS-DECOUPEES
!       MAILX  : LISTE DES NUMEROS DES MAILLES SOUS-DECOUPEES
!       NSETOT : NOMBRE TOTAL DE SOUS-ELEMENT
!       NNNTOT : NOMBRE TOTAL DE NOUVEAUX NOEUDS
!       NCOTOT : LONGUEUR DE LA CONNECTIVITE DES NOUVEAUX NOEUDS
!       LOGRMA : LONGUEUR DES NOUVEAUX GROUP_MA
!       LISTGR : LISTE DES GROUPES CONTENANT CHAQUE MAILLE
!
    integer(kind=8) :: nbma, isepma, jcesd, jcesl, iad, ima
    integer(kind=8) :: nbman, nbmac, nbmax, ngr, igr, j1, n1, nbelt, iel, nsetot
    integer(kind=8) :: imac, imax, jmac, jmax, nnntot, ncotot, nse, n, nbgma
    integer(kind=8) ::  ndime, iret1, jlogma
    character(len=8) :: massmo, typma
    character(len=19) :: ces, ligrel
    character(len=24) :: sepmai, liel
    integer(kind=8), pointer :: tmdim(:) => null()
    integer(kind=8), pointer :: tma(:) => null()
    integer(kind=8), pointer :: typm(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
    parameter(massmo='NON')
!
!
    call jemarq()
!
    call dismoi('NB_MA_MAILLA', malini, 'MAILLAGE', repi=nbma)
    call jeveuo(malini//'.TYPMAIL', 'L', vi=tma)
!
!     TABLEAU D'ENTIERS DIMENSIONNÉ AU NOMBRE DE MAILLE DU MAILLAGE
!     INITIAL, INDIQUANT L'APPARTENANCE DES MAILLES A UN DES 2 GROUPES :
!         - 0 SI LA MAILLE N'EST PAS AFFECTEE D'UN MODELE -> MAILC
!         - 1 SI LA MAILLE N'EST PAS SOUS-DECOUPEE        -> MAILC
!         - 2 SI LA MAILLE EST SOUS-DECOUPEE              -> MAILX
    sepmai = '&XPOSEP.SEPMAI'
    call wkvect(sepmai, 'V V I', nbma, isepma)
!
!     DANS SEPMAI, MISE A -1 DE TOUTES LES MAILLES DU MODELE
!     LES MAILLES SANS MODELE RESTENT A 0
!     DE PLUS, LES MAILLES POI1 SONT CONSIDEREES COMME SANS MODELE
!     ET NE SERONT PAS POST TRAITEES
    ligrel = mo//'.MODELE'
    liel = ligrel//'.LIEL'
    call jelira(liel, 'NMAXOC', ngr)
    do igr = 1, ngr
        call jeveuo(jexnum(liel, igr), 'L', j1)
        call jelira(jexnum(liel, igr), 'LONMAX', n1)
        nbelt = n1-1
        do iel = 1, nbelt
            ima = zi(j1-1+iel)
            call jenuno(jexnum('&CATA.TM.NOMTM', tma(ima)), typma)
            if (typma .eq. 'POI1') goto 100
            zi(isepma-1+ima) = -1
        end do
100     continue
    end do
!
    ces = '&&XPOSEP.TOPOSE.LON'
    call celces(mo//'.TOPOSE.LON', 'V', ces)
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESV', 'L', vi=cesv)
    call jeveuo(ces//'.CESL', 'L', jcesl)
!
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=tmdim)
    call jeveuo(malini//'.TYPMAIL', 'L', vi=typm)
!
!     CREATION DE LA LISTE DES GROUPES CONTENANT CHAQUE MAILLE
    call ligrma(malini, listgr)
!
!     CREATION DU VECTEUR DIMENTIONNANT LA TAILLE DES NOUVEAUX GROUP_MA
    call jeexin(malini//'.GROUPEMA', iret1)
    nbgma = 0
    if (iret1 .gt. 0) call jelira(malini//'.GROUPEMA', 'NUTIOC', nbgma)
    if (nbgma .gt. 0) call wkvect(logrma, 'V V I', nbgma, jlogma)
!
    nbmac = 0
    nbmax = 0
    nsetot = 0
    nnntot = 0
    ncotot = 0
!
!     BOUCLE SUR LES MAILLES
    do ima = 1, nbma
        call jenuno(jexnum('&CATA.TM.NOMTM', tma(ima)), typma)
!
!       SI MASSMO = 'NON' : ON ZAPPE LES MAILLES SANS MODELE
        if (massmo .eq. 'NON' .and. zi(isepma-1+ima) .eq. 0) goto 200
!
!       RECUPERATION DE NSE
        call cesexi('C', jcesd, jcesl, ima, 1, &
                    1, 1, iad)
        if (iad .ne. 0) then
            nse = cesv(iad)
        else
            nse = 0
        end if
!
!       SI NSE DIFFERENT DE 0 ALORS MAILLE SOUS DECOUPEE
!       SINON, MAILLE CLASSIQUE
        if (nse .ne. 0) then
            zi(isepma-1+ima) = 2
            nbmax = nbmax+1
!         DIMENSION TOPOLOGIQUE DE LA MAILLE
            ndime = tmdim(typm(ima))
!
!         AUGMENTATION DE NSETOT AVEC LE NOMBRE DE NSE SUR LA MAILLE
            nsetot = nsetot+nse
!         AUGMENTATION DU NOMBRE DE NOUVEAUX NOEUDS (NNNTOT)
            nnntot = nnntot+cesv(iad+2)
!         AUGMENTATION DU NOMBRE DE NOEUDS DANS LA CONNECTIVITE TOT
            if (ismali(typma)) then
                ncotot = ncotot+nse*(ndime+1)
            else
                if (ndime .eq. 1) ncotot = ncotot+nse*3
                if (ndime .eq. 2) ncotot = ncotot+nse*6
                if (ndime .eq. 3) ncotot = ncotot+nse*10
            end if
!
!         AUGMENTATION DE LA TAILLE DES GROUP_MA
            call xpogma(nbgma, nse, listgr, ima, jlogma)
!
        else
            zi(isepma-1+ima) = 1
            nbmac = nbmac+1
!         N : NOMBRE DE NOEUDS DE LA MAILLE
            call jelira(jexnum(malini//'.CONNEX', ima), 'LONMAX', n)
            ncotot = ncotot+n
!
!         AUGMENTATION DE LA TAILLE DES GROUP_MA
            call xpogma(nbgma, 1, listgr, ima, jlogma)
!
        end if
!
200     continue
    end do
!
!     NOMBRE DE MAILLES NON TRAITEES
    nbman = nbma-nbmac-nbmax
!
    call utmess('I', 'XFEM_7', sk='NON TRAITEES', si=nbman)
    call utmess('I', 'XFEM_7', sk='CLASSIQUES', si=nbmac)
    call utmess('I', 'XFEM_7', sk='X-FEM', si=nbmax)
!
    imac = 0
    imax = 0
!
!     CREATION DES 2 GROUPES DE MAILLES DE SORTIE
    if (nbmac .ne. 0) call wkvect(mailc, 'V V I', nbmac, jmac)
    if (nbmax .ne. 0) call wkvect(mailx, 'V V I', nbmax, jmax)
    do ima = 1, nbma
        if (zi(isepma-1+ima) .eq. 1) then
            imac = imac+1
            zi(jmac-1+imac) = ima
        else if (zi(isepma-1+ima) .eq. 2) then
            imax = imax+1
            zi(jmax-1+imax) = ima
        end if
    end do
!
    call jedetr(sepmai)
!
    call detrsd('CHAM_ELEM_S', ces)
!
    call jedema()
end subroutine
