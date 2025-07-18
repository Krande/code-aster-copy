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

subroutine pjma2p(ndim, moa2, ma2p, corres)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
! ----------------------------------------------------------------------
! COMMANDE PROJ_CHAMP / METHODE='ECLA_PG'
!
! BUT :  CREER UN MAILLAGE (MA2P) DONT LES NOEUDS SONT POSITIONNES SUR
!        LES POINTS DE GAUSS D'UN MODELE (MOA2).
! REMARQUE : ON UTILISE L'OPTION COOR_ELGA CE QUI CORRESPOND EN GENERAL
!            A LA FAMILLE DE POINTS DE GAUS "RIGI"
! ----------------------------------------------------------------------
! IN NDIM : 2/3 : DIMENSION DES MAILLES A PROJETER
! IN MOA2 : MODELE "2"
! IN/JXOUT MA2P : MAILLAGE 2 PRIME (OBTENU A PARTIR DES PG DU MODELE 2)
! IN/JXVAR : ON COMPLETE LA SD_CORRESP_2_MAILLA AVEC L'OBJET .PJEL
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/calc_coor_elga.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utflmd.h"
#include "asterfort/utmamo.h"
#include "asterfort/wkvect.h"
!
    character(len=16) :: corres
    character(len=8) :: ma2p, moa2
    integer(kind=8) :: ndim
! ----------------------------------------------------------------------
    integer(kind=8) :: ntgeo, ipo, ipg, nuno2
    integer(kind=8) ::  nbno2p, nno2, ino2p
    integer(kind=8) ::  j1, ipoi1, ipy5, ipy13
    integer(kind=8) :: nbma, nbpt, nbcmp, nbmamo
    integer(kind=8) :: ima, ipt, icmp, iad, iadime
    integer(kind=8) ::  jdimt, jpo2, nbtrou, jlitr
    integer(kind=8) :: jcesd, jcesl, iatypm
    character(len=8) ::  mail2
    character(len=19) :: chamg, ces, chgeom, ligrel
    character(len=24) :: coodsc, limato, litrou
    real(kind=8) :: xmoy(3), rayo
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
! ----------------------------------------------------------------------
    call jemarq()
!
!     -- RECUPERATION DU NOM DU MAILLAGE 2
    call dismoi('NOM_MAILLA', moa2, 'MODELE', repk=mail2)
    call jeveuo(mail2//'.TYPMAIL', 'L', vi=typmail)
!
!     -- RECUPERATION DU CHAMP DE COORDONNEES DU MAILLAGE 2
    chgeom = mail2//'.COORDO'
!
!
!     -- ON REDUIT LE LIGREL DE MOA2 SUR MAILLES DE DIMENSION NDIM :
    ligrel = '&&PJMA2P.LIGREL'
    limato = '&&PJMA2P.LIMATOT'
    litrou = '&&PJMA2P.LITROU'
!     -- ON NE CONSERVE QUE LES MAILLES AFFECTEES :
    call utmamo(moa2, nbmamo, limato)
!     -- ON NE CONSERVE QUE LES MAILLES DE DIMENSION NDIM :
    call utflmd(mail2, limato, nbmamo, ndim, ' ', &
                nbtrou, litrou)
    ASSERT(nbtrou .gt. 0)
    call jeveuo(litrou, 'L', jlitr)
    call exlim1(zi(jlitr), nbtrou, moa2, 'V', ligrel)
    call jedetr(limato)
    call jedetr(litrou)
!
!
!     1.  CALCUL DU CHAMP DE COORDONNEES DES ELGA (CHAMG):
!     -------------------------------------------------------
    chamg = '&&PJMA2P.PGCOOR'
!
    call calc_coor_elga(moa2, ligrel, chgeom, chamg)
!
!     -- TRANSFORMATION DE CE CHAMP EN CHAM_ELEM_S
    ces = '&&PJMA2P.PGCORS'
    call celces(chamg, 'V', ces)
!
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESL', 'L', jcesl)
    call jeveuo(ces//'.CESV', 'E', vr=cesv)
    nbma = zi(jcesd-1+1)
!
!     2.1 MODIFICATION DES COORDONNEES DE CERTAINS PG (PYRAM/FPG27)
!         CAR CES POINTS DE GAUSS SONT EN "DEHORS" DES PYRAMIDES
!     ----------------------------------------------------------------
    call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM5'), ipy5)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM13'), ipy13)
    do ima = 1, nbma
        if (typmail(ima) .eq. ipy5 .or. typmail(ima) .eq. ipy13) then
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            if (nbpt .eq. 27) then
                do icmp = 1, 3
                    xmoy(icmp) = 0.d0
                end do
!           -- XMOY : CENTRE DE LA PYRAMIDE :
                do ipt = 1, 15
                    do icmp = 1, 3
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    1, icmp, iad)
                        ASSERT(iad .gt. 0)
                        xmoy(icmp) = xmoy(icmp)+cesv(iad)
                    end do
                end do
                do icmp = 1, 3
                    xmoy(icmp) = xmoy(icmp)/15
                end do
!
!           -- ON "RAMENE" LES 12 DERNIERS PG VERS LE CENTRE (10%) :
                do ipt = 16, 27
                    do icmp = 1, 3
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    1, icmp, iad)
                        ASSERT(iad .gt. 0)
                        rayo = cesv(iad)-xmoy(icmp)
                        cesv(iad) = cesv(iad)-0.6d0*rayo
                    end do
                end do
            end if
        end if
    end do
!
!
!     2. CALCUL DE NBNO2P : NOMBRE DE NOEUDS (ET DE MAILLES) DE MA2P
!        CALCUL DE '.PJEF_EL'
!     ----------------------------------------------------------------
    nbno2p = 0
    do ima = 1, nbma
        call jeveuo(jexnum('&CATA.TM.TMDIM', typmail(ima)), 'L', jdimt)
        if (zi(jdimt) .eq. ndim) then
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            if (nbpt .eq. 0) cycle
            nbno2p = nbno2p+nbpt
        end if
    end do
!
!     ON CREE UN TABLEAU, POUR CHAQUE JPO2, ON STOCKE DEUX VALEURS :
!      * LA PREMIERE VALEUR EST LE NUMERO DE LA MAILLE
!      * LA DEUXIEME VALEUR EST LE NUMERO DU PG DANS CETTE MAILLE
    call wkvect(corres//'.PJEF_EL', 'V V I', nbno2p*2, jpo2)
!
    ipo = 1
    do ima = 1, nbma
        call jeveuo(jexnum('&CATA.TM.TMDIM', typmail(ima)), 'L', jdimt)
        if (zi(jdimt) .eq. ndim) then
            nbpt = zi(jcesd-1+5+4*(ima-1)+1)
            if (nbpt .eq. 0) cycle
            do ipg = 1, nbpt
                zi(jpo2-1+ipo) = ima
                zi(jpo2-1+ipo+1) = ipg
                ipo = ipo+2
            end do
        end if
    end do
!
!     3. CREATION DU .DIME DU NOUVEAU MAILLAGE
!        IL Y A AUTANT DE MAILLES QUE DE NOEUDS
!        TOUTES LES MAILLES SONT DES POI1
!     --------------------------------------------------
    call wkvect(ma2p//'.DIME', 'V V I', 6, iadime)
    zi(iadime-1+1) = nbno2p
    zi(iadime-1+3) = nbno2p
    zi(iadime-1+6) = 3
!
!
!
!     5. CREATION DU .CONNEX ET DU .TYPMAIL DU NOUVEAU MAILLAGE
!     ----------------------------------------------------------
    call jecrec(ma2p//'.CONNEX', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbno2p)
    call jeecra(ma2p//'.CONNEX', 'LONT', nbno2p, ' ')
    call jeveuo(ma2p//'.CONNEX', 'E', vi=connex)
!
    call wkvect(ma2p//'.TYPMAIL', 'V V I', nbno2p, iatypm)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ipoi1)
!
    nuno2 = 0
    do ima = 1, nbno2p
        zi(iatypm-1+ima) = ipoi1
        nno2 = 1
        call jecroc(jexnum(ma2p//'.CONNEX', ima))
        call jeecra(jexnum(ma2p//'.CONNEX', ima), 'LONMAX', nno2)
        nuno2 = nuno2+1
        connex(nuno2) = nuno2
    end do
!
!
!
!
!
!     -- CREATION DE COORDO.VALE DU NOUVEAU MAILLAGE
!     --------------------------------------------------
    call wkvect(ma2p//'.COORDO    .VALE', 'V V R', 3*nbno2p, j1)
!
    ino2p = 0
    do ima = 1, nbma
        nbpt = zi(jcesd-1+5+4*(ima-1)+1)
        nbcmp = zi(jcesd-1+5+4*(ima-1)+3)
        if (nbpt .eq. 0) goto 160
        call jeveuo(jexnum('&CATA.TM.TMDIM', typmail(ima)), 'L', jdimt)
!
        if (zi(jdimt) .eq. ndim) then
            ASSERT(nbcmp .ge. 3)
            do ipt = 1, nbpt
                ino2p = ino2p+1
                do icmp = 1, 3
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                1, icmp, iad)
                    if (iad .gt. 0) then
                        zr(j1-1+3*(ino2p-1)+icmp) = cesv(iad)
                    end if
                end do
            end do
        end if
160     continue
    end do
    ASSERT(ino2p .eq. nbno2p)
!
!
!     -- CREATION DU .DESC DU NOUVEAU MAILLAGE
!     --------------------------------------------------
    coodsc = ma2p//'.COORDO    .DESC'
!
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)
    call jecreo(coodsc, 'V V I')
    call jeecra(coodsc, 'LONMAX', 3)
    call jeecra(coodsc, 'DOCU', cval='CHGO')
    call jeveuo(coodsc, 'E', iad)
    zi(iad) = ntgeo
    zi(iad+1) = -3
    zi(iad+2) = 14
!
    call detrsd('CHAM_ELEM', chamg)
    call detrsd('CHAM_ELEM_S', ces)
!
    call jedema()
end subroutine
