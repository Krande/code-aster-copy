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

subroutine pjmasp(moa2, masp, corres, noca)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
! ----------------------------------------------------------------------
! COMMANDE PROJ_CHAMP / METHODE='SOUS_POINT'
!
! BUT :  CREER UN MAILLAGE (MASP) DONT LES NOEUDS SONT POSITIONNES SUR
!        LES SOUS-POINTS DE GAUSS D'UN MODELE (MOA2) POUR CHAQUE
!        FAMILLE DE POINTS DE LA LISTE MATER.
! ----------------------------------------------------------------------
! IN MOA2 : MODELE "2"
! IN/JXOUT MASP : MAILLAGE 2 PRIME (OBTENU A PARTIR DES PG DU MODELE 2)
! IN/JXVAR : ON CREE L'OBJET CORRES.PJEF_SP
! ----------------------------------------------------------------------
    character(len=16) :: corres
    character(len=8) :: masp, moa2, noca
! ----------------------------------------------------------------------
    integer(kind=8) :: ntgeo, ipo, ipg, nuno2
    integer(kind=8) ::  nbnosp, nno2, ino2p
    integer(kind=8) ::  j1, ipoi1
    integer(kind=8) :: nbma, nbpt, nbsp, nbcmp
    integer(kind=8) :: ima, ipt, isp, icmp, iad, iadime
    integer(kind=8) :: jtypma, jpo2
    integer(kind=8) :: jcesd, jcesl, iatypm
    integer(kind=8) :: nchi, nbpgmx, nbspmx
    character(len=8) ::  mail2, lpain(6)
    character(len=19) :: chamg, ces, chgeom, ligrel
    character(len=24) :: coodsc
    character(len=24) :: lchin(6)
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
! ----------------------------------------------------------------------
    call jemarq()
!
!
!     -- RECUPERATION DU NOM DU MAILLAGE 2
    call dismoi('NOM_MAILLA', moa2, 'MODELE', repk=mail2)
    call jeveuo(mail2//'.TYPMAIL', 'L', jtypma)
!
!     -- RECUPERATION DU CHAMP DE COORDONNEES DU MAILLAGE 2
    chgeom = mail2//'.COORDO'
!
    call dismoi('NOM_LIGREL', moa2, 'MODELE', repk=ligrel)
!     1.  CALCUL DU CHAMP DE COORDONNEES DES ELGA (CHAMG):
!     -------------------------------------------------------
!
    nchi = 6
    lchin(1) = chgeom(1:19)
    lpain(1) = 'PGEOMER'
    lchin(2) = noca//'.CARORIEN'
    lpain(2) = 'PCAORIE'
    lchin(3) = noca//'.CAFIBR'
    lpain(3) = 'PFIBRES'
    lchin(4) = noca//'.CANBSP'
    lpain(4) = 'PNBSP_I'
    lchin(5) = noca//'.CARCOQUE'
    lpain(5) = 'PCACOQU'
    lchin(6) = noca//'.CARGEOPO'
    lpain(6) = 'PCAGEPO'
    chamg = '&&PJMASP.PGCOOR'
    call cesvar(noca, ' ', ligrel, chamg)
    call calcul('S', 'COOR_ELGA_MATER', ligrel, nchi, lchin, &
                lpain, 1, chamg, 'PCOOPGM', 'V', &
                'OUI')
!
!     -- TRANSFORMATION DE CE CHAMP EN CHAM_ELEM_S
    ces = '&&PJMASP.PGCORS'
    call celces(chamg, 'V', ces)
!
!
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESL', 'L', jcesl)
    call jeveuo(ces//'.CESV', 'E', vr=cesv)
    nbma = zi(jcesd-1+1)
!
!
!     2. CALCUL DE NBNOSP : NOMBRE DE NOEUDS (ET DE MAILLES) DE MASP
!        CALCUL DE '.PJEF_SP'
!     ----------------------------------------------------------------
    nbnosp = 0
!
!
    nbpgmx = zi(jcesd-1+3)
    nbspmx = zi(jcesd-1+4)
!
!     NBMA*NBPGMX*NBSPMX*3 = NB MAX DE MAILLES * NB DE PG MAX  *
!                              NB DE SP MAX * 3
!     ON CREE UN TABLEAU, POUR CHAQUE JPO2, ON STOCKE TROIS VALEURS :
!      * LA PREMIERE VALEUR EST LE NUMERO DE LA MAILLE
!      * LA DEUXIEME VALEUR EST LE NUMERO DU PG DANS CETTE MAILLE
!      EX : LE PG 3 DE LA FAMILLE 2 DE LA LISTE MATER AURA LE NUMERO
!      DE PG IPG = NB DE PG DE LA FAMILLE 1 + 3
!      * LA TROISIEME VALEUR EST LE NUMERO DU SOUS-POINT
!
    call wkvect(corres//'.PJEF_SP', 'V V I', nbma*nbpgmx*nbspmx*3, jpo2)
!
    ipo = 1
    do ima = 1, nbma
        nbpt = zi(jcesd-1+5+4*(ima-1)+1)
        nbsp = zi(jcesd-1+5+4*(ima-1)+2)
!          IF (NBPT.LT.1) GOTO 100
        if (nbsp .lt. 1) goto 100
        do ipg = 1, nbpt
            do isp = 1, nbsp
                zi(jpo2-1+ipo) = ima
                zi(jpo2-1+ipo+1) = ipg
                zi(jpo2-1+ipo+2) = isp
                ipo = ipo+3
            end do
        end do
        nbnosp = nbnosp+nbpt*nbsp
100     continue
    end do
!
!     3. CREATION DU .DIME DU NOUVEAU MAILLAGE
!        IL Y A AUTANT DE MAILLES QUE DE NOEUDS
!        TOUTES LES MAILLES SONT DES POI1
!     --------------------------------------------------
    call wkvect(masp//'.DIME', 'V V I', 6, iadime)
    zi(iadime-1+1) = nbnosp
    zi(iadime-1+3) = nbnosp
    zi(iadime-1+6) = 3
!
!
!
!     5. CREATION DU .CONNEX ET DU .TYPMAIL DU NOUVEAU MAILLAGE
!     ----------------------------------------------------------
    call jecrec(masp//'.CONNEX', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbnosp)
    call jeecra(masp//'.CONNEX', 'LONT', nbnosp, ' ')
    call jeveuo(masp//'.CONNEX', 'E', vi=connex)
!
    call wkvect(masp//'.TYPMAIL', 'V V I', nbnosp, iatypm)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ipoi1)
!
    nuno2 = 0
    do ima = 1, nbnosp
        zi(iatypm-1+ima) = ipoi1
        nno2 = 1
        call jecroc(jexnum(masp//'.CONNEX', ima))
        call jeecra(jexnum(masp//'.CONNEX', ima), 'LONMAX', nno2)
        nuno2 = nuno2+1
        connex(nuno2) = nuno2
    end do
!
!
!
!     -- CREATION DE COORDO.VALE DU NOUVEAU MAILLAGE
!     --------------------------------------------------
    call wkvect(masp//'.COORDO    .VALE', 'V V R', 3*nbnosp, j1)
!
    ino2p = 0
    do ima = 1, nbma
        nbpt = zi(jcesd-1+5+4*(ima-1)+1)
        nbsp = zi(jcesd-1+5+4*(ima-1)+2)
        nbcmp = zi(jcesd-1+5+4*(ima-1)+3)
        if (nbpt .eq. 0) goto 160
        ASSERT(nbcmp .eq. 3)
        do ipt = 1, nbpt
            do isp = 1, nbsp
                ino2p = ino2p+1
                do icmp = 1, 3
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        zr(j1-1+3*(ino2p-1)+icmp) = cesv(iad)
                    end if
                end do
            end do
        end do
!
160     continue
    end do
    ASSERT(ino2p .eq. nbnosp)
!
!
!     -- CREATION DU .DESC DU NOUVEAU MAILLAGE
!     --------------------------------------------------
    coodsc = masp//'.COORDO    .DESC'
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
