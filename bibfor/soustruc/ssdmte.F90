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
subroutine ssdmte(mag)
    implicit none
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: mag
! ----------------------------------------------------------------------
!     BUT:
!        - TERMINER LE TRAITEMENT
!          DES COMMANDES DEFI_MAILLAGE ET CONC_MAILLAGE.
!        - CREER LES OBJETS :
!            BASE GLOBALE : .COORDO , .NOMNOE
!        - MODIFIER LES OBJETS :
!            BASE GLOBALE : .SUPMAIL, .GROUPENO ET .CONNEX
!            POUR TENIR COMPTE DES NOEUDS CONFONDUS.
!
!     IN:
!        MAG : NOM DU MAILLAGE RESULTAT.
!
    character(len=8) :: nomacr, nomnoe
    aster_logical :: recom
    character(len=19) :: coordo
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
!     NBNOPH : NOMBRE DE NOEUDS PHYSIQUES DE MAG "AVANT"
!              (AVANT DE CONFONDRE CERTAINS NOEUDS)
!     NBNOP2 : NOMBRE DE NOEUDS PHYSIQUES DE MAG "APRES"
!              (APRES AVOIR CONFONDU CERTAINS NOEUDS)
!     NBNOCO : NOMBRE DE NOEUDS CONFONDUS.
!              (NBNOCO=NBNOP2-NBNOPH)
!     NBNOLA : NOMBRE DE NOEUDS "LAGRANGE" DE MAG.
!     NBNOT2 : NOMBRE DE NOEUDS TOTAL DE MAG "APRES"
!              (APRES AVOIR CONFONDU CERTAINS NOEUDS)
!              (NBNOT2= NBNOLA+NBNOP2)
!
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i2coex, iacoex, iadesc
    integer(kind=8) :: iagno, ianeno
    integer(kind=8) :: iasupm, iatypl, iavale, ibid, ico, igeomr, igno
    integer(kind=8) :: ilcoex, ima, ino, iret, isma, jno, k
    integer(kind=8) :: kno, nbgno, nbma, nbno, nbnoco, nbnoe, nbnoet
    integer(kind=8) :: nbnogn, nbnol, nbnola, nbnop2, nbnoph, nbnot2, nbsma
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: noeud_conf(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
    character(len=8), pointer :: nomnoe_2(:) => null()
    real(kind=8), pointer :: coordo_2(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: dime_2(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(mag//'.DIME', 'E', vi=dime)
    nbnoph = dime(1)
    nbnola = dime(2)
    nbma = dime(3)
    nbsma = dime(4)
    call jeveuo(mag//'.COORDO_2', 'L', vr=coordo_2)
    call jeveuo(mag//'.NOEUD_CONF', 'E', vi=noeud_conf)
    call jeveuo(mag//'.NOMNOE_2', 'L', vk8=nomnoe_2)
!
    if (nbsma .gt. 0) call jeveuo(mag//'.DIME_2', 'L', vi=dime_2)
    if (nbsma .gt. 0) call jeveuo(mag//'.NOMACR', 'L', vk8=vnomacr)
!
!
!     -- ON COMPTE LES NOEUDS PHYSIQUES REELLEMENT CONSERVES :
!     -------------------------------------------------------
    ico = 0
    do ino = 1, nbnoph
        if (noeud_conf(ino) .eq. ino) ico = ico+1
    end do
    nbnop2 = ico
    nbnot2 = nbnop2+nbnola
    nbnoco = nbnoph-nbnop2
!
    call jecreo(mag//'.NOMNOE', 'G N K8')
    call jeecra(mag//'.NOMNOE', 'NOMMAX', nbnot2)
!
!
!     -- CREATION DE .TYPL :
!     ----------------------
    if (nbnola .gt. 0) then
        call wkvect(mag//'.TYPL', 'G V I', nbnola, iatypl)
        do isma = 1, nbsma
            nomacr = vnomacr(isma)
            call jeveuo(nomacr//'.CONX', 'L', vi=conx)
            call jeveuo(jexnum(mag//'.SUPMAIL', isma), 'L', iasupm)
            nbnoe = dime_2(4*(isma-1)+1)
            nbnol = dime_2(4*(isma-1)+2)
            nbnoet = nbnoe+nbnol
            do i = 1, nbnoet
                ino = zi(iasupm-1+i)
                if (ino .gt. nbnoph) then
                    zi(iatypl-1+ino-nbnoph) = conx(3*(i-1)+3)
                end if
            end do
        end do
    end if
!
!
!     -- CREATION DU CHAMP .COORDO :
!     ------------------------------
    coordo = mag//'.COORDO'
!
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), igeomr)
    call wkvect(coordo//'.DESC', 'G V I', 3, iadesc)
    call jeecra(coordo//'.DESC', 'DOCU', ibid, 'CHGO')
    zi(iadesc-1+1) = igeomr
!     -- TOUJOURS 3 COMPOSANTES X, Y ET Z
    zi(iadesc-1+2) = -3
!     -- 14 = 2**1 + 2**2 + 2**3
    zi(iadesc-1+3) = 14
!
    call wkvect(coordo//'.VALE', 'G V R', 3*nbnop2, iavale)
!     -- NOM DES NOEUDS PHYSIQUES (ET LEUR COORDONNEES) :
    ico = 0
    do 3, ino = 1, nbnoph
        jno = noeud_conf(ino)
        if (ino .ne. jno) goto 3
        ico = ico+1
        if (nomnoe_2(ino) .ne. ' ') then
            nomnoe = nomnoe_2(ino)
        else
            nomnoe = 'N?'
            call codent(ico, 'G', nomnoe(2:8))
        end if
        call jecroc(jexnom(mag//'.NOMNOE', nomnoe))
        do k = 1, 3
            zr(iavale-1+3*(ico-1)+k) = coordo_2(3*(ino-1)+k)
        end do
3   end do
!     -- NOM DES NOEUDS DE LAGRANGE :
    nomnoe = '&?'
    do 4, ino = 1, nbnola
        call codent(ino, 'G', nomnoe(2:8))
        call jecroc(jexnom(mag//'.NOMNOE', nomnoe))
4   end do
!
!
!     -- ON OTE LA "RECURSIVITE" DE .NOEUD_CONF:
!     ------------------------------------------
5   continue
    recom = .false.
    do 6, ino = 1, nbnoph
        jno = noeud_conf(ino)
        if (jno .ne. ino) then
            ASSERT(jno .le. ino)
            kno = noeud_conf(jno)
            if (kno .ne. jno) then
                noeud_conf(ino) = kno
                recom = .true.
            end if
        end if
6   end do
    if (recom) goto 5
!
!
!     -- ON COMPACTE LES NUMEROS DES NOEUDS CONSERVES:
!     ------------------------------------------------
    call wkvect(mag//'.NENO', 'V V I', nbnoph, ianeno)
    ico = 0
    do 7, ino = 1, nbnoph
        jno = noeud_conf(ino)
        if (jno .eq. ino) then
            ico = ico+1
            zi(ianeno-1+ino) = ico
        end if
7   end do
!
!
!     -- MODIFICATION DES OBJETS POUR TENIR COMPTE DE .NOEUD_CONF:
!     -------------------------------------------------------------
!
!     -- MODIFICATION DE .CONNEX:
!     ---------------------------
    if (nbma .gt. 0) then
        call jeveuo(mag//'.CONNEX', 'E', iacoex)
        call jeveuo(jexatr(mag//'.CONNEX', 'LONCUM'), 'L', ilcoex)
    end if
    do ima = 1, nbma
        nbno = zi(ilcoex-1+ima+1)-zi(ilcoex-1+ima)
        i2coex = iacoex-1+zi(ilcoex-1+ima)
        do i = 1, nbno
            ino = zi(i2coex-1+i)
            if (ino .le. nbnoph) then
                jno = zi(ianeno-1+noeud_conf(ino))
                zi(i2coex-1+i) = jno
            else
                zi(i2coex-1+i) = ino-nbnoco
            end if
        end do
    end do
!
!     -- MODIFICATION DE .SUPMAIL:
!     ----------------------------
    do isma = 1, nbsma
        call jeveuo(jexnum(mag//'.SUPMAIL', isma), 'E', iasupm)
        nbnoe = dime_2(4*(isma-1)+1)
        nbnol = dime_2(4*(isma-1)+2)
        nbnoet = nbnoe+nbnol
        do i = 1, nbnoet
            ino = zi(iasupm-1+i)
            if (ino .le. nbnoph) then
                jno = zi(ianeno-1+noeud_conf(ino))
                zi(iasupm-1+i) = jno
            else
                zi(iasupm-1+i) = ino-nbnoco
            end if
        end do
    end do
!
!     -- MODIFICATION DE .GROUPENO:
!     ----------------------------
    call jeexin(mag//'.GROUPENO', iret)
    if (iret .gt. 0) then
        call jelira(mag//'.GROUPENO', 'NUTIOC', nbgno)
        do 9, igno = 1, nbgno
            call jeveuo(jexnum(mag//'.GROUPENO', igno), 'E', iagno)
            call jelira(jexnum(mag//'.GROUPENO', igno), 'LONUTI', nbnogn)
            do i = 1, nbnogn
                ino = zi(iagno-1+i)
                jno = zi(ianeno-1+noeud_conf(ino))
                zi(iagno-1+i) = jno
            end do
9           continue
            end if
!
!
!     -- REMISE A JOUR DEFINITIVE DU NOMBRE DE NOEUDS PHYSIQUES:
!     ----------------------------------------------------------
            dime(1) = nbnop2
!
!
            call jedema()
            end subroutine
