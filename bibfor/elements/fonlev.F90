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

subroutine fonlev(resu, noma, nbnoff)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: resu, noma
    integer(kind=8) :: nbnoff
! FONCTION REALISEE:
!
!     VERIFICATION DES LEVRES ET DE LEUR CONNEXITE
!
!     ENTREES:
!        RESU       : NOM DU CONCEPT RESULTAT DE L'OPERATEUR
!        NOMA       : NOM DU MAILLAGE
!        NBNOFF     : NOMBRE DE NOEUDS EN FOND DE FISSURE
!
!     -----------------------------------------------------------------
!
    integer(kind=8) :: jmai1, jadr, jnoe1, jmai2, jmaii, jjj, iatyma
    integer(kind=8) ::   iamase, ityp
    integer(kind=8) :: igr, ngr, i, j, k, ibid, k2, j2
    integer(kind=8) :: nbmai, indice
    integer(kind=8) :: nn, compta, nbmas1, nbmas2, nbmal
    integer(kind=8) :: iret, iret1, iret2, jjj2
    character(len=4) :: typma
    character(len=6) :: nompro
    character(len=8) :: maille, type, noeug, typmcl(2), motcle(2)
    character(len=9) :: typlev(2), motfac, valk(2)
    character(len=24) :: nomobj, grouma, conec, trav, trav2
    character(len=8), pointer :: inf(:) => null()
    character(len=8), pointer :: sup(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
!
    nompro = 'FONLEV'
!
!     ------------------------------------------------------------------
!     INITIALISATION DE VARIABLES
!     ------------------------------------------------------------------
    grouma = noma//'.GROUPEMA'
    conec = noma//'.CONNEX'
    call jeveuo(noma//'.TYPMAIL', 'L', iatyma)
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES DEUX LEVRES
!     CELLES-CI SONT TRAITEES DE LA MEME MANIERE
!     ------------------------------------------------------------------
    typlev(1) = 'LEVRE_SUP'
    typlev(2) = 'LEVRE_INF'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
    do indice = 1, 2
        motfac = typlev(indice)
!
!       EVALUATION DU NOMBRE DE MAILLES ET CONSTRUCTION DU VECTEUR DE
!       MAILLES DE TRAVAIL
!
        trav = '&&'//nompro//'.'//motfac
        trav2 = '&&'//nompro//'.'//motfac//'2'
        call reliem(' ', noma, 'NO_MAILLE', motfac, 1, &
                    2, motcle, typmcl, trav, nbmal)
        if (nbmal .eq. 0) goto 999
        call wkvect(trav2, 'V V K24', nbmal, jjj2)
        call jeveuo(trav, 'L', jjj)
!
!     ------------------------------------------------------------------
!     VERIFICATION DE L'EXISTENCE DES GROUPES DE MAILLES RENSEIGNEES
!     ET CALCUL DU NOMBRE TOTAL DE MAILLES
!     ------------------------------------------------------------------
!
        call getvtx(motfac, 'GROUP_MA', iocc=1, nbval=nbmal, vect=zk24(jjj2), &
                    nbret=ngr)
!
!
! ---   ALLOCATION D'UN PREMIER OBJET DE TRAVAIL
!       DANS LEQUEL SERONT STOCKES LES MAILLES AVANT DE S'ASSURER QU'IL
!       N'Y A PAS DUPPLICATION
!
        call wkvect('&&'//nompro//'.MAIL', 'V V K8', nbmal, jmai1)
        jmaii = jmai1
!       ----------------------------------------------------------------
!      VERIFICATION POUR LES MAILLES DE LA LEVRE COURANTE
!      SI ON A UN SEUL NOEUD ALORS ELLES SONT DE TYPE SEG
!      SI ON A PLUSIEURS NOEUDS ALORS ELLES SONT DE TYPE QUAD OU TRIA
!      ET CALCUL DU NOMBRE TOTAL DE MAILLES DE LA LEVRE COURANTE
!      ----------------------------------------------------------------
!       SI GROUP_MA
        do igr = 1, ngr
!
            call jelira(jexnom(grouma, zk24(jjj2-1+igr)), 'LONMAX', nbmai)
            call jeveuo(jexnom(grouma, zk24(jjj2-1+igr)), 'L', jadr)
!
!
            do i = 1, nbmai
                maille = int_to_char8(zi(jadr-1+i))
                ibid = char8_to_int(maille)
                ityp = iatyma-1+ibid
                call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
                typma = type(1:4)
                if (((typma .ne. 'QUAD') .and. (typma .ne. 'TRIA')) .and. (nbnoff .gt. 1)) then
                    valk(1) = type
                    valk(2) = motfac
                    call utmess('F', 'RUPTURE0_65', nk=2, valk=valk)
                elseif ((typma(1:3) .ne. 'SEG') .and. (nbnoff .eq. 1)) &
                    then
                    valk(1) = type
                    valk(2) = motfac
                    call utmess('F', 'RUPTURE0_75', nk=2, valk=valk)
                else
                    zk8(jmai1) = maille
                    jmai1 = jmai1+1
                end if
!
            end do
!
        end do
!
! --- VERIFICATION QU'IL N Y A PAS DUPLICATION DES ENTITES ET STOCKAGE
!     ON MET 'O' SI L'ENTITE EST DUPPLIQUEE
!
!       ALLOCATION DU VECTEUR .LEVRESUP.MAIL ET .LEVREINF.MAIL
        nomobj = resu//'.LEVRE'//motfac(7:9)//'.MAIL'
        call wkvect(nomobj, 'G V K8', nbmal, jmai2)
        k2 = 1
        do i = 1, nbmal-1
            if (zk8(jmaii-1+i) .ne. '0') then
                zk8(jmai2-1+k2) = zk8(jmaii-1+i)
                k2 = k2+1
                do j = i+1, nbmal
                    if (zk8(jmaii-1+i) .eq. zk8(jmaii-1+j)) then
                        zk8(jmaii-1+j) = '0'
                        j2 = i
                    end if
                end do
            end if
        end do
        if (zk8(jmaii-1+nbmal) .ne. '0') then
            zk8(jmai2-1+k2) = zk8(jmaii-1+nbmal)
            k2 = k2+1
        end if
        k2 = k2-1
!
        if (k2 .ne. nbmal) then
            valk(1) = motfac
            valk(2) = zk8(jmaii-1+j2)
            call utmess('E', 'RUPTURE0_70', nk=2, valk=valk)
        end if
!
! --- VERIFICATION COHERENCE LEVRE SUP / FOND
!
        call jeexin(resu//'.FOND.NOEU', iret)
        if (iret .ne. 0) then
            call jelira(resu//'.FOND.NOEU', 'LONUTI', nbnoff)
            call jeveuo(resu//'.FOND.NOEU', 'L', jnoe1)
        else
            ASSERT(.FALSE.)
        end if
        if (nbnoff .gt. 1) then
            do i = 1, nbnoff
                compta = 0
                do j = 1, nbmal
                    maille = int_to_char8(zi(jadr-1+j))
                    ibid = char8_to_int(maille)
                    ityp = iatyma-1+ibid
                    call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
                    call dismoi('NBNO_TYPMAIL', type, 'TYPE_MAILLE', repi=nn)
                    if ((type(1:5) .ne. 'QUAD8') .and. (type(1:5) .ne. 'TRIA3') .and. &
                        (type(1:5) .ne. 'QUAD4') .and. (type(1:5) .ne. 'TRIA6') .and. &
                        (type(1:5) .ne. 'QUAD9') .and. (type(1:5) .ne. 'TRIA7')) then
                        valk(1) = type(1:5)
                        valk(2) = motfac
                        call utmess('F', 'RUPTURE0_65', nk=2, valk=valk)
                    end if
                    call jeveuo(jexnum(conec, ibid), 'L', iamase)
                    noeug = int_to_char8(zi(iamase))
                    do k = 1, nn
                        noeug = int_to_char8(zi(iamase-1+k))
                        if (noeug .eq. zk8(jnoe1-1+i)) then
                            compta = compta+1
                            goto 610
                        end if
                    end do
                end do
                if (compta .eq. 0) then
                    valk(1) = zk8(jnoe1-1+i)
                    valk(2) = motfac
                    call utmess('F', 'RUPTURE0_72', nk=2, valk=valk)
                end if
610             continue
            end do
        end if
!
!
! --- DESTRUCTION DES OBJETS DE TRAVAIL
!
        call jedetr(trav)
        call jedetr(trav2)
!
        call jedetr('&&'//nompro//'.MAIL')
!
    end do
! ----------------------------------------------------------
!    COMPARAISON LEVRE SUP / LEVRE INF AFIN DE S'ASSURER
!    QU'ELLES N'ONT PAS DE MAILLES EN COMMUN
! ----------------------------------------------------------
    call jeexin(resu//'.LEVRESUP.MAIL', iret1)
    call jeexin(resu//'.LEVREINF.MAIL', iret2)
    if ((iret1 .ne. 0) .and. (iret2 .ne. 0)) then
        call jeveuo(resu//'.LEVRESUP.MAIL', 'L', vk8=sup)
        call jeveuo(resu//'.LEVREINF.MAIL', 'L', vk8=inf)
        call jelira(resu//'.LEVRESUP.MAIL', 'LONMAX', nbmas1)
        call jelira(resu//'.LEVRESUP.MAIL', 'LONMAX', nbmas2)
        do i = 1, nbmas1
            do j = 1, nbmas2
                if (sup(i) .eq. inf(j)) then
                    call utmess('F', 'RUPTURE0_73', sk=sup(i))
                end if
            end do
        end do
    end if
999 continue
!
    call jedema()
!
end subroutine
