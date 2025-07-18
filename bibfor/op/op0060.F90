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
subroutine op0060()
    implicit none
!
!      OPERATEUR :     DEFI_FOND_FISS
!
!-----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/fonbas2.h"
#include "asterfort/fonfis2.h"
#include "asterfort/fonimp.h"
#include "asterfort/foninf2.h"
#include "asterfort/fonlev.h"
#include "asterfort/fonmai2.h"
#include "asterfort/fonnoe2.h"
#include "asterfort/fonnof2.h"
#include "asterfort/fonnor2.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: iadr1, ifm, niv
    integer(kind=8) :: nbocc, nbnoff
    integer(kind=8) :: ibas, ibid, iocc, idon, idonn, ifonoe, ndonn, nvenor
    integer(kind=8) :: iret1, iret, irets, jnorm
    integer(kind=8) :: n1, n2
    character(len=6) :: nompro
    character(len=8) :: resu, noma, typfon, confin, typmp, typm
    character(len=9) :: entit(4)
    character(len=13) :: motcl(4)
    character(len=16) :: typres, oper
    character(len=19) :: basnof, basloc, cnxinv, lnno, ltno
    character(len=24) :: valk(2), entnom, abscur, fonoeu, absfon, coorfond
! DEB-------------------------------------------------------------------
!
    call jemarq()
    nompro = 'OP0060'
!
    call infniv(ifm, niv)
!
! ---  RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getres(resu, typres, oper)
!
! ---  RECUPERATIONS RELATIVES AU MAILLAGE
!      -----------------------------------
!
    call getvid(' ', 'MAILLAGE', scal=noma, nbret=nbocc)
    ASSERT(.not. isParallelMesh(noma))
!
! ---  RECUPERATION DE LA CONNECTIVITE INVERSE
!
    cnxinv = '&&'//nompro//'.CNXINV'
    call cncinv(noma, [ibid], 0, 'V', cnxinv)
!
!
!     ---------------------------------------------------------------
!     RECUPERATION DU TYPE DE FOND : OUVERT OU FERME
!     ---------------------------------------------------------------
!
    call getfac('FOND_FISS', nbocc)
    do iocc = 1, nbocc
!
        call getvtx('FOND_FISS', 'TYPE_FOND', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvtx('FOND_FISS', 'TYPE_FOND', iocc=iocc, scal=typfon, nbret=n1)
        else
            typfon = 'OUVERT'
        end if
!
!       DANS LE CAS D'UN FOND FERME, ON INTERDIT LES MCS NOEUD ET GROUP_NO
        if (typfon .eq. 'FERME') then
            call getvtx('FOND_FISS', 'GROUP_NO', iocc=iocc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call utmess('F', 'RUPTURE0_98')
            end if
        end if
!
!     ---------------------------------------------------------------
!     VERIFICATION DE L'EXISTANCE DES ENTITES DU MAILLAGE RENSEIGNEES
!     ET CONSTRUCTION DE VECTEURS DE TRAVAIL POUR CHACUNE D'ELLES
!     ---------------------------------------------------------------
!
        entit(1) = '.GROUPENO'
        entit(2) = '.GROUPEMA'
        entit(3) = '.GROUPENO'
        motcl(1) = 'GROUP_NO'
        motcl(2) = 'GROUP_MA'
        motcl(3) = 'GROUP_NO_ORIG'
        if (typfon .eq. 'OUVERT') then
            entit(4) = '.GROUPENO'
            motcl(4) = 'GROUP_NO_EXTR'
            ndonn = 4
        else if (typfon .eq. 'FERME') then
            entit(4) = '.GROUPEMA'
            motcl(4) = 'GROUP_MA_ORIG'
            ndonn = 4
        else
            ASSERT(.FALSE.)
        end if
        do idonn = 1, ndonn
            call getvtx('FOND_FISS', motcl(idonn), iocc=iocc, nbval=0, nbret=n1)
            n1 = -n1
            if (n1 .gt. 0) then
                call wkvect('&&'//nompro//'.'//motcl(idonn), 'V V K24', n1, iadr1)
                call getvtx('FOND_FISS', motcl(idonn), iocc=iocc, nbval=n1, vect=zk24(iadr1), &
                            nbret=n2)
                do idon = 1, n1
                    entnom = zk24(iadr1-1+idon)
                    call jenonu(jexnom(noma//entit(idonn), entnom), ibid)
                    if (ibid .eq. 0) then
                        valk(1) = entnom
                        valk(2) = motcl(idonn)
                        call utmess('F', 'RUPTURE0_7', nk=2, valk=valk)
                    end if
                end do
            end if
        end do
!
!
!       ---------------------------------------------------------------
!       CONSTRUCTION DE FOND DE FISSURE
!       ---------------------------------------------------------------
!
!        SI LE MOT CLE FACTEUR EST GROUP_NO ( cas 2D)
!        ----------------------------------------
!
        call jeexin('&&'//nompro//'.GROUP_NO', iret1)
        if (iret1 .ne. 0) then
            call jeveuo(noma//'.DIME', 'L', iret1)
!!          LE MOT-CLE GROUP_NO EST UNIQUEMENT AUTORISE EN DIMENSION 2
            if (zi(iret1-1+6) .eq. 2) then
                call fonnoe2(resu, noma, nompro, nbnoff, typmp)
            else
                call utmess('F', 'RUPTURE1_9')
            end if
        end if
!
!        SI LE MOT CLE FACTEUR EST GROUP_MA (cas 3D)
!        ----------------------------------------
!
        call jeexin('&&'//nompro//'.GROUP_MA', iret1)
        if (iret1 .ne. 0) then
            call jeveuo(noma//'.DIME', 'L', iret1)
!!          LE MOT-CLE GROUP_MA EST UNIQUEMENT AUTORISE EN DIMENSION 3
            if (zi(iret1-1+6) .eq. 3) then
                call fonmai2(resu, noma, typfon, iocc, nbnoff, &
                             typm)
            else
                call utmess('F', 'RUPTURE1_8')
            end if
        end if
!
!
!       DESTRUCTION DES VECTEURS DE TRAVAIL
!       ----------------------------------------
        do idonn = 1, ndonn
            call jeexin('&&'//nompro//'.'//motcl(idonn), iret)
            if (iret .ne. 0) call jedetr('&&'//nompro//'.'//motcl(idonn))
        end do
!
    end do
!
!
!     ---------------------------------------------------------------
!     VERIFICATION DES DONNEES SUR LES LEVRES ET LES VECTEURS
!     ---------------------------------------------------------------
!
!
!     TRAITEMENT DES LEVRES: LEVRE_SUP ET LEVRE_INF
!     ----------------------------------------
!
    call fonlev(resu, noma, nbnoff)
!
!
!     TRAITEMENT DU MOT-CLE NORMALE
!     ----------------------------------------
!
    call getvr8(' ', 'NORMALE', nbval=0, nbret=nvenor)
    if (nvenor .ne. 0) then
        nvenor = -nvenor
        call wkvect(resu//'.NORMALE', 'G V R8', 3, jnorm)
        call getvr8(' ', 'NORMALE', nbval=3, vect=zr(jnorm), nbret=nvenor)
    end if
!   OBJET CONTENANT LA BASE LOCALE EN CHAQUE NOEUD DU FOND DE FISSURE
!   OBJET QUI N'EXISTE QUE DANS DEFI_FOND_FISS
    basnof = resu//'.BASNOF'
!
    call fonnor2(resu, noma, cnxinv, typm, basnof)
!
    call jedetr(cnxinv)
!
!     ---------------------------------------------------------------
!     CREATION DU VECTEUR .ABSFOND CONTENANT LES ABSCISSES
!     CURVILIGNES DES NOEUDS DU FOND DE FISSURE
!     ---------------------------------------------------------------
!
!     VECTEUR CONTENANT LES NOMS DES NOEUDS DU FOND DE FISSURE
!     ----------------------------------------
    call jeexin(resu//'.FOND.NOEU', iret)
    if (iret .ne. 0) then
        fonoeu = resu//'.FOND.NOEU'
        call jeveuo(fonoeu, 'L', ifonoe)
        if (typfon .eq. 'FERME') then
            ASSERT(zk8(ifonoe+1-1) .eq. zk8(ifonoe+nbnoff-1))
        end if
    else
        ASSERT(.FALSE.)
    end if
!
!   CALCUL DE L'ABSCISSE CURVILIGNE ET DES COORDONNES DES NOEUDS DU FOND
    absfon = resu//'.ABSFON'
    coorfond = resu//'.COORFOND'
!
    call fonfis2(noma, nbnoff, fonoeu, absfon, coorfond)
!
!     ---------------------------------------------------------------
!     CREATION DE LA BASE LOCALE, DES LEVEL SETS ET DE L'ABCISSE
!     CURVILIGNE EN CHAQUE NOEUD
!     ---------------------------------------------------------------
!
!     LA BASE LOCALE ET DES LEVEL SETS SONT CALCULEES EN CHAQUE NOEUD
!     QUE SI L'OBJET BASENOF EXISTE DEJA
    call jeexin(basnof, ibas)
    if (ibas .ne. 0) then
        abscur = resu//'.ABSCUR'
        basloc = resu//'.BASLOC'
        lnno = resu//'.LNNO'
        ltno = resu//'.LTNO'
        call fonbas2(noma, basnof, typm, fonoeu, coorfond, &
                     nbnoff, absfon, basloc, abscur, lnno, &
                     ltno)
    end if
!
!     ---------------------------------------------------------------
!     EXTRACTION DES NOEUDS DES LEVRES SUR DIRECTON NORMALE
!     ---------------------------------------------------------------
!
    call getvtx(' ', 'CONFIG_INIT', scal=confin, nbret=ibid)
    if (confin .eq. 'COLLEE') then
        call jeexin(resu//'.LEVRESUP.MAIL', irets)
!
        if (irets .ne. 0) then
!
            call fonnof2(resu, noma, typfon, nbnoff, basnof)
        end if
    end if
!
!     ---------------------------------------------------------------
!     STOCKAGE D'INFOS UTILES DANS LA SD EN SORTIE
!     1- SYMETRIE (OUI, NON)
!     2- CONFIG_INIT (COLLE, DECOLLE)
!     3- TYPE_FOND (OUVERT, FERME)
!     4- NOM DU MAILLAGE
!     5- TYPE DE MAILLE EN FOND DE FISSURE (SEG2 ou SEG3)
!     6- NOM DU GROUPE DE MAILLE POUR LA LEVRE SUP
!     7- NOM DU GROUPE DE MAILLE POUR LA LEVRE INF
!     ---------------------------------------------------------------
!
    call foninf2(resu, typm, typfon, noma)
!
!     ---------------------------------------------------------------
!     IMPRESSIONS SI INFO=2
!     ---------------------------------------------------------------
!
    if (niv .eq. 2) then
        call fonimp(resu)
    end if
!
    call jedetr('&&OP0060.COORFOND')
!
    call jedema()
end subroutine
