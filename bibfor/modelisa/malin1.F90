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

subroutine malin1(motfaz, chargz, iocc, indmot, lisnoz, &
                  lonlis)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*) :: motfaz, chargz, lisnoz
!
!     CREATION DU VECTEUR DE K8 DE NOM LISNOZ ET DE LONGUEUR
!     LONLIS.
!     CE VECTEUR CONTIENT LA LISTE DES NOMS DES NOEUDS DEFINIS
!     PAR LES MOTS-CLES : GROUP_MA OU MAILLE
!     APRES LE MOT-FACTEUR LIAISON_ELEM.
!     CETTE LISTE NE CONTIENT QU'UNE OCCURENCE DES NOEUDS.
!
! IN       : MOTFAZ : MOT-CLE FACTEUR 'LIAISON_ELEM'
! IN       : CHARGZ : NOM D'UNE SD CHARGE
! IN       : IOCC   : NUMERO D'OCCURENCE DU MOT-FACTEUR
! IN       : INDMOT : INDICE = 0 --> TRAITEMENT DES MOTS-CLES
!                                    'GROUP_MA' OU 'MAILLE'
!                            = 1 --> TRAITEMENT DES MOTS-CLES
!                                     'GROUP_MA_1' OU 'MAILLE_1'
!                            = 2 --> TRAITEMENT DES MOTS-CLES
!                                     'GROUP_MA_2' OU 'MAILLE_2
! OUT      : LISNOZ : NOM DE LA LISTE DES NOEUDS
! OUT      : LONLIS : LONGUEUR DE LA LISTE DES NOEUDS
! ----------------------------------------------------------------------
!
    character(len=8) :: charge
    character(len=8) :: noma, nomnoe, nomail
    character(len=16) :: momail, mogrma
    character(len=16) :: motfac
    character(len=24) :: grmama, lisnoe
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, idim1, idim2, idimax, igr, ima
    integer(kind=8) :: in1, indlis, indmot, indnoe, ino, iocc, jdes
    integer(kind=8) :: jgro, jlist, lonlis, m, ier
    integer(kind=8) :: n1, n2, nbma, nbmail, ng, ngr, nliai
    integer(kind=8) :: nmai, numail
    character(len=24), pointer :: trav1(:) => null()
    character(len=8), pointer :: trav2(:) => null()
    integer(kind=8), pointer :: trav3(:) => null()
    aster_logical :: lcolle, lcolle2
!-----------------------------------------------------------------------
    call jemarq()
    charge = chargz
    motfac = motfaz
    lisnoe = lisnoz
!
    if (indmot .eq. 0) then
        momail = 'MAILLE'
        mogrma = 'GROUP_MA'
    else if (indmot .eq. 1) then
        momail = 'MAILLE_1'
        mogrma = 'GROUP_MA_1'
    else if (indmot .eq. 2) then
        momail = 'MAILLE_2'
        mogrma = 'GROUP_MA_2'
    end if
!
    call getfac(motfac, nliai)
    if (nliai .eq. 0) goto 999
!
    call dismoi('NOM_MAILLA', charge, 'CHARGE', repk=noma)
    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    lcolle2 = .false.
    call jeexin(noma//'.NOMMAI', ier)
    if (ier .ne. 0) then
        lcolle2 = .true.
    end if
!
    grmama = noma//'.GROUPEMA'
!
    idimax = 0
    idim1 = 0
    idim2 = 0
!
!     -- CALCUL DE IDIM1=NB_NOEUD/MAILLE*NB_MAILLE/GROUP_MA*NB_GROUP_MA
!        ET VERIFICATION DE L'APPARTENANCE DES GROUP_MA
!        AUX GROUP_MA DU MAILLAGE
!        -------------------------------------------------------
    call getvtx(motfac, mogrma, iocc=iocc, nbval=0, nbret=ng)
    if (ng .ne. 0) then
        ng = -ng
        AS_ALLOCATE(vk24=trav1, size=ng)
        call getvem(noma, 'GROUP_MA', motfac, mogrma, iocc, &
                    ng, trav1, ngr)
        do igr = 1, ngr
            call jeveuo(jexnom(grmama, trav1(igr)), 'L', jgro)
            call jelira(jexnom(grmama, trav1(igr)), 'LONUTI', nbmail)
            do m = 1, nbmail
                numail = zi(jgro-1+m)
                nomail = int_to_char8(numail, lcolle2, noma, 'MAILLE')
                ibid = char8_to_int(nomail, lcolle2, noma, 'MAILLE')
                call jelira(jexnum(noma//'.CONNEX', ibid), 'LONMAX', n1)
                idim1 = idim1+n1
            end do
        end do
    end if
!
!     -- CALCUL DE IDIM2=NB_NOEUD/MAILLE*NB_MAILLE DE LISTE DE MAILLES
!        ET VERIFICATION DE L'APPARTENANCE DES MAILLES
!        AUX MAILLES DU MAILLAGE
!        -------------------------------------------------------
    call getvtx(motfac, momail, iocc=iocc, nbval=0, nbret=nbma)
    if (nbma .ne. 0) then
        nbma = -nbma
        AS_ALLOCATE(vk8=trav2, size=nbma)
        call getvem(noma, 'MAILLE', motfac, momail, iocc, &
                    nbma, trav2, nmai)
        do ima = 1, nmai
            ibid = char8_to_int(trav2(ima), lcolle2, noma, 'MAILLE')
            call jelira(jexnum(noma//'.CONNEX', ibid), 'LONMAX', n2)
            idim2 = idim2+n2
        end do
    end if
!
!     -- IDIMAX = MAJORANT DE LA LONGUEUR DE LA LISTE DE NOEUDS
!    ----------------------------------------------------------
    idimax = idim1+idim2
!
!     -- ALLOCATION DU TABLEAU DES NOMS DE NOEUDS
!    ----------------------------------------------
    call wkvect(lisnoe, 'V V K8', idimax, jlist)
!
    indnoe = 0
!
    call getvtx(motfac, mogrma, iocc=iocc, nbval=0, nbret=ng)
    if (ng .ne. 0) then
        ng = -ng
        call getvtx(motfac, mogrma, iocc=iocc, nbval=ng, vect=trav1, &
                    nbret=ngr)
        do igr = 1, ngr
            call jeveuo(jexnom(grmama, trav1(igr)), 'L', jgro)
            call jelira(jexnom(grmama, trav1(igr)), 'LONUTI', nbmail)
            do m = 1, nbmail
                numail = zi(jgro-1+m)
                nomail = int_to_char8(numail, lcolle2, noma, 'MAILLE')
                ibid = char8_to_int(nomail, lcolle2, noma, 'MAILLE')
                call jeveuo(jexnum(noma//'.CONNEX', ibid), 'L', jdes)
                call jelira(jexnum(noma//'.CONNEX', ibid), 'LONMAX', n1)
                do ino = 1, n1
                    nomnoe = int_to_char8(zi(jdes+ino-1), lcolle, noma, 'NOEUD')
                    indnoe = indnoe+1
                    zk8(jlist+indnoe-1) = nomnoe
                end do
            end do
        end do
    end if
!
    call getvtx(motfac, momail, iocc=iocc, nbval=0, nbret=nbma)
    if (nbma .ne. 0) then
        nbma = -nbma
        call getvtx(motfac, momail, iocc=iocc, nbval=nbma, vect=trav2, &
                    nbret=nmai)
        do ima = 1, nmai
            ibid = char8_to_int(trav2(ima), lcolle2, noma, 'MAILLE')
            call jeveuo(jexnum(noma//'.CONNEX', ibid), 'L', jdes)
            ibid = char8_to_int(trav2(ima), lcolle2, noma, 'MAILLE')
            call jelira(jexnum(noma//'.CONNEX', ibid), 'LONMAX', n2)
            do ino = 1, n2
                nomnoe = int_to_char8(zi(jdes+ino-1), lcolle, noma, 'NOEUD')
                indnoe = indnoe+1
                zk8(jlist+indnoe-1) = nomnoe
            end do
        end do
    end if
!
!     -- ELIMINATION DES REDONDANCES EVENTUELLES DES NOEUDS
!        DE LA LISTE
!    -------------------------------------------------------------
    AS_ALLOCATE(vi=trav3, size=idimax)
!
    do ino = 1, idimax
        do in1 = ino+1, idimax
            if (zk8(jlist+in1-1) .eq. zk8(jlist+ino-1)) then
                trav3(in1) = 1
            end if
        end do
    end do
!
    indlis = 0
!
    do ino = 1, idimax
        if (trav3(ino) .eq. 0) then
            indlis = indlis+1
            zk8(jlist+indlis-1) = zk8(jlist+ino-1)
        end if
    end do
!
    lonlis = indlis
!
    AS_DEALLOCATE(vk24=trav1)
    AS_DEALLOCATE(vk8=trav2)
    AS_DEALLOCATE(vi=trav3)
!
999 continue
    call jedema()
end subroutine
