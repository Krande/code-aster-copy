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
subroutine elimun(noma, nomo, motfac, nzocu, nbgdcu, &
                  compcu, nopono, nolino, lisnoe, poinoe, &
                  nnoco)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/exiscp.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/palino.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: noma, nomo
    character(len=16) :: motfac
    integer(kind=8) :: nzocu
    character(len=24) :: nbgdcu
    character(len=24) :: compcu
    character(len=24) :: nopono
    character(len=24) :: nolino
    character(len=24) :: lisnoe
    character(len=24) :: poinoe
    integer(kind=8) :: nnoco
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (LIAISON_UNILATER - LECTURE)
!
! ELIMINATION AU SEIN DE CHAQUE SURFACE DE CONTACT POTENTIELLE DES
! NOEUDS ET MAILLES REDONDANTS. MODIFICATION DES POINTEURS ASSOCIES.
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NOMO   : NOM DU MODELE
! IN  MOTFAC : MOT-CLEF FACTEUR POUR LIAISON UNILATERALE
! IN  NZOCU  : NOMBRE DE ZONES DE LIAISON_UNILATERALE
! IN  NBGDCU : NOM JEVEUX DE LA SD INFOS POINTEURS GRANDEURS
! IN  COMPCU : NOM JEVEUX DE LA SD CONTENANT LES GRANDEURS DU MEMBRE
!              DE GAUCHE
! IN  NOPONO : NOM DE L'OBJET CONTENANT LE VECTEUR D'INDIRECTION
! IN  NOLINO : NOM DE L'OBJET CONTENANT LA LISTE DES NOEUDS
! IN  POINOE : NOM DE L'OBJET CONTENANT LE VECTEUR D'INDIRECTION
!               DES NOEUDS APRES NETTOYAGE
! IN  LISNOE : NOM DE L'OBJET CONTENANT LES NOEUDS APRES NETTOYAGE
! IN  NBNOE  : NOMBRE DE NOEUDS DANS LA LISTE RESULTANTE
!                VAUT NBNOE = NBTOT-NBSUP
! I/O NNOCO  : NOMBRE DE TOTAL DE NOEUDS POUR TOUTES LES OCCURRENCES
!
!
!
!
    integer(kind=8) :: jdebut, juelim, jdecal, jdecat
    integer(kind=8) :: nbelim
    character(len=8) :: k8bla, cmp, nomnoe
    integer(kind=8) :: i, j, icmp, izone, ino, numno1, numno2
    integer(kind=8) :: nbno, nbsup, nb, nbcmp, ntsup
    integer(kind=8) :: n1, n2, n3
    integer(kind=8) :: jnl, jnp, jpoi, jnoe
    character(len=24) :: nelim
    integer(kind=8) :: jelim
    integer(kind=8) :: jnbgd, jcmpg
    integer(kind=8) :: exist(1)
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES SD
!
    call jeveuo(nolino, 'E', jnl)
    call jeveuo(nopono, 'L', jnp)
    call jeveuo(nbgdcu, 'L', jnbgd)
    call jeveuo(compcu, 'L', jcmpg)
!
    ntsup = 0
    k8bla = ' '
    nelim = '&&ELIMUN.ELIM'
!
! --- CREATION DU POINTEUR
!
    call wkvect(poinoe, 'V V I', nzocu+1, jpoi)
    zi(jpoi) = 1
!
    do izone = 1, nzocu
!
! --- VECTEUR CONTENANT LES NOEUDS ZONE
!
        nbno = zi(jnp+izone)-zi(jnp+izone-1)
        jdebut = zi(jnp+izone-1)
        n1 = 0
        n2 = 0
        n3 = 0
!
! --- ELIMINATION DES PURS DOUBLONS
!
        do i = 1, nbno
            numno1 = zi(jnl-2+jdebut+i)
            if (numno1 .ne. 0) then
                do j = i+1, nbno
                    numno2 = zi(jnl-2+jdebut+j)
                    if ((numno1 .eq. numno2) .and. (numno2 .ne. 0)) then
                        zi(jnl-2+jdebut+j) = 0
                        n1 = n1+1
                    end if
                end do
            end if
        end do
!
! --- RECUPERATION INFOS SANS_NOEUD, SANS_GROUP_NO
!
        call palino(noma, motfac, 'SANS_GROUP_NO', 'SANS_NOEUD', izone, &
                    nelim)
        call jeveuo(nelim, 'L', juelim)
        nbelim = zi(juelim)
!
! --- ELIMINATION DES SANS_GROUP_NO, SANS_NOEUD
!
        call jelira(nelim, 'LONMAX', nbelim)
        call jeveuo(nelim, 'E', jelim)
        do i = 1, nbelim
            numno1 = zi(juelim-1+i)
            if (numno1 .ne. 0) then
                do j = 1, nbno
                    numno2 = zi(jnl-2+jdebut+j)
                    if ((numno1 .eq. numno2) .and. (numno2 .ne. 0)) then
                        zi(jnl-2+jdebut+j) = 0
                        n2 = n2+1
                    end if
                end do
            end if
        end do
!
! --- ELIMINATION DES NOEUDS NE COMPORTANT AUCUNE DES GRANDEURS
!
        nbcmp = zi(jnbgd+izone)-zi(jnbgd+izone-1)
        jdecat = zi(jnbgd+izone-1)
!
        do ino = 1, nbno
            numno1 = zi(jnl-2+jdebut+ino)
            if (numno1 .ne. 0) then
                nomnoe = int_to_char8(numno1)
                nb = 0
                do icmp = 1, nbcmp
!
                    cmp = zk8(jcmpg-1+jdecat+icmp-1)
                    call exiscp(cmp, k8bla, nomo, 1, 'NUM', &
                                k8bla, [numno1], exist)
                    if (exist(1) .eq. 0) then
                        nb = nb+1
                    end if
                end do
                if (nb .eq. nbcmp) then
                    zi(jnl-2+jdebut+ino) = 0
                    n3 = n3+1
                end if
            end if
        end do
!
! --- NOMBRE DE NOEUDS A SUPPRIMER
!
        nbsup = n1+n2+n3
        ntsup = ntsup+nbsup
!
! --- MAJ VECTEUR POINTEUR INDIRECT (POINOE)
!
        zi(jpoi+izone) = zi(jpoi+izone-1)+nbno-nbsup
        if (nbno .eq. nbsup) then
            call utmess('F', 'UNILATER_48')
        end if
    end do
!
! --- CREATION DU VECTEUR RESULTANT
!
    nnoco = nnoco-ntsup
    call wkvect(lisnoe, 'V V I', nnoco, jnoe)
!
! --- ELIMINATION EFFECTIVE DES NOEUDS
!
    jdecal = 0
    do izone = 1, nzocu
        nbno = zi(jnp+izone)-zi(jnp+izone-1)
        jdebut = zi(jnp+izone-1)
        do ino = 1, nbno
            numno1 = zi(jnl-2+jdebut+ino)
            if (numno1 .ne. 0) then
                zi(jnoe+jdecal) = numno1
                jdecal = jdecal+1
            end if
        end do
    end do
!
    call jedetr(nelim)
!
    call jedema()
!
end subroutine
