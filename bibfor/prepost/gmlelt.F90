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
subroutine gmlelt(igmsh, maxnod, nbtyma, nbmail, nbnoma, &
                  nuconn, versio, nbgrou)
! aslint: disable=
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/iunifi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: igmsh, maxnod, nbtyma, nbmail, nbnoma(nbtyma), nuconn(19, 32)
    integer(kind=8) :: versio
    integer(kind=8), intent(inout) :: nbgrou
!
!      GMLELT --   LECTURE DES NUMEROS DES ELEMENTS, DE LEUR TYPE,
!                  DE LEUR NUMERO DE GROUPE, DU NOMBRE DE LEURS
!                  CONNECTIVITES ET DE LEURS CONNECTIVITES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    IGMSH          IN    I         UNITE LOGIQUE DU FICHIER GMSH
!    MAXNOD         IN    I         NOMBRE MAXIMUM DE NOEUDS POUR
!                                   UNE MAILLE DONNEE
!    NBTYMA         IN    I         NOMBRE  DE TYPES DE MAILLES
!    NBMAIL         OUT   I         NOMBRE TOTAL DE MAILLES
!    NBNOMA         IN    I         NOMBRE DE NOEUDS DE LA MAILLE
!                                    POUR UN TYPE DE MAILLE DONNEE
!    NUCONN         IN    I         PASSAGE DE LA NUMEROTATION DES NDS
!                                     D'UNE MAILLE : ASTER -> GMSH
!    VERSIO         IN    I         VERSION DU FICHIER GMSH
!    NBGROU         IN    I      NUMBER OF GROUPS
!
! ......................................................................
!
!
!
!
    character(len=8) :: k8bid
    integer(kind=8) :: imes, nbmxte, nbtag, i, ij, k, icurgr
    integer(kind=8) :: indgro, ima, ibid, ityp, ino, node
    integer(kind=8) :: jnuma, jtypma, jnbnma, jnoma, jnbmag, jnbtym, jdime
    integer(kind=8) :: jindma, jtag, jgr, total_length, nb_nonempty_groups
    aster_logical :: same_dim
!
    parameter(nbmxte=19)
    integer(kind=8) :: nbno(nbmxte)
    integer(kind=8) :: dime(nbmxte)
    integer(kind=8), pointer :: noeuds(:) => null()
    data nbno/2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27,&
     &                      18, 14, 1, 8, 20, 15, 13/
    data dime/1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3,&
     &                      3, 3, 0, 2, 3, 3, 3/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATION :
!     --------------
    k8bid = '        '
!
! --- RECUPERATION DES NUMEROS D'UNITE LOGIQUE :
!     ----------------------------------------
    imes = iunifi('MESSAGE')
!
! --- LECTURE DU NOMBRE D'ELEMENTS :
!     ----------------------------
    read (igmsh, '(I10)') nbmail
!
! --- CREATION DE VECTEURS DE TRAVAIL :
!     -------------------------------
    call jedetr('&&PREGMS.NUMERO.MAILLES')
    call jedetr('&&PREGMS.TYPE.MAILLES')
    call jedetr('&&PREGMS.GROUPE.MAILLES')
    call jedetr('&&PREGMS.NBNO.MAILLES')
    call jedetr('&&PREGMS.CONNEC.MAILLES')
    call jedetr('&&PREGMS.NBMA.GROUP_MA')
    call jedetr('&&PREGMS.NBTYP.MAILLES')
    call jedetr('&&PREGMS.LISTE.GROUP_MA')
    call jedetr('&&PREGMS.TAGS')
!
! ---   VECTEUR DES NUMEROS DES MAILLES
    call wkvect('&&PREGMS.NUMERO.MAILLES', 'V V I', nbmail, jnuma)
! ---   VECTEUR DU TYPE DES MAILLES
    call wkvect('&&PREGMS.TYPE.MAILLES', 'V V I', nbmail, jtypma)
! ---   VECTEUR DU NOMBRE DE CONNECTIVITES DES MAILLES
    call wkvect('&&PREGMS.NBNO.MAILLES', 'V V I', nbmail, jnbnma)
! ---   VECTEUR DES CONNECTIVITES DES MAILLES
    call wkvect('&&PREGMS.CONNEC.MAILLES', 'V V I', maxnod*nbmail, jnoma)
! ---   VECTEUR DU NOMBRE DE MAILLES POUR UN GROUPE DE MAILLES
    if (nbgrou > 0) then
        call wkvect('&&PREGMS.NBMA.GROUP_MA', 'V V I', nbgrou, jnbmag)
        call jeveuo('&&PREGMS.NUMERO.GROUP_MA', 'L', jindma)
        if (versio == 2) call jeveuo('&&PREGMS.DIME.GROUP_MA', 'L', jdime)
    end if
! ---   VECTEUR DU NOMBRE DE MAILLES PAR TYPE DE MAILLES
    call wkvect('&&PREGMS.NBTYP.MAILLES', 'V V I', nbtyma, jnbtym)
! --- INDICATION DE DESTRUCTION DES NOEUDS
    call jeveuo('&&PREGMS.DETR.NOEUDS', 'E', vi=noeuds)
! --- TAGS POUR LE FORMAT VERSION 2 :
! --- COLLECTION DE NBMAIL VECTEURS (UN VECTEUR PAR MAILLE)
    call jecrec('&&PREGMS.TAGS', 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbmail)
!
! --- LECTURE DES ENREGISTREMENTS RELATIFS AUX MAILLES ET AFFECTATION
! --- DES VECTEURS DE TRAVAIL :
!     -----------------------
    k = 0
    ij = 0
!
! --- ICURGR : NUMERO DU GROUPE GMSH
!     NBGROU : NBRE DE GROUPES TROUVES
!     INDGRO : INDICE DU GROUPE
    icurgr = 0
    indgro = 0
    total_length = 0
    nb_nonempty_groups = 0
    do ima = 1, nbmail
        if (versio .eq. 1) then
            nbtag = 1
            !
        else if (versio .eq. 2) then
            read (igmsh, *) ibid, ibid, nbtag
            backspace (igmsh)
            !
        else
            ASSERT(.false.)
        end if
!
        call jecroc(jexnum('&&PREGMS.TAGS', ima))
        call jeecra(jexnum('&&PREGMS.TAGS', ima), 'LONMAX', nbtag)
        call jeveuo(jexnum('&&PREGMS.TAGS', ima), 'E', jtag)
!
        if (versio .eq. 1) then
!
            read (igmsh, *) zi(jnuma+ima-1), zi(jtypma+ima-1), zi(jtag), &
                ibid, zi(jnbnma+ima-1), (zi(jnoma+ij+k-1), k=1, &
                                         zi(jnbnma+ima-1))
!
        else if (versio .eq. 2) then
!
            read (igmsh, *) zi(jnuma+ima-1), zi(jtypma+ima-1), nbtag, &
                (zi(jtag-1+k), k=1, nbtag), &
                (zi(jnoma+ij+k-1), k=1, nbno(zi(jtypma+ima-1)))
!
            zi(jnbnma+ima-1) = nbno(zi(jtypma+ima-1))
!
        else
            ASSERT(.false.)
        end if
!
!      INDICATION DES NOEUDS QUI NE SONT PAS ORPHELINS
        ityp = zi(jtypma+ima-1)
        do ino = 1, nbnoma(ityp)
            node = zi(jnoma+ij+nuconn(ityp, ino)-1)
            noeuds(node+1) = 1
        end do
!
!       RECHERCHE DE L'INDICE DU GROUPE DE LA MAILLE
!       SI LE GROUPE N'EST PAS DANS PHYSICAL GROUPS, IL EST IGNORE
!       LA DIMENSION DES MAILLES DOIT CORRESPONDRE
        do i = 1, nbgrou
            do k = 1, nbtag
                if (zi(jtag+k-1) .ne. zi(jindma+i-1)) then
                    cycle
                end if
                same_dim = ASTER_TRUE
                if (versio == 2) then
                    same_dim = ASTER_FALSE
                    if (dime(ityp) .eq. zi(jdime-1+i)) same_dim = ASTER_TRUE
                end if
                if (same_dim) then
                    if (zi(jnbmag+i-1) == 0) nb_nonempty_groups = nb_nonempty_groups+1
                    zi(jnbmag+i-1) = zi(jnbmag+i-1)+1
                    total_length = total_length+1
                    ! Exit because a tag may be present twice in GMSH
                    exit
                end if
            end do
        end do
        ij = ij+zi(jnbnma+ima-1)
        zi(jnbtym+zi(jtypma+ima-1)-1) = zi(jnbtym+zi(jtypma+ima-1)-1)+1
!
    end do
!
    if (nbgrou .ne. 0) then
!
! --- CREATION DE LA COLLECTION DES GROUPES DE MAILLES :
!     ------------------------------------------------
!
!       LES GROUPES SONT SUPPOSES DISJOINTS: LONGUEUR TOTALE MAX NBMAIL
        call jecrec('&&PREGMS.LISTE.GROUP_MA', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                    nbgrou)
        call jeecra('&&PREGMS.LISTE.GROUP_MA', 'LONT', total_length+nbgrou-nb_nonempty_groups)
!
        do i = 1, nbgrou
            call jeecra(jexnum('&&PREGMS.LISTE.GROUP_MA', i), 'LONMAX', max(zi(jnbmag+i-1), 1))
            zi(jnbmag+i-1) = 0
        end do
!
! --- AFFECTATION DES OBJETS RELATIFS AUX GROUPES DE MAILLES :
!     ------------------------------------------------------
        k = 0
! --- ICURGR : NUMERO DU GROUPE GMSH
!     NBGROU : NBRE DE GROUPES TROUVES
!     INDGRO : INDICE DU GROUPE
        icurgr = 0
        indgro = 0
        do ima = 1, nbmail
            !
            call jelira(jexnum('&&PREGMS.TAGS', ima), 'LONMAX', nbtag)
            call jeveuo(jexnum('&&PREGMS.TAGS', ima), 'L', jtag)
            ityp = zi(jtypma+ima-1)
            !
            ! RECHERCHE DU GROUPE AUQUEL LA MAILLE APPARTIENT
            do i = 1, nbgrou
                call jeveuo(jexnum('&&PREGMS.LISTE.GROUP_MA', i), 'E', jgr)
                do k = 1, nbtag
                    if (zi(jtag+k-1) .ne. zi(jindma+i-1)) then
                        cycle
                    end if
                    same_dim = ASTER_TRUE
                    if (versio == 2) then
                        same_dim = ASTER_FALSE
                        if (dime(ityp) .eq. zi(jdime-1+i)) same_dim = ASTER_TRUE
                    end if
                    if (same_dim) then
                        zi(jnbmag+i-1) = zi(jnbmag+i-1)+1
                        zi(jgr+zi(jnbmag+i-1)-1) = zi(jnuma+ima-1)
                        ! Exit because a tag may be present twice in GMSH
                        exit
                    end if
                end do
            end do
!
        end do
!
    end if
!
    write (imes, *) 'NOMBRE DE MAILLES : ', nbmail
!
    call jedema()
!
! ============================ FIN DE LA ROUTINE ======================
end subroutine
