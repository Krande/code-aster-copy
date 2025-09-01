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

subroutine irmail(form, ifi, versio, noma, lmod, &
                  ligrelz, infmai, formar, lfichUniq, nosdfu)
!
    implicit none
!
!     BUT: ECRITURE DU MAILLAGE AU FORMAT RESULTAT, IDEAS, ENSIGHT, MED,
!     ENTREE:
!        FORM     : FORMAT DES IMPRESSIONS: IDEAS, ENSIGHT, ...
!        IFI      : UNITE LOGIQUE D'IMPRESSION
!        VERSIO   : VERSION IDEAS 4 OU 5 PAR DEFAUT 5
!        NOMA     : NOM UTILISATEUR DU MAILLAGE A ECRIRE
!        LMOD     : LOGIQUE INDIQUANT SI IMPRESSION MODELE OU MAILLAGE
!                 .TRUE. MODELE
!        LIGREL   : NOM DU LIGREL DU CHAMP ' ' SI SEULEMENT MAILLAGE
!        INFMAI   : POUR LE FORMAT MED, NIVEAU DES INFORMATIONS A IMPRIMER
!        FORMAR   : FORMAT REEL OU COMPLEXE
!        LFICHUNIQ: ASTER LOGICAL, FICHIER UNIQUE
!        NOSDFU   : NOM STRUCTURE DONNEE
!     ------------------------------------------------------------------
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/iradhs.h"
#include "asterfort/irmare.h"
#include "asterfort/irmasu.h"
#include "asterfort/irmgms.h"
#include "asterfort/irmhdf.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
!---------------- ARGUMENTS --------------------------------------------
    integer(kind=8)             :: versio, infmai
    aster_logical       :: lmod
    character(len=8)    :: noma
    character(len=19)   :: ligrelz
    character(len=16)   :: formar
    character(len=*)    :: form
    aster_logical, optional :: lfichUniq
    character(len=8), optional :: nosdfu
!---------------- VARIABLES LOCALES ------------------------------------
!
    integer(kind=8) :: ifi, igm, ign
    integer(kind=8) :: ima, ino, iret
    integer(kind=8) :: jnogm, jnogn
    integer(kind=8) :: jnomai, jnonoe, jpoin
    integer(kind=8) :: jtitr, jtypl
!
    integer(kind=8) :: lon1, maxnod, nbgrm, nbgrn
    integer(kind=8) :: nbmai, nbnoe, nbtitr, ndim, repi
!
    aster_logical :: lmasu, lgmsh
!
    character(len=6) :: nopaje
    character(len=8) :: nosdf2
    character(len=11) :: nojgrp
    character(len=19) :: ligrel
    character(len=80)       :: titmai
    real(kind=8), pointer   :: vale(:) => null()
    integer(kind=8), pointer        :: connex(:) => null()
    integer(kind=8), pointer        :: codegra(:) => null()
    integer(kind=8), pointer        :: codephd(:) => null()
    integer(kind=8), pointer        :: codephy(:) => null()
    integer(kind=8), pointer        :: permuta(:) => null()
    integer(kind=8), pointer        :: typmail(:) => null()
    aster_logical :: lfu
!     ------------------------------------------------------------------
!
    call jemarq()
    if (.not. present(lfichUniq)) then
        lfu = .false._1
        nosdf2 = ' '
    else
        lfu = lfichUniq
        ASSERT(present(nosdfu))
        nosdf2 = nosdfu
    end if
!
    ligrel = ' '
    if (lmod) then
        ligrel = ligrelz
    end if
!
    if (ligrel .ne. ' ') then
        call dismoi('DIM_GEOM', ligrel, 'LIGREL', repi=repi)
!       avec repi   =   1  : 1D
!                   =   2  : 2D
!                   =   3  : 3D
!                   = 120  : 1D+2D     MELANGE
!                   = 103  : 1D+3D     MELANGE
!                   =  23  : 2D+3D     MELANGE
!                   = 123  : 1D+2D+3D  MELANGE
        select case (repi)
        case (1)
            ndim = 1
        case (2, 120)
            ndim = 2
        case (3, 103, 23, 123)
            ndim = 3
        end select
    else
        call dismoi('DIM_GEOM_B', noma, 'MAILLAGE', repi=ndim)
    end if
!
!     --- RECUPERATION DU NOMBRE DE MAILLES
    call jelira(noma//'.TYPMAIL', 'LONMAX', nbmai)
!
!     --- NBNOE = NOMBRE DE NOEUDS DU MAILLAGE (RECUPERATION VALEUR)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoe)
!
!     --- RECUPERATION DES VECTEURS COORDONNEES DES NOEUDS JCOOR
!                      DU  VECTEUR DES CONNECTIVITES
!                      DU  POINTEUR SUR LES CONNECTIVITES
!                      DU  POINTEUR SUR LES TYPES DE MAILLE
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jpoin)
    call jeveuo(noma//'.TYPMAIL        ', 'L', vi=typmail)
!
!     --- CONSTITUTION DU TITRE (SUR PLUSIEURS LIGNES EVENTUELLEMENT)
    call jeexin(noma//'           .TITR', iret)
    if (iret .gt. 0) then
        call jeveuo(noma//'           .TITR', 'L', jtitr)
        call jelira(noma//'           .TITR', 'LONMAX', nbtitr)
    else
        nbtitr = 1
        call wkvect(noma//'           .TITR', 'V V K80', nbtitr, jtitr)
        zk80(jtitr) = 'MAILLAGE RESTREINT'
    end if
!
!
!       - DESTRUCTION PUIS ALLOCATION DE ZONES DE TRAVAIL
    call jeexin('&&IRMAIL.NOMMAI', iret)
    if (iret .ne. 0) call jedetr('&&IRMAIL.NOMMAI')
    call jeexin('&&IRMAIL.NOMNOE', iret)
    if (iret .ne. 0) call jedetr('&&IRMAIL.NOMNOE')
    call wkvect('&&IRMAIL.NOMMAI', 'V V K8', nbmai, jnomai)
    call wkvect('&&IRMAIL.NOMNOE', 'V V K8', nbnoe, jnonoe)
!       - RECUPERATION DES NOMS DES MAILLES
    do ima = 1, nbmai
        zk8(jnomai-1+ima) = int_to_char8(ima)
    end do
!       - RECUPERATION DES NOMS DES NOEUDS
    do ino = 1, nbnoe
        zk8(jnonoe-1+ino) = int_to_char8(ino)
    end do
!       - TEST EXISTENCE DE GROUPES DE NOEUDS
    if (lfu) then
        nopaje = 'NOMMAX'
        nojgrp = '.PAR_GRPNOE'
    else
        nopaje = 'NUTIOC'
        nojgrp = '.GROUPENO'
    end if
    call jeexin(noma//nojgrp, iret)
    if (iret .ne. 0) then
!         - RECUPERATION DU NOMBRE ET DES NOMS DES GROUPES DE NOEUDS
        call jelira(noma//nojgrp, nopaje, nbgrn)
        if (nbgrn .ne. 0) then
            call wkvect('&&IRMAIL.NOMGRNO', 'V V K24', nbgrn, jnogn)
            do ign = 1, nbgrn
                call jenuno(jexnum(noma//nojgrp, ign), zk24(jnogn-1+ign))
            end do
        else
!           - SI PAS DE GROUPE DE NOEUDS - NOMBRE DE GROUPES = 0
            nbgrn = 0
        end if
    else
        jnogn = 1
        nbgrn = 0
    end if
!       - TEST EXISTENCE DE GROUPES DE MAILLE
    if (lfu) then
        nojgrp = '.PAR_GRPMAI'
    else
        nojgrp = '.GROUPEMA'
    end if
    call jeexin(noma//nojgrp, iret)
    if (iret .ne. 0) then
!         - RECUPERATION DU NOMBRE ET DES NOMS DES GROUPES DE MAILLES
        call jelira(noma//nojgrp, nopaje, nbgrm)
        if (nbgrm .ne. 0) then
            call wkvect('&&IRMAIL.NOMGRMA', 'V V K24', nbgrm, jnogm)
            do igm = 1, nbgrm
                call jenuno(jexnum(noma//nojgrp, igm), zk24(jnogm-1+igm))
            end do
        else
            nbgrm = 0
        end if
    else
        jnogm = 1
        nbgrm = 0
    end if
    if (lmod) then
!       - IMPRESSION DU MODELE
!         --> ON RECUPERE LE TYPE D'ELEMENT FINI DES MAILLES
        call jeveuo(ligrel//'.TYFE', 'L', jtypl)
    else
        jtypl = 1
    end if
!
    if (form .eq. 'RESULTAT') then
!       - TRAITEMENT DU FORMAT 'RESULTAT'
        call irmare(ifi, ndim, nbnoe, vale, nbmai, &
                    connex, zi(jpoin), noma, typmail, zi(jtypl), &
                    lmod, zk80(jtitr), nbtitr, nbgrn, nbgrm, &
                    zk8(jnomai), zk8(jnonoe), formar)
!
    else if (form .eq. 'ASTER') then
!       - TRAITEMENT DU FORMAT 'ASTER'
        call irmare(ifi, ndim, nbnoe, vale, nbmai, &
                    connex, zi(jpoin), noma, typmail, zi(jtypl), &
                    lmod, zk80(jtitr), nbtitr, nbgrn, nbgrm, &
                    zk8(jnomai), zk8(jnonoe), formar)
!
    else if (form .eq. 'MED') then
!       - TRAITEMENT DU FORMAT ECHANGE DE DONNEES 'MED'
        if (lfu) then
            call irmhdf(ifi, ndim, nbnoe, vale, nbmai, &
                        connex, zi(jpoin), noma, typmail, zk80(jtitr), &
                        nbtitr, nbgrn, zk24(jnogn), nbgrm, zk24(jnogm), &
                        zk8(jnomai), zk8(jnonoe), infmai, lfu, nosdf2)
        else
            call irmhdf(ifi, ndim, nbnoe, vale, nbmai, &
                        connex, zi(jpoin), noma, typmail, zk80(jtitr), &
                        nbtitr, nbgrn, zk24(jnogn), nbgrm, zk24(jnogm), &
                        zk8(jnomai), zk8(jnonoe), infmai)
        end if
!
    else if (form .eq. 'GMSH') then
!       - TRAITEMENT DU FORMAT 'GMSH'
!         ON REGARDE SI LE MAILLAGE EST UN MAILLAGE GMSH (LGMSH)
        lgmsh = .false.
        call jeexin(noma//'           .TITR', iret)
        if (iret .ne. 0) then
            call jeveuo(noma//'           .TITR', 'L', jtitr)
            call jelira(noma//'           .TITR', 'LONMAX', nbtitr)
            if (nbtitr .ge. 1) then
                titmai = zk80(jtitr-1+1)
                if (titmai(10:31) .eq. 'AUTEUR=INTERFACE_GMSH') then
                    lgmsh = .true.
                end if
            end if
        end if
        call irmgms(ifi, ndim, nbnoe, noma, nbgrm, &
                    zk8(jnonoe), lgmsh, versio)
!
    else if (form(1:5) .eq. 'IDEAS') then
!       - TRAITEMENT FORMAT 'IDEAS'
!         ON REGARDE SI LE MAILLAGE EST UN MAILLAGE SUPERTAB (LMASU)
        lmasu = .false.
        call jeexin(noma//'           .TITR', iret)
        if (iret .ne. 0) then
            call jeveuo(noma//'           .TITR', 'L', jtitr)
            call jelira(noma//'           .TITR', 'LONMAX', nbtitr)
            if (nbtitr .ge. 1) then
                titmai = zk80(jtitr-1+1)
                if (titmai(10:31) .eq. 'AUTEUR=INTERFACE_IDEAS') then
                    lmasu = .true.
                end if
            end if
        end if
!       - SOUS PROGRAMME : TRAITER LES ADHERENCES SUPERTAB
        call iradhs(versio)
        call jeveuo('&&IRADHS.CODEGRA', 'L', vi=codegra)
        call jeveuo('&&IRADHS.CODEPHY', 'L', vi=codephy)
        call jeveuo('&&IRADHS.CODEPHD', 'L', vi=codephd)
        call jeveuo('&&IRADHS.PERMUTA', 'L', vi=permuta)
        call jelira('&&IRADHS.PERMUTA', 'LONMAX', lon1)
        maxnod = permuta(lon1)
        call irmasu(ifi, ndim, nbnoe, vale, nbmai, &
                    connex, zi(jpoin), typmail, zi(jtypl), codegra, &
                    codephy, codephd, permuta, maxnod, lmod, &
                    noma, nbgrn, zk24(jnogn), nbgrm, zk24(jnogm), &
                    lmasu, zk8(jnomai), zk8(jnonoe), versio)
!       - DESTRUCTION ZONE ALLOUEE POUR GPES DE NOEUDS SI ELLE EXISTE
        call jeexin('&&IRMASU.NOMGRNO', iret)
        if (iret .ne. 0) then
            call jedetr('&&IRMASU.NOMGRNO')
        end if
!       - DESTRUCTION ZONE ALLOUEE POUR GPES DE MAILLES SI ELLE EXISTE
        call jeexin('&&IRMASU.NOMGRMA', iret)
        if (iret .ne. 0) then
            call jedetr('&&IRMASU.NOMGRMA')
        end if
        call jedetr('&&IRADHS.PERMUTA')
        call jedetr('&&IRADHS.CODEGRA')
        call jedetr('&&IRADHS.CODEPHY')
        call jedetr('&&IRADHS.CODEPHD')
        call jedetr('&&IRMAIL.NOMMAI')
        call jedetr('&&IRMAIL.NOMNOE')
!
    end if
!
! --- MENAGE
    call jedetr('&&IRMAIL.NOMMAI')
    call jedetr('&&IRMAIL.NOMNOE')
    call jedetr('&&IRMAIL.NOMGRMA')
    call jedetr('&&IRMAIL.NOMGRNO')
!
    call jedema()
end subroutine
