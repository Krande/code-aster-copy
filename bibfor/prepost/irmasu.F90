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
subroutine irmasu(ifc, ndim, nno, coordo, nbma, &
                  connex, point, typma, typel, codgra, &
                  codphy, codphd, permut, maxnod, lmod, &
                  noma, nbgrn, nogn, nbgrm, nogm, &
                  lmasu, nomai, nonoe, versio)
! aslint: disable=W1504
    implicit none
!
!     BUT :   ECRITURE DU MAILLAGE SUR FICHIER UNIVERSEL SUPERTAB
!     ENTREE:
!       NOMA   : NOM DU MAILLAGE
!       IFC    : NUMERO D'UNITE LOGIQUE DU FICHIER UNIVERSEL
!       NDIM   : DIMENSION DU PROBLEME (2  OU 3)
!       NNO    : NOMBRE DE NOEUDS DU MAILLAGE
!       COORDO : VECTEUR DES COORDONNEES DES NOEUDS
!       NBMA   : NOMBRE DE MAILLES DU MAILLAGE
!       CONNEX : CONNECTIVITES
!       POINT  : POINTEUR DANS LES CONNECTIVITES
!       TYPMA  : TYPES DES MAILLES
!       TYPEL  : TYPES DES ELEMENTS
!       CODGRA : CODE GRAPHIQUE SUPERTAB ASSOCIE A CHAQUE TYPE DE MAILLE
!       CODPHY : CODE PHYSIQUE SUPERTAB ASSOCIE A CHAQUE TYPE D'ELEMENT
!       PERMUT : TABLEAU DES PERMUTATIONS DES NOEUDS DE CHAQUE MAILLE
!       MAXNOD : NOMBRE MAXI DE NOEUDS POUR UN TYPE_MAILLE
!       LMOD   : LOGIQUE INDIQUANT SI IMPRESSION MODELE OU MAILLAGE
!                 .TRUE. MODELE
!       VERSIO : NUMERO DE VERSION SUPERTAB 4 OU 5, 5 PAR DEFAUT
!                (EN VERSION 5 CREATION DU DATASET 775)
!       LMASU  : INDIQUE SI MAILLAGE ISSU IDEAS .TRUE.
!       NOMAI  : NOMS DES MAILLES
!       NONOE  : NOMS DES NOEUDS
!       TOUT CE QUI SUIT CONCERNE LES GROUPES:
!          NBGRN: NOMBRE DE GROUPES DE NOEUDS
!          NOGN : NOM DES GROUPES DE NOEUDS
!          NBGRM: NOMBRE DE GROUPES DE MAILLES
!          NOGM : NOM DES GROUPES DE MAILLES
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxliis.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: noma
!     ------------------------------------------------------------------
    character(len=8) :: nomai(*), nonoe(*), nomtm
    character(len=24) :: nogn(*), nogm(*)
    real(kind=8) :: coordo(*), r(3), bidon, rcsf
    integer(kind=8) :: maxnod, connex(*), typma(*), point(*), typel(*)
    integer(kind=8) :: nodsup(32), nodast(32), permut(maxnod, *), codgra(*), codphy(*)
    integer(kind=8) :: icodno, icodma, versio, codphd(*)
    integer(kind=8) :: itri7, iqua9, iseg4, ihex27, ipen18, it15, ipen21
    aster_logical :: lmasu, lpout, lmod
! ---------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iagrma, iagrno, ica, icc, icf, icod
    integer(kind=8) :: icod1, icod2, icol, icold, icole, icoll, icolp
    integer(kind=8) :: icolt, icou, ics, icss, icst, icsv, ier
    integer(kind=8) :: ifc, igm, igm2, ign, ima, imail, imas
    integer(kind=8) :: imat, inbo, ino, ioff, iphy, ipoin, iret
    integer(kind=8) :: isup, itseg2, itype, j, jm, jmagr, jn
    integer(kind=8) :: jnogr, k, l, m, nbgrm, nbgrn
    integer(kind=8) :: nbm, nbm2, nbma, nbn, ndim, nno, nnoe
!
!-----------------------------------------------------------------------
    call jemarq()
    lpout = .false.
    itri7 = 0
    iqua9 = 0
    iseg4 = 0
    ihex27 = 0
    ipen18 = 0
    ipen21 = 0
    it15 = 0
!
!     RECHERCHE DE LA PRESENCE DE POUTRES
!
    do imail = 1, nbma
!
        if (lmod) then
            icod = codphy(typel(imail))
        else
            icod = codphd(typma(imail))
        end if
        if (icod .eq. 21 .or. icod .eq. 22 .or. icod .eq. 23 .or. icod .eq. 24) then
            lpout = .true.
            goto 2
        end if
    end do
2   continue
!
!     ECRITURE D'UN DATASET 775 (PROPRIETES DES POUTRES) BIDON
!
    if (versio .eq. 5 .and. lpout) then
        write (ifc, '(A)') '    -1'
        write (ifc, '(A)') '   775   %PROPRIETES BIDON SECTION POUTRES'
        icst = 1
        icss = 0
        icsv = 0
        write (ifc, '(3I10)') icst, icss, icsv
        write (ifc, '(A)') 'BEAM1'
        bidon = 0.0d0
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        write (ifc, '(4(1PE13.6))') bidon, bidon, bidon, bidon
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        write (ifc, '(6(1PE13.6))') bidon, bidon, bidon, bidon, bidon, bidon
        icolp = 11
        icoll = 7
        icole = 8
        icold = 14
        icol = 1
        icolt = 10
        write (ifc, '(6I10)') icolp, icoll, icole, icold, icol, icolt
        icf = 0
        ica = 45
        ics = 1
        icc = 11
        rcsf = 1.0d0
        write (ifc, '(4I10,1PE13.6)') icf, ica, ics, icc, rcsf
        write (ifc, '(A)') '    -1'
    end if
!
!     --- ECRITURE DES COORDONNEES DES NOEUDS
!
    write (ifc, '(A)') '    -1'
    if (versio .eq. 5) then
        write (ifc, '(A)') '  2411'
    else if (versio .eq. 4) then
        write (ifc, '(A)') '    15   %NOEUDS'
    end if
!     I : NUMERO DES NOEUDS
!     J : SYSTEME DE COORDONNEES
!     K : ""           ""
!     L : COULEUR
    j = 1
    k = 1
    l = 11
    ndim = 3
!     --- BOUCLE SUR TOUS LES NOEUDS DU MAILLAGE
    do i = 1, nno
!       - ON ECRIT TOUJOURS 3 COORDONNEES (PARAMETRE NDIM ECRASE)
        do m = 1, ndim
            r(m) = coordo(3*(i-1)+m)
        end do
!       - RECUPERATION DU NUMERO DE NOEUDS (RECHERCHE SI MAILLAGE IDEAS)
        if (lmasu) then
            call lxliis(nonoe(i) (1:8), ino, ier)
        else
            ino = i
        end if
        if (versio .eq. 5) then
            write (ifc, fmt='(4I10)') ino, j, k, l
            write (ifc, fmt='(3E25.17)') (r(m), m=1, ndim)
        else if (versio .eq. 4) then
            write (ifc, fmt='(4I10,3E13.6)') ino, j, k, l, (r(m), m=1, ndim)
        end if
    end do
    write (ifc, '(A)') '    -1'
!
!     --- ECRITURE DES CONNECTIVITES
!
    write (ifc, '(A)') '    -1'
    if (versio .eq. 5) then
        write (ifc, '(A)') '  2412'
    else if (versio .eq. 4) then
        write (ifc, '(A)') '    71   %ELEMENTS'
    end if
    imat = 1
    icou = 6+imat
    do ima = 1, nbma
        imat = 1
        itype = typma(ima)
!       - RECUPERATION DU NOMBRE DE NOEUDS DE L'ELEMENT
        ipoin = point(ima)
        nnoe = point(ima+1)-ipoin
!CC
        call jenuno(jexnum('&CATA.TM.NOMTM', itype), nomtm)
        if (nomtm .eq. 'HEXA27') then
            if (ihex27 .eq. 0) then
                call utmess('I', 'PREPOST2_78')
            end if
            ihex27 = 1
            nnoe = nnoe-7
        else if (nomtm .eq. 'TRIA7') then
            if (itri7 .eq. 0) then
                call utmess('I', 'PREPOST2_79')
            end if
            itri7 = 1
            nnoe = nnoe-1
        else if (nomtm .eq. 'PENTA18') then
            if (ipen18 .eq. 0) then
                call utmess('I', 'PREPOST2_85')
            end if
            ipen18 = 1
            nnoe = nnoe-3
        else if (nomtm .eq. 'TETRA15') then
            if (it15 .eq. 0) then
                call utmess('I', 'PREPOST2_86')
            end if
            it15 = 1
            nnoe = nnoe-5
        else if (nomtm .eq. 'PENTA21') then
            if (ipen21 .eq. 0) then
                call utmess('I', 'PREPOST2_87')
            end if
            ipen21 = 1
            nnoe = 15
        else if (nomtm .eq. 'QUAD9') then
            if (iqua9 .eq. 0) then
                call utmess('I', 'PREPOST2_80')
            end if
            iqua9 = 1
            nnoe = nnoe-1
        else if (nomtm .eq. 'SEG4') then
            if (iseg4 .eq. 0) then
                call utmess('I', 'PREPOST2_81')
            end if
            iseg4 = 1
            nnoe = nnoe-2
!JMP
            call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), itseg2)
            itype = itseg2
        end if
!CC
!       - IMAT = 2 CORRESPOND A UN ELEMENT REDUIT A UN NOEUD
        if (nnoe .eq. 1) imat = 2
!       - PERMUTATION DES NOEUDS SELON LE TYPE DE MAILLE
        do j = 1, nnoe
            nodast(j) = connex(ipoin-1+j)
            isup = permut(j, itype)
            nodsup(isup) = nodast(j)
        end do
!       - CODE GRAPHIQUE DE LA MAILLE
        icod1 = codgra(itype)
!       - SI LMOD =.TRUE. ON IMPRIME LE MODELE SINON LE MAILLAGE
        if (lmod) then
            if (typel(ima) .eq. 0) then
                icod2 = 0
            else
                icod2 = codphy(typel(ima))
                iphy = typel(ima)
            end if
        else
            icod2 = codphd(typma(ima))
!JMPESSAI          ICOD2=CODPHD(ITYPE)
            iphy = itype
        end if
!
!       - ELEMENTS NON DISPONIBLES DANS IDEAS
        if (codphd(typma(ima)) .eq. 6000) then
            call utmess('A', 'PREPOST2_82')
        else if (codphd(typma(ima)) .eq. 6001) then
            call utmess('A', 'PREPOST2_83')
        end if
        ASSERT((icod2 .le. 10000) .and. (icod2 .ge. 0))
        if (icod1 .ne. 0 .and. icod2 .ne. 0) then
            if (lmasu) then
                call lxliis(nomai(ima) (1:8), imas, ier)
                do j = 1, nnoe
                    call lxliis(nonoe(nodsup(j)) (1:8), nodsup(j), ier)
                end do
            else
                imas = ima
            end if
            if (versio .eq. 5) then
                iphy = 1
                write (ifc, fmt='(6I10)') imas, icod2, iphy, imat, icou, &
                    nnoe
!           - ENREGISTREMENT SPECIAL POOUR LES POUTRES
                if (icod2 .eq. 21 .or. icod2 .eq. 22 .or. icod2 .eq. 23 .or. icod2 .eq. 24) then
                    inbo = 0
                    ioff = 1
                    write (ifc, fmt='(3I10)') inbo, ioff, ioff
                end if
            else if (versio .eq. 4) then
                write (ifc, fmt='(7I10)') imas, icod1, icod2, iphy, imat, &
                    icou, nnoe
            end if
            write (ifc, fmt='(8I10)') (nodsup(j), j=1, nnoe)
        end if
    end do
    write (ifc, '(A)') '    -1'
!
!     --- ECRITURE DES GROUPES DE NOEUDS ET DE MAILLES
    write (ifc, '(A)') '    -1'
    write (ifc, '(A)') '  2477'
    icodma = 8
    icodno = 7
    call jeexin('&&IRMASU.NOGR', iret)
    if (iret .ne. 0) call jedetr('&&IRMASU.NOGR')
!     --- TRAITEMENT DES GROUPES DE NOEUDS
    do ign = 1, nbgrn
        call jeveuo(jexnum(noma//'.GROUPENO', ign), 'L', iagrno)
        call jelira(jexnum(noma//'.GROUPENO', ign), 'LONUTI', nbn)
        call wkvect('&&IRMASU.NOGR', 'V V I', nbn, jnogr)
        write (ifc, fmt='(8I10)') ign, 0, 0, 0, 0, 0, 0, nbn
        write (ifc, fmt='(A)') nogn(ign)
        do jn = 1, nbn
            zi(jnogr-1+jn) = zi(iagrno-1+jn)
            if (lmasu) then
                call lxliis(nonoe(zi(jnogr-1+jn)) (1:8), zi(jnogr-1+jn), ier)
            end if
        end do
        write (ifc, fmt='(8I10)') (icodno, zi(jnogr-1+jn), 0, 0, jn=1, nbn)
        call jedetr('&&IRMASU.NOGR')
    end do
!     --- TRAITEMENT DES GROUPES DE MAILLES
    igm2 = 0
    call jeexin('&&IRMASU.MAGR', iret)
    if (iret .ne. 0) call jedetr('&&IRMASU.MAGR')
    do igm = 1, nbgrm
!       - RECUPERATION DU NOMBRE DE MAILLES DU GROUPE
!         ET DU POINTEUR SURLE GROUPE DE MAILLES
        call jeveuo(jexnum(noma//'.GROUPEMA', igm), 'L', iagrma)
        call jelira(jexnum(noma//'.GROUPEMA', igm), 'LONUTI', nbm)
        call wkvect('&&IRMASU.MAGR', 'V V I', nbm, jmagr)
        nbm2 = 0
        do jm = 1, nbm
!         - MAILLE PRISE EN COMPTE SI TYPE D'ELEMENT DIFFERENT DE 0
            if (lmod) then
                if (typel(zi(iagrma-1+jm)) .eq. 0) goto 756
            end if
            itype = typma(zi(iagrma-1+jm))
            call jenuno(jexnum('&CATA.TM.NOMTM', itype), nomtm)
            if ((nomtm .eq. 'PYRAM5') .or. (nomtm .eq. 'PYRAM13')) goto 756
            nbm2 = nbm2+1
            zi(jmagr-1+nbm2) = zi(iagrma-1+jm)
            if (lmasu) then
!           - SI MAILLAGE IDEAS RECUPERATION DU NUMERO DE MAILLE
                call lxliis(nomai(zi(jmagr-1+nbm2)) (1:8), zi(jmagr-1+nbm2), ier)
            end if
756         continue
        end do
        if (nbm2 .ne. 0) then
            igm2 = igm2+1
            write (ifc, fmt='(8I10)') nbgrn+igm2, 0, 0, 0, 0, 0, 0, nbm2
            write (ifc, fmt='(A)') nogm(igm)
            write (ifc, fmt='(8I10)') (icodma, zi(jmagr-1+jm), 0, 0, jm=1, &
                                       nbm2)
        end if
        call jedetr('&&IRMASU.MAGR')
    end do
    write (ifc, '(A)') '    -1'
    call jedema()
end subroutine
