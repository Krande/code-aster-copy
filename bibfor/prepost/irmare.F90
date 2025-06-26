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
subroutine irmare(ifc, ndim, nno, coordo, nbma, &
                  connex, point, noma, typma, typel, &
                  lmod, titre, nbtitr, nbgrn, nbgrm, &
                  nomai, nonoe, formar)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=80) :: titre(*)
    character(len=8) :: nomai(*), nonoe(*), noma
    character(len=16) :: formar
    real(kind=8) :: coordo(*)
    integer(kind=8) :: connex(*), typma(*), point(*), typel(*), ifc, nbtitr
    aster_logical :: lmod
!
!
!     BUT :   ECRITURE DU MAILLAGE AU FORMAT ASTER
!     ENTREE:
!       IFC    : NUMERO D'UNITE LOGIQUE DU FICHIER ASTER
!       NDIM   : DIMENSION DU PROBLEME (2  OU 3)
!       NNO    : NOMBRE DE NOEUDS DU MAILLAGE
!       COORDO : VECTEUR DES COORDONNEES DES NOEUDS
!       NBMA   : NOMBRE DE MAILLES DU MAILLAGE
!       CONNEX : CONNECTIVITES
!       POINT  : POINTEUR DANS LES CONNECTIVITES
!       NOMAT  : NOM DU MAILLAGE
!       TYPMA  : TYPES DES MAILLES
!       TYPEL  : TYPES DES ELEMENTS
!       LMOD   : LOGIQUE INDIQUANT SI IMPRESSION MODELE OU MAILLAGE
!                 .TRUE. : ON N'IMPRIME QUE LES MAILLES DU MODELE
!       TITRE  : TITRE ASSOCIE AU MAILLAGE
!       TOUT CE QUI SUIT CONCERNE LES GROUPES:
!          NBGRN: NOMBRE DE GROUPES DE NOEUDS
!          NBGRM: NOMBRE DE GROUPES DE MAILLES
!          NOMAI: NOMS DES MAILLES
!          NONOE: NOMS DES NOEUDS
!
!     ------------------------------------------------------------------
! ---------------------------------------------------------------------
!
    character(len=8) :: nomm
    character(len=10) :: format
    character(len=24) :: nomgr
    character(len=50) :: fmt
!
!     ECRITURE DU TITRE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iagrma, iagrno, ico, ifin, igm, ign
    integer(kind=8) :: ima, ino, ipo, ipoin, it, itype, itypp
    integer(kind=8) :: j, jm, jmai, jn, k, nbfois, nbgrm
    integer(kind=8) :: nbgrn, nbm, nbma, nbn, nbrest, ndim, nno
    integer(kind=8) :: nnoe
!-----------------------------------------------------------------------
    call jemarq()
    format = formar
    fmt = '(1X,A8,1X,'//format//',1X,'//format//',1X,'//format//')'
    write (ifc, *) 'TITRE'
    do it = 1, nbtitr
        write (ifc, '(A)') titre(it)
    end do
    write (ifc, *) 'FINSF'
    write (ifc, *) '%'
!
!
!     ECRITURE DES NOEUDS
!     --------------------
    if (ndim .eq. 3) then
        write (ifc, *) 'COOR_3D'
    else if (ndim .eq. 2) then
        write (ifc, *) 'COOR_2D'
    else if (ndim .eq. 1) then
        write (ifc, *) 'COOR_1D'
    else
        call utmess('F', 'PREPOST2_77')
    end if
    do ino = 1, nno
        write (ifc, fmt) "N"//nonoe(ino) (1:7), (coordo(3*(ino-1)+j), j=1, ndim)
    end do
!
!
!     ECRITURE DES MAILLES
!     ---------------------
    itypp = 0
    ifin = 0
    do ima = 1, nbma
        itype = typma(ima)
        ipoin = point(ima)
        nnoe = point(ima+1)-ipoin
!
!        -- SI LMOD =.TRUE. ON IMPRIME LE MODELE SINON LE MAILLAGE
        if (lmod) then
            if (typel(ima) .eq. 0) then
                goto 21
            end if
        end if
        if (itype .ne. itypp) then
            call jenuno(jexnum('&CATA.TM.NOMTM', itype), nomm)
            write (ifc, *) 'FINSF'
            write (ifc, *) '%'
            itypp = itype
            write (ifc, *) nomm
            ifin = 1
        end if
        nbfois = nnoe/7
        nbrest = nnoe-7*nbfois
        if (nbfois .ge. 1) then
            write (ifc, 1003) "M"//nomai(ima) (1:7), &
                ("N"//nonoe(connex(ipoin-1+i)) (1:7), i=1, 7)
            ico = 8
            do i = 2, nbfois
                write (ifc, 1004) ("N"//nonoe(connex(ipoin-1+k)) (1:7), k=ico, ico+6)
                ico = ico+7
            end do
            if (nbrest .ne. 0) then
                write (ifc, 1004) ("N"//nonoe(connex(ipoin-1+i)) (1:7), i=ico, nnoe)
            end if
        else
            write (ifc, 1003) "M"//nomai(ima) (1:7), &
                ("N"//nonoe(connex(ipoin-1+i)) (1:7), i=1, nnoe)
        end if
21      continue
    end do
    if (ifin .eq. 1) then
        write (ifc, *) 'FINSF'
        write (ifc, *) '%'
    end if
!
!
!     ECRITURE DES GROUPES DE NOEUDS
!     -------------------------------
    do ign = 1, nbgrn
        call jenuno(jexnum(noma//'.GROUPENO', ign), nomgr)
        call jelira(jexnum(noma//'.GROUPENO', ign), 'LONUTI', nbn)
        write (ifc, *) 'GROUP_NO'
        write (ifc, *) nomgr
        if (nbn .gt. 0) then
            call jeveuo(jexnum(noma//'.GROUPENO', ign), 'L', iagrno)
            write (ifc, '(7(1X,A8))') ("N"//nonoe(zi(iagrno-1+jn)) (1:7), jn=1, nbn)
        end if
        write (ifc, *) 'FINSF'
        write (ifc, *) '%'
    end do
!
!
!     ECRITURE DES GROUPES DE MAILLES
!     --------------------------------
    do igm = 1, nbgrm
        call jenuno(jexnum(noma//'.GROUPEMA', igm), nomgr)
        call jelira(jexnum(noma//'.GROUPEMA', igm), 'LONUTI', nbm)
        write (ifc, *) 'GROUP_MA'
        write (ifc, *) nomgr
        if (nbm .gt. 0) then
            call jeveuo(jexnum(noma//'.GROUPEMA', igm), 'L', iagrma)
            call wkvect('&&IRMARE.NOMAI', 'V V K8', max(nbm, 1), jmai)
            ipo = 0
            do jm = 1, nbm
                if (lmod) then
                    if (typel(zi(iagrma-1+jm)) .eq. 0) goto 756
                end if
                zk8(jmai-1+ipo+1) = "M"//nomai(zi(iagrma-1+jm)) (1:7)
                ipo = ipo+1
756             continue
            end do
            if (ipo .ne. 0) then
                write (ifc, '(7(1X,A8))') (zk8(jmai-1+jm), jm=1, ipo)
            end if
            call jedetr('&&IRMARE.NOMAI')
        end if
        write (ifc, *) 'FINSF'
        write (ifc, *) '%'
    end do
!
!
    write (ifc, *) 'FIN'
!
1003 format(8(1x, a8))
1004 format(9x, 7(1x, a8))
    call jedema()
end subroutine
