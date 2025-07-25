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
subroutine gmgnre(noma, nbnoto, litrav, listma, nbma, &
                  listno, nbno, selez)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: noma
    character(len=*) :: selez
    character(len=16) :: selec
    integer(kind=8) :: nbma, nbnoto, nbno, listma(*), listno(*), litrav(*)
! ----------------------------------------------------------------------
!     BUT: REMPLIR LA LISTE DE NOEUD SOUS-JACENTE A LA LISTE DE MAILLE
!
!     IN: NOMA   : NOM DU MAILLAGE
!         NBNOTO : NOMBRE DE NOEUDS TOTAL DU MAILLAGE.
!         LISTMA : LISTE DES NUMEROS DE MAILLES A TRAITER.
!           NBMA : NOMBRE DE MAILLES DANS LA LISTE.
!         LITRAV : VECTEUR DE TRAVAIL.
!         SELEC  :  SELECTION DES NOEUDS (TOUS, SOMMET, MILIEU, CENTRE)
!
!     OUT:
!         LISTNO : LISTE DES NOEUDS TROUVES
!          NBNO  : NOMBRE DE NOEUDS TROUVE.
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: typm, notyma(19)
    integer(kind=8) :: posini, posfin, sel, nutyma
    integer(kind=8) :: pini(3, 19), pfin(3, 19)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacnex, ima, ino, nbnoma, numno
    integer(kind=8), pointer :: typmail(:) => null()
!-----------------------------------------------------------------------
    data notyma/'POI1',&
     &              'SEG2', 'SEG3',&
     &              'TRIA3', 'TRIA6', 'TRIA7',&
     &              'QUAD4', 'QUAD8', 'QUAD9',&
     &              'TETRA4', 'TETRA10',&
     &              'PENTA6', 'PENTA15', 'PENTA18',&
     &              'PYRAM5', 'PYRAM13',&
     &              'HEXA8', 'HEXA20', 'HEXA27'/
!
!
!
!
    data pini/1, 0, 0,&
     &            1, 0, 0,&
     &            1, 3, 0,&
     &            1, 0, 0,&
     &            1, 4, 0,&
     &            1, 4, 7,&
     &            1, 0, 0,&
     &            1, 5, 0,&
     &            1, 5, 9,&
     &            1, 0, 0,&
     &            1, 5, 0,&
     &            1, 0, 0,&
     &            1, 7, 0,&
     &            1, 7, 16,&
     &            1, 0, 0,&
     &            1, 6, 0,&
     &            1, 0, 0,&
     &            1, 9, 0,&
     &            1, 9, 21/
!
    data pfin/1, 0, 0,&
     &            2, 0, 0,&
     &            2, 3, 0,&
     &            3, 0, 0,&
     &            3, 6, 0,&
     &            3, 6, 7,&
     &            4, 0, 0,&
     &            4, 8, 0,&
     &            4, 8, 9,&
     &            4, 0, 0,&
     &            4, 10, 0,&
     &            6, 0, 0,&
     &            6, 15, 0,&
     &            6, 15, 18,&
     &            5, 0, 0,&
     &            5, 13, 0,&
     &            8, 0, 0,&
     &            8, 20, 0,&
     &            8, 20, 27/
!
!
!
!     -- ON PARCOURE LA LISTE DES MAILLES ET ON COCHE LES NOEUDS
!     -- DANS LITRAV:
!
    call jemarq()
    selec = selez
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
    if (selec .eq. 'TOUS') sel = 0
    if (selec .eq. 'SOMMET') sel = 1
    if (selec .eq. 'MILIEU') sel = 2
    if (selec .eq. 'CENTRE') sel = 3
!
    do i = 1, nbnoto
        litrav(i) = 0
    end do
!
    do i = 1, nbma
        ima = listma(i)
        call jeveuo(jexnum(noma//'.CONNEX', ima), 'L', iacnex)
        call jelira(jexnum(noma//'.CONNEX', ima), 'LONMAX', nbnoma)
!
        if (sel .eq. 0) then
            posini = 1
            posfin = nbnoma
        else
            call jenuno(jexnum('&CATA.TM.NOMTM', typmail(ima)), typm)
            do nutyma = 1, 18
                if (typm .eq. notyma(nutyma)) then
                    posini = pini(sel, nutyma)
                    posfin = pfin(sel, nutyma)
                    goto 20
                end if
            end do
            call utmess('F', 'MODELISA4_68', sk=typm)
20          continue
            if (posfin .eq. 0) goto 2
        end if
!
        do ino = posini, posfin
            numno = zi(iacnex-1+ino)
            litrav(numno) = litrav(numno)+1
        end do
2       continue
    end do
!
!     -- ON COMPTE LES NOEUDS COCHES ET ON LES RECOPIE DANS LISTNO:
!
    nbno = 0
    do i = 1, nbnoto
        if (litrav(i) .gt. 0) then
            nbno = nbno+1
            listno(nbno) = i
        end if
    end do
!
    call jedema()
end subroutine
