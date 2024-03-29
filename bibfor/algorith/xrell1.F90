! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine xrell1(tabl_node, nb_edge, node_sele, nb_node_sele, sdline_crack, &
                  tabai, l_ainter)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: samuel.geniaut at edf.fr
! aslint: disable=W1306
!
    integer, intent(in) :: nb_edge
    integer, intent(in) :: nb_node_sele
    integer, intent(in) :: tabl_node(3, nb_edge)
    integer, intent(inout) :: node_sele(nb_edge)
    character(len=14), intent(in) :: sdline_crack
    aster_logical, intent(in) :: l_ainter
    character(len=19) :: tabai
!
! --------------------------------------------------------------------------------------------------
!
! XFEM - Contact definition
!
! Create list of linear relations (algorithm 1)
!
! --------------------------------------------------------------------------------------------------
!
! (VOIR BOOK VI 15/07/05)
!    - CREATION DES RELATIONS DE LIAISONS ENTRE LAGRANGE
!
! --------------------------------------------------------------------------------------------------
!
! In  sdline_crack   : name of datastructure of linear relations for crack
! In  nb_edge        : number of cut edges
! In  tabl_node      : table of nodes for edges (middle et vertex nodes)
! In  node_sele      : selected nodes
! In  nb_node_sele   : number of selected nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, in, dimeq, ia, ext, libre, k, eq(100), tabno(nb_edge, 3), ie
    integer :: liseqt(nb_edge, 2), nreleq, jlis1, jtabai, ier
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    do i = 1, nb_edge
        do j = 1, 3
            tabno(i, j) = tabl_node(j, i)
        end do
    end do
    nreleq = 0
!
! --- CREATION DU TABLEAU TEMPORAIRE DES RELATION D'EGALITE : LISEQT
!
    do i = 1, nb_node_sele
        in = node_sele(i)
        dimeq = 0
        do ia = 1, nb_edge
            do j = 1, 2
!           ON CHERCHE LES ARETES EMANANTES
                if (tabno(ia, j) .eq. in) then
                    ext = tabno(ia, 3-j)
!             ON REGARDE SI L'AUTRE EXTREMITE EST LIBRE
                    libre = 1
                    do k = 1, nb_node_sele
                        if (ext .eq. node_sele(k)) libre = 0
                    end do
                    if (libre .eq. 1) then
                        dimeq = dimeq+1
                        eq(dimeq) = ia
                    end if
                end if
            end do
        end do
        do ie = 1, dimeq
            nreleq = nreleq+1
            liseqt(nreleq, 1) = tabno(eq(ie), 1)
            liseqt(nreleq, 2) = tabno(eq(ie), 2)
        end do
    end do
!
! --- Stockage du tableau temporaire: liaison dans tabai (5ème composante de TOPOFAC.AI)
!
    if (l_ainter) then
        call jeexin(tabai, ier)
        ASSERT(ier .eq. 0)
        if (nreleq .ne. 0) then
            call wkvect(tabai, 'V V I', nreleq*3, jtabai)
            do i = 1, nreleq
                zi(jtabai-1+3*(i-1)+1) = liseqt(i, 1)
                zi(jtabai-1+3*(i-1)+2) = liseqt(i, 2)
                zi(jtabai-1+3*(i-1)+3) = 1
            end do
        end if
    end if
!
! --- STOCKAGE DE LISEQT
!
    if (nreleq .gt. 0) then
        call wkvect(sdline_crack, 'G V I', nreleq*2, jlis1)
        do ie = 1, nreleq
            zi(jlis1-1+2*(ie-1)+1) = liseqt(ie, 1)
            zi(jlis1-1+2*(ie-1)+2) = liseqt(ie, 2)
        end do
    end if
!
    call jedema()
end subroutine
