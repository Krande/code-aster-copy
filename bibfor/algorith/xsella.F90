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

subroutine xsella(crack, nb_node_mesh, nb_edge, tabl_node, node_sele, &
                  nb_node_sele)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: samuel.geniaut at edf.fr
! aslint: disable=W1306
!
    character(len=8), intent(in) :: crack
    integer, intent(in) :: nb_node_mesh
    integer, intent(in) :: nb_edge
    integer, intent(in) :: tabl_node(3, nb_edge)
    integer, intent(inout) :: node_sele(nb_edge)
    integer, intent(out) :: nb_node_sele
!
! --------------------------------------------------------------------------------------------------
!
! XFEM - Contact definition
!
! Selection of edges
!
! --------------------------------------------------------------------------------------------------
!
! (VOIR BOOK VI 15/07/05) :
!    - SELECTION DES NOEUDS POUR LES LAGRANGES
!
! --------------------------------------------------------------------------------------------------
!
! In  crack          : name of crack
! In  nb_node_mesh   : number of (physical) nodes in mesh
! In  nb_edge        : number of cut edges
! In  tabl_node      : table of nodes for edges (middle et vertex nodes)
! IO  node_sele      : selected nodes
! Out nb_node_sele   : number of selected nodes
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j, k, ik, cptno, deja, noeud(2*nb_edge), tabdir(nb_edge, 2)
    integer :: scorno(2*nb_edge), scorar(nb_edge), ia, max, bestar, no1, no2
    integer :: bestno, noconn, narcas, liarca(nb_edge), nar, ni, nj
    integer :: tabno(nb_edge, 3), ii, cpt
    real(kind=8) :: scorn2(2*nb_edge), scora2(nb_edge), li, lj, maxr
    character(len=19) :: cnsln
    real(kind=8), pointer :: cnsv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
! --- TRANSFO. CHAM_NO -> CHAM_NO_S DE LA LEVEL SET NORMALE
!
    cnsln = '&&XSELLA.CNSLN'
    call cnocns(crack(1:8)//'.LNNO', 'V', cnsln)
    call jeveuo(cnsln(1:19)//'.CNSV', 'L', vr=cnsv)
!
! --- INITIALISATIONS
!
    nar = nb_edge
    do i = 1, nar
        do j = 1, 3
            tabno(i, j) = tabl_node(j, i)
        end do
    end do
    cptno = 0
    scorno(1:2*nar) = 0
    noeud(1:2*nar) = 0
    cpt = 0
!
! --- CALCUL DE SCORNO : NB D'ARETES CONNECTESS AU NOEUD
! --- CALCUL DE SCORN2 : VALEUR ABSOLUE DE LA LEVEL SET NORMALE
!
    do i = 1, nar
        do j = 1, 2
            deja = 0
            do k = 1, cpt
                if (tabno(i, j) .eq. noeud(k)) then
                    deja = 1
                    ik = k
                end if
            end do
            if (deja .eq. 0) then
                cpt = cpt+1
                noeud(cpt) = tabno(i, j)
                tabdir(i, j) = cpt
                scorn2(cpt) = abs(cnsv(tabno(i, j)))
            else
                tabdir(i, j) = ik
            end if
            scorno(tabdir(i, j)) = scorno(tabdir(i, j))+1
        end do
    end do
!
! --- BOUCLE TANT QUE TOUTES LES ARETES NE SONT PAS CASSEES
!
    do ii = 1, nb_node_mesh
!
        if (nar .eq. 0) goto 300
!
! --- CALCUL SCORAR : DIFF DE SCORE POUR CHAQUE NOEUD
! --- CALCUL SCOAR2 : RAPPORT DES LEVEL SETS
!
        do ia = 1, nar
            ni = scorno(tabdir(ia, 1))
            nj = scorno(tabdir(ia, 2))
            scorar(ia) = abs(ni-nj)
            li = abs(cnsv(tabno(ia, 1)))
            lj = abs(cnsv(tabno(ia, 2)))
            if (ni .gt. nj) scora2(ia) = li/(li+lj)
            if (ni .lt. nj) scora2(ia) = lj/(li+lj)
            if (ni .eq. nj) scora2(ia) = min(li, lj)/(li+lj)
        end do
!
! --- MEILLEURE ARETE : PICK MEILLEUR NOEUD
!
        max = -1
        maxr = -1.d0
        do ia = 1, nar
            if ((scorar(ia) .gt. max) .or. (scorar(ia) .eq. max .and. scora2(ia) .ge. maxr)) then
                max = scorar(ia)
                maxr = scora2(ia)
                bestar = ia
            end if
        end do
!
        no1 = tabdir(bestar, 1)
        no2 = tabdir(bestar, 2)
        bestno = no1
        if ((scorno(no2) .gt. scorno(no1)) .or. &
            (scorno(no2) .eq. scorno(no1) .and. scorn2(no2) .lt. scorn2(no1))) then
            bestno = no2
        end if
        cptno = cptno+1
        node_sele(cptno) = noeud(bestno)
!
! --- UPDATE SCORE DES NOEUDS
!
        do i = 1, nar
            do j = 1, 2
                if (tabdir(i, j) .eq. bestno) then
                    scorno(bestno) = scorno(bestno)-1
                    noconn = tabdir(i, 3-j)
                    scorno(noconn) = scorno(noconn)-1
                end if
            end do
        end do
!
! --- ON LISTE LES ARETES CONNECTEES A�BESTNO PUIS ON LES CASSE
!
        narcas = 0
        do ia = 1, nar
            if (tabdir(ia, 1) .eq. bestno .or. tabdir(ia, 2) .eq. bestno) then
                narcas = narcas+1
                liarca(narcas) = ia
            end if
        end do
!
! --- FEINTE : ON LES SUPPRIME EN PARTANT PAR LA FIN !
!
        do ia = narcas, 1, -1
!          TABNO(LIARCA(IA),:)=[]
!          TABDIR(LIARCA(IA),:)=[]
            do i = liarca(ia), nar-1
                tabno(i, 1) = tabno(i+1, 1)
                tabno(i, 2) = tabno(i+1, 2)
                tabno(i, 3) = tabno(i+1, 3)
                tabdir(i, 1) = tabdir(i+1, 1)
                tabdir(i, 2) = tabdir(i+1, 2)
            end do
            tabno(nar, 1) = 0
            tabno(nar, 2) = 0
            tabno(nar, 3) = 0
            tabdir(nar, 1) = 0
            tabdir(nar, 2) = 0
            nar = nar-1
        end do
    end do
300 continue
!
! --- NOMBRE DE PICKED NODES
!
    nb_node_sele = cptno
!
! --- DESTRUCTION DES OBJETS TEMPORAIRES
!
    call detrsd('CHAM_NO_S', cnsln)
!
end subroutine
