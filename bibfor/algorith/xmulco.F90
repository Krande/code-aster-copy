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
!
subroutine xmulco(contac, ddls, ddlc, ddlm, iaint, &
                  ifiss, jheano, vstnc, lact, lcalel, &
                  lelim, ndim, nfh, nfiss, ninter, &
                  nlact, nno, nnol, nnom, nnos, &
                  pla, typma)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/xlacti.h"
#include "asterfort/xplmat.h"
! IN TYPMA : TYPE DE MAILLE
! IN NINTER : NOMBRE DE POINTS D'INTERSECTION
! IN IAINT : ADRESSE TOPOFAC.AI POUR LA FISSURE COURANTE
! OUT LACT : LISTE LAMBDA ACTIFS
! OUT NLACT : NOMBRE LAMBDA ACTIFS
! IN NFISS : NOMBRE DE FISSURES VUES PAR L ELEMENT
! IN JSTNC : ADRESSE TABLEAU FISSURE ACTIVE OU NON
! IN NFH  : NOMBRE ENRICHISSEMENTS HEAVISIDE POUR L ELEMENT
! IN IFISS : INDICE DE LA FISSURE
! IN JHEANO : HEAVISIDE ACTIF AU NOEUD?
! IN CONTAC : P1P1 OU P2P1
! OUT NNOL : NOMBRE DE NOEUDS PORTEURS DE LAGRANGES
! IN NNO : NOMBRE TOT DE NOEUDS
! IN NNOS : NOMBRE NOEUDS SOMMETS
! IN NDIM : DIMENSION
! IN DDLC : NOMBRE DE DDLS DE CONTACT PAR NOEUD
! IN DDLC : NOMBRE DE DDLS PAR NOEUD SOMMET
! IN DDLM : NOMBRE DE DDLS PAR NOEUD MILIEU
! IN NNOM : NOMBRE DE NOEUDS MILIEUX
! IN LCALEL : SI TE DE CALCUL ELEMENTAIRE
! OUT PLA : PLACES DE LAGRANGES DE CONTACT DANS LA MATRICE
! OUT LELIM : YA-T-IL DES LAGRANGES DE CONTACT
!
! --- LISTE DES LAMBDAS ACTIFS
!
    integer :: contac, ddlc, ddlm, iaint, ifiss, jheano, vstnc(*)
    integer :: lact(8), ndim, nfh, nfiss, ninter, nlact
    integer :: nno, nnol, nnom, nnos, pla(27)
    integer :: pli, i, ddls
    aster_logical :: lelim, lcalel
    character(len=8) :: typma
! ----------------------------------------------------------------------
!
    call xlacti(typma, ninter, iaint, lact, nlact)
    if (lcalel) then
        if (nlact .lt. nnos) lelim = .true.
        if (nfiss .eq. 1) then
            do i = 1, nnos
                if (lact(i) .eq. 0) vstnc(i) = 0
            end do
        else
            do i = 1, nnos
                if (lact(i) .eq. 0) vstnc((i-1)*nfh+zi(jheano-1+(i-1)*nfiss+ifiss)) = 0
            end do
        end if
    end if
! --- NOMBRE DE LAMBDAS ET LEUR PLACE DANS LA MATRICE
    if (contac .eq. 1) nnol = nno
    if (contac .eq. 2) nnol = nno
    if (contac .eq. 3) nnol = nnos
    do i = 1, nnol
        call xplmat(ddls, ddlc, ddlm, nnos, nnom, &
                    i, pli)
        if (nfiss .eq. 1) then
            pla(i) = pli
        else
            pla(i) = pli+ndim*(zi(jheano-1+(i-1)*nfiss+ifiss)-1)
        end if
!
    end do
end subroutine
