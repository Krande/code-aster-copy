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

subroutine xmulhm(contac, ddls, ddlc, ddlm, jaint, ifiss, &
                  jheano, vstnc, lact, lcalel, lelim, &
                  nfh, nfiss, ninter, &
                  nlact, nno, nnol, nnom, nnos, &
                  pla, pos, typma, jstano)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
! ======================================================================
! person_in_charge: daniele.colombo at ifpen.fr
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
! IN DDLC : NOMBRE DE DDLS DE CONTACT PAR NOEUD
! IN DDLM :
! IN NNOM : NOMBRE DE NOEUDS MILIEUX
! IN LCALEL : SI TE DE CALCUL ELEMENTAIRE
! OUT PLA : PLACES DE LAGRANGES DE CONTACT DANS LA MATRICE
! OUT LELIM : YA-T-IL DES LAGRANGES DE CONTACT
!
! --- LISTE DES LAMBDAS ACTIFS
!
#include "asterfort/conare.h"
#include "asterfort/xplmat.h"
#include "asterfort/xxmmvd.h"
    integer :: contac, ddlc, ddlm, jaint, ifiss, jheano, vstnc(*)
    integer :: lact(16), nfh, nfiss, ninter, nlact(2)
    integer :: nno, nnol, nnom, nnos, pla(27), vit(16), nli
    integer :: pli, i, ddls, ino, iar, nvit, ino1, ino2
    integer :: ar(12, 3), nbar, zxain, pos(16)
    integer, optional, intent(in) :: jstano
    aster_logical :: lelim, lcalel
    character(len=8) :: typma
! ----------------------------------------------------------------------
!
    do ino = 1, 16
        lact(ino) = 0
        vit(ino) = 0
    end do
    nlact(1:2) = 0
    call conare(typma, ar, nbar)
    zxain = xxmmvd('ZXAIN')
!
! --- ON ACTIVE LES NOEUDS CONNECTES AUX POINTS D'INTERSECTION
!
    do nli = 1, ninter
        iar = int(zr(jaint-1+zxain*(nli-1)+1))
        ino = int(zr(jaint-1+zxain*(nli-1)+2))
        nvit = int(zr(jaint-1+zxain*(nli-1)+5))
        if (ino .gt. 0) then
            lact(ino) = nli
        else if (iar .gt. 0) then
            ino1 = ar(iar, 1)
            ino2 = ar(iar, 2)
            if (nvit .eq. 1) then
                lact(ino1) = nli
                vit(ino1) = 1
                lact(ino2) = nli
                vit(ino2) = 1
            else
                if (vit(ino1) .eq. 0) lact(ino1) = nli
                if (vit(ino2) .eq. 0) lact(ino2) = nli
            end if
        end if
    end do
! --- ON COMPTE LE NOMBRE DE NOEUDS ACTIFS
    do ino = 1, 8
        if (lact(ino) .ne. 0) nlact(1) = nlact(1)+1
    end do
!
    if (lcalel) then
        if (nlact(1) .lt. nnos .or. nlact(2) .lt. nnos) lelim = .true.
        if (nfiss .eq. 1) then
            do i = 1, nnos
                if (lact(i) .eq. 0) vstnc(i) = 0
            end do
        else
            if (ninter .gt. 0) then
                do i = 1, nnos
                    if (lact(i) .eq. 0) then
                        vstnc((i-1)*nfh+zi(jheano-1+(i-1)*nfiss+ifiss)) = 0
                    end if
                end do
            else
                do i = 1, nnos
                    if (lact(i) .eq. 0 .and. zi(jstano-1+(i-1)*nfiss+ifiss) .eq. 0) then
                        vstnc((i-1)*nfh+zi(jheano-1+(i-1)*nfiss+ifiss)) = 0
                    end if
                end do
            end if
!
        end if
    end if
! --- NOMBRE DE LAMBDAS ET LEUR PLACE DANS LA MATRICE
    if (contac .eq. 1) nnol = nno
    if (contac .eq. 2) nnol = nnos
    if (contac .eq. 3) nnol = nnos
    do i = 1, nnol
        call xplmat(ddls, ddlc, ddlm, &
                    nnos, nnom, i, pli)
        if (nfiss .eq. 1) then
            pla(i) = pli
            pos(i) = 1
        else
            pla(i) = pli+ddlc/nfh*(zi(jheano-1+(i-1)*nfiss+ifiss)-1)
            pos(i) = zi(jheano-1+(i-1)*nfiss+ifiss)
        end if
!
    end do
end subroutine
