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
subroutine avsieq(nbordr, vwork, tdisp, kwork, sommw, &
                  tspaq, i, jvsieq)
!
    implicit none
#include "jeveux.h"
#include "asterfort/fgequi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    integer(kind=8) :: nbordr, tdisp, kwork, sommw, tspaq, i, jvsieq
    real(kind=8) :: vwork(tdisp)
!, vsieq(nbordr)
! ----------------------------------------------------------------------
! BUT: CALCULER LA CONTRAINTE EQUIVALENTE
! ----------------------------------------------------------------------
! ARGUMENTS :
!  NBORDR : IN   I  : NOMBRE DE NUMEROS D'ORDRE.
!  VWORK  : IN   R  : VECTEUR DE TRAVAIL CONTENANT
!                     L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                     ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                     DU <<PAQUET>> DE MAILLES.
!  TDISP  : IN   I  : TAILLE DU VECTEUR DE TRAVAIL.
!  KWORK  : IN   I  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET
!                               MAILLES OU LE 1ER NOEUD DU PAQUET DE
!                               NOEUDS;
!                     KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                               MAILLES OU LE IEME NOEUD DU PAQUET
!                               DE NOEUDS.
!  SOMMW  : IN   I  : SOMME DES POINTS DE GAUSS OU DES NOEUDS DES N
!                     MAILLES PRECEDANT LA MAILLE COURANTE.
!  TSPAQ  : IN   I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                     OU DE NOEUDS COURANT.
!  I      : IN   I  : IEME POINT DE GAUSS OU IEME NOEUD.
!  VSIEQ : OUT  I  : ADDRESS VECTEUR CONTENANT  A
!                     TOUS LES INSTANTS.
! ----------------------------------------------------------------------
    integer(kind=8) :: iordr, adrs, j, decal
    real(kind=8) :: sig(6), equi(17)
!     ------------------------------------------------------------------
!
!234567                                                              012
!
    call jemarq()
!
    decal = 18
    do iordr = 1, nbordr
        adrs = (iordr-1)*tspaq+kwork*sommw*decal+(i-1)*decal
!
        do j = 1, 6
            sig(j) = vwork(adrs+j)
        end do
!
        call fgequi(sig, 'SIGM', 3, equi)
        zr(jvsieq+iordr) = equi(1)
!
    end do
!
    call jedema()
!
end subroutine
