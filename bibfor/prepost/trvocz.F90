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

subroutine trvocz(ndim, nbvmas, livmas, jrepmo, &
                  jcelds, jcelvs, trxmoy)
    implicit none
    integer(kind=8), intent(in) :: ndim, nbvmas, livmas(nbvmas)
    integer(kind=8), intent(in) :: jrepmo, jcelds, jcelvs
    real(kind=8), intent(out) :: trxmoy
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fgequi.h"
#include "asterfort/r8inir.h"
!
! person_in_charge: sam.cuvilliez at edf.fr
!
! operateur POST_VOIS_CZM : calcul (en 2D et 3D) de la moyenne du taux
! de triaxialite des contraintes dans les elements massifs voisins (1
! ou 2) d'un element cohesif
!
! ----------------------------------------------------------------------
!
! in  nbvmas : nombre (1 ou 2) de mailles voisines dans le massif
! in  livmas : liste des mailles voisines dans le massif
! in  jrepmo : adresse de l'objet '.REPE' du ligrel associe au modele
! in  jcelds : adresse de l'objet '.CELD' du cham_elem SIEF_ELGA
! in  jcelvs : adresse de l'objet '.CELV' du cham_elem SIEF_ELGA
!
! out trxmoy : moyenne du taux de triaxialite
!              -> trxmoy = somme_voisins(somme_pg(triax)/npg)/nvoisins
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: numavo, jgrel, jelem, iainfo, iaval, nlong, nbpg, nbsig
    integer(kind=8) :: idecpg, j, kpg, lgcat
    real(kind=8) :: trxtmp, equi(17)
!
! ----------------------------------------------------------------------
!
! --- initialisations
!
    call r8inir(17, 0.d0, equi, 1)
!
    trxmoy = 0.d0
!
    if (ndim .eq. 2) then
        nbsig = 4
    else if (ndim .eq. 3) then
        nbsig = 6
    else
        ASSERT(.false.)
    end if
!
! -------------------------------------------------------------------
!   boucle sur les mailles voisines
! -------------------------------------------------------------------
!
    do j = 1, nbvmas
!
        trxtmp = 0.d0
!
        numavo = livmas(j)
!
!       recuperation dans l'objet '.REPE' du numero de grel associe
!       a la maille numavo, et de sa position dansce grel
        jgrel = zi(jrepmo-1+2*(numavo-1)+1)
        jelem = zi(jrepmo-1+2*(numavo-1)+2)
!
!       adresse dans l'objet '.CELD' des infos concernant ce grel
        iainfo = zi(jcelds-1+4+jgrel)
!
!       adresse dans l'objet '.CELV' de la premiere valeur du champ
!       de contrainte pour l'element porte par numavo
        iaval = zi(jcelds-1+iainfo+4+4*(jelem-1)+4)
!
!       nombre de points de Gauss
        lgcat = zi(jcelds-1+iainfo+3)
        nlong = zi(jcelds-1+iainfo+4+4*(jelem-1)+3)
        ASSERT(lgcat .eq. nlong)
        nbpg = nlong/nbsig
!
! -------------------------------------------------------------------
!       boucle sur les points de Gauss
! -------------------------------------------------------------------
!
        do kpg = 1, nbpg
!
!           taux de triaxialite
            idecpg = jcelvs-1+iaval+nbsig*(kpg-1)
            call fgequi(zr(idecpg), 'SIGM_DIR', ndim, equi)
            trxtmp = trxtmp+equi(17)
!
        end do
!
!       moyenne arithmetique (sur les PG) du taux de triaxialite
        trxtmp = trxtmp/nbpg
        trxmoy = trxmoy+trxtmp
!
    end do
!
!   moyenne arithmetique (sur EF voisins) du taux de triaxialite
    trxmoy = trxmoy/nbvmas
!
end subroutine
