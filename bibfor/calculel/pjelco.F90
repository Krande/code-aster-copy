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

subroutine pjelco(moa1, moa2, cham1, corres, base)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
!     COMMANDE:  PROJ_CHAMP /  METHODE='ECLA_PG'
! BUT : CALCULER LA STRUCTURE DE DONNEE CORRESP_2_MAILLA
!       DANS LE CAS OU IL Y A UN CHAM_ELGA A TRAITER (CHAM1)
! ----------------------------------------------------------------------
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/dismoi.h"
#include "asterfort/pjefco.h"
#include "asterfort/pjma1p.h"
#include "asterfort/pjma2p.h"
#include "asterfort/utmess.h"
    character(len=8) :: moa1, moa2
    character(len=16) :: corres
    character(len=19) :: cham1
    character(len=1) :: base
    character(len=8) :: ma1p, ma2p
    integer(kind=8) :: ndim, ndim1, ndim2
!     ----------------------------------------------
!
    ASSERT(base .eq. 'V')
!
!
!     -- CALCUL DE NDIM :
    call dismoi('DIM_GEOM', moa1, 'MODELE', repi=ndim1)
    call dismoi('DIM_GEOM', moa2, 'MODELE', repi=ndim2)
    ASSERT(ndim1 .eq. ndim2)
    ndim = ndim1
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
    call utmess('I', 'CALCULEL3_28', si=ndim)
!
!
!     CREATION DU MAILLAGE 1 PRIME (MA1P)
!     REMPLISSAGE DU .PJEF_MP DANS LA SD CORRES
!     QUI EST LE NOM DU MAILLAGE 1 PRIME
!     ----------------------------------------------
    ma1p = '&&PJELC1'
    call pjma1p(moa1, ma1p, cham1, corres)
    call cargeo(ma1p)
!
!
!     CREATION DU MAILLAGE 2 PRIME (MA2P)
!     REMPLISSAGE DU .PJEF_EL DANS LA SD CORRES
!     QUI EST UN TABLEAU REFERENCANT, POUR CHAQUE ELGA,
!     SON NUMERO ET LE NUMERO DE LA MAILLE A LAQUELLE IL APPARTIENT
!     ----------------------------------------------
    ma2p = '&&PJELC2'
    call pjma2p(ndim, moa2, ma2p, corres)
!
!     -- APPEL A LA ROUTINE "USUELLE" PJEFCO
!        AVEC LES DEUX MAILLAGES PRIME
!     ----------------------------------------------
    call pjefco(ma1p, ma2p, corres, 'V')
!
end subroutine
