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

subroutine te0496(option, nomte)
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL POUR TOUS LES ÉLÉMENTS DU NOMBRE
!       - DE VARIABLES INTERNES
!       - DE SOUS-POINTS
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/lteatt.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!
    integer(kind=8) :: nbLayer, npgh, nbSect, nbFibr, nbvari, jcompo, jdcel, jnbsp
    integer(kind=8) :: itab(2), iret
!
! --------------------------------------------------------------------------------------------------
!
!     PCOMPOR contient 1 seule cmp : NBVARI
!     Si PCOMPOR n'est pas fourni  : NBVARI = 0
!     Les autres infos viennent de PNBSP_I
! --------------------------------------------------------------------------------------------------
    call jevech('PDCEL_I', 'E', jdcel)
    call tecach('ONO', 'PCOMPOR', 'L', iret, nval=2, itab=itab)
    if (itab(1) .ne. 0) then
        ASSERT(itab(2) .eq. 1)
        jcompo = itab(1)
        read (zk16(jcompo-1+1), '(I16)') nbvari
    else
        nbvari = 0
    end if
!
!   Par défaut pour tous les éléments
    zi(jdcel-1+1) = 1
    zi(jdcel-1+2) = nbvari
!
!   PNBSP_I : infos nécessaires au calcul du nombre de sous-points.
!   Si le champ n'est pas donné ==> valeur par défaut
    call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jnbsp)
    if (jnbsp .eq. 0) goto 999
!
!   Cas des éléments "coque épaisse" (multi-couche)
    if ((nomte .eq. 'MEC3QU9H') .or. (nomte .eq. 'MEC3TR7H') .or. (nomte .eq. 'MECXSE3')) then
        nbLayer = zi(jnbsp-1+1)
        npgh = 3
        zi(jdcel-1+1) = npgh*nbLayer
!
!   Cas des éléments "DKT"
    else if ((nomte .eq. 'MEDKQU4') .or. (nomte .eq. 'MEDKTR3') .or. &
             (nomte .eq. 'MEDSQU4') .or. (nomte .eq. 'MEDSTR3') .or. &
             (nomte .eq. 'MEQ4QU4') .or. (nomte .eq. 'MET3TR3')) then
        nbLayer = zi(jnbsp-1+1)
        npgh = 3
        zi(jdcel-1+1) = npgh*nbLayer

!   Cas des éléments "TUYAU"
    else if (lteatt('TUYAU', 'OUI')) then
        nbLayer = zi(jnbsp-1+1)
        nbSect = zi(jnbsp-1+2)
        zi(jdcel-1+1) = (2*nbSect+1)*(2*nbLayer+1)
!
!   Cas des éléments de poutre "multi-fibres"
    else if (nomte .eq. 'MECA_POU_D_EM') then
        nbFibr = zi(jnbsp-1+1)
        zi(jdcel-1+1) = nbFibr
!
    else if (nomte .eq. 'MECA_POU_D_TGM') then
        nbFibr = zi(jnbsp-1+1)
        zi(jdcel-1+1) = nbFibr
!
    else if (nomte .eq. 'MECA_POU_D_SQUE') then
        nbFibr = zi(jnbsp-1+1)
        zi(jdcel-1+1) = nbFibr
    end if
!
999 continue
end subroutine
