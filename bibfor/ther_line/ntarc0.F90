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
subroutine ntarc0(result, model, materfield, caraElem, listLoadResu, &
                  para, nume_store, time_curr)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rssepa.h"
!
    character(len=8), intent(in) :: result, model, materfield, caraElem
    integer(kind=8), intent(in) :: nume_store
    real(kind=8), intent(in) :: time_curr
    real(kind=8), intent(in) :: para(*)
    character(len=24), intent(in) :: listLoadResu
!
! --------------------------------------------------------------------------------------------------
!
! THER_* - Storing
!
! Storing parameters
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESULT : NOM UTILISATEUR DU CONCEPT RESULTAT
! IN  MODELE : NOM DU MODELE
! IN  MATE   : CHAMP DE MATERIAU
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! IN  PARA   : PARAMETRES DU CALCUL
!               (1) THETA
!               (2) DELTAT
! IN  SDCRIT : VALEUR DES CRITERES DE CONVERGENCE
! IN  LISCH2 : NOM DE LA SD INFO CHARGE POUR STOCKAGE DANS LA SD
! IN  NUMARC : NUMERO D'ARCHIVAGE
! IN  INSTAN : VALEUR DE L'INSTANT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jvPara
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Store time
!
    call rsadpa(result, 'E', 1, 'INST', nume_store, 0, sjv=jvPara)
    zr(jvPara) = time_curr
!
! - Store some parameters
!
    call rsadpa(result, 'E', 1, 'PARM_THETA', nume_store, 0, sjv=jvPara)
    zr(jvPara) = para(1)
!
! - Store others
!
    call rssepa(result, nume_store, model, materfield, caraElem, listLoadResu)
!
    call jedema()
end subroutine
