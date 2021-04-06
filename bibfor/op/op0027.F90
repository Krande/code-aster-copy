! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: tanguy.mathieu at edf.fr
!
subroutine op0027()
!
use calcG_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cgComputeGtheta.h"
#include "asterfort/cgComputeMatrix.h"
#include "asterfort/cgComputeTheta.h"
#include "asterfort/cgExportTableG.h"
#include "asterfort/cgVerification.h"
#include "asterfort/deprecated_algom.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "jeveux.h"
! --------------------------------------------------------------------------------------------------
!
!      OPERATEUR :     CALC_H
!
!      BUT:CALCUL DU TAUX DE RESTITUTION D'ENERGIE PAR LA METHODE THETA
!          CALCUL DES FACTEURS D'INTENSITE DE CONTRAINTES
!
!---------------------------------------------------------------------------------------------------
!
    type(CalcG_field) :: cgField
    type(CalcG_theta) :: cgTheta
    type(CalcG_study) :: cgStudy
    type(CalcG_table) :: cgTable
!
    integer           :: i_opt, i_nume
!---------------------------------------------------------------------------------------------------
    call jemarq()
    call infmaj()
    call deprecated_algom('CALC_H')
!
! Fiches concernées par le chantier (A supprimer à la fin)
! A Faire: , #27931, #30288
!
!-- Initialisation des champs et des paramètres
    call cgField%initialize()
    call cgTheta%initialize()
    call cgTable%initialize(cgField, cgTheta)
    call cgStudy%initialize(cgField%result_in, cgField%list_nume(1))
!
!-- Calcul de la courbure
    if (cgField%ndim == 3) then
        call cgTheta%compute_curvature(cgStudy%model)
    endif
!
!-- Verification (A nettoyer)
    call cgVerification(cgField, cgTheta, cgStudy)
!
!-- Compute Theta factors
    call cgComputeTheta(cgField, cgTheta)
!
! --- Compute A Matrix from equation A*G(s)=g(theta)
!
    call cgComputeMatrix(cgField, cgTheta)
!
!-- Loop on option
    do i_nume = 1, cgField%nb_nume
!
        call cgStudy%initialize(cgField%result_in, cgField%list_nume(i_nume))
!
! ----  Récupération des champs utiles pour l'appel à calcul
        call cgStudy%getField(cgField%result_in)
!
! ----  Maillage similaire sd_fond_fissure et sd_resu
        ASSERT(cgTheta%mesh == cgStudy%mesh)
!
        do i_opt = 1, cgField%nb_option
!
            call cgStudy%setOption(cgField%list_option(i_opt), cgField%isModeMeca())
!
            call cgStudy%getParameter(cgField%result_in)
!
!---------- Calcul de G(theta) pour les éléments 2D/3D option G et K
            call cgComputeGtheta(cgField, cgTheta, cgStudy, cgTable)
        end do
!
        call cgTable%save(cgField, cgTheta, cgStudy)
!
    end do
!
! --- Cleaning
    call cgField%clean()
!
!-- Création de la table container
    call cgExportTableG(cgField, cgTheta, cgTable)
!
    call jedema()
!
end subroutine
