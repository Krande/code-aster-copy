! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
#include "asterfort/ccbcop.h"
#include "asterfort/cgComputeGtheta.h"
#include "asterfort/cgComputeTheta.h"
#include "asterfort/cgExportTableG.h"
#include "asterfort/cgTableG.h"
#include "asterfort/cgVerification.h"
#include "asterfort/deprecated_algom.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
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
!
    integer           :: iopt, istore, jopt, nbropt
    character(len=19) :: lisopt
!---------------------------------------------------------------------------------------------------
    call jemarq()
    call infmaj()
    call deprecated_algom('CALC_H')
!
! Fiches concernées par le chantier (A supprimer à la fin)
! A Faire: #29573, #27931, #29703, #30288
!
!-- Initialisation des champs et des paramètres
    call cgField%initialize()
    call cgTheta%initialize()
    call cgStudy%initialize(cgField%result_in, cgField%list_nume(1))
!
!-- Calcul de la courbure
    if (cgField%ndim == 3) then
        call cgTheta%compute_courbature(cgStudy%model)
    endif
!
!-- Verification (A nettoyer)
    call cgVerification(cgField, cgTheta)
!
!-- Compute Theta factors
    call cgComputeTheta(cgField, cgTheta)
!
!-- ELAS INCR
    if (cgField%l_incr) then

        lisopt = '&&OP0027.LISOPT'
        nbropt = 2
!
        call wkvect(lisopt, 'V V K16', nbropt, jopt)
        zk16(jopt) = 'VARI_ELNO'
        zk16(jopt+1) = 'EPSP_ELNO'
!
        call ccbcop(cgField%result_in, cgField%result_out, cgField%list_nume_name,&
                    cgField%nb_nume, lisopt, nbropt)
    endif
!
!-- Loop on option
    do iopt = 1, cgField%nb_option
!
        call cgStudy%setOption(cgField%list_option(iopt))
!
        if (cgField%isModeMeca()) then
            if (cgStudy%option .eq. 'K') then
                cgStudy%l_modal = ASTER_TRUE
            else
                call utmess('F', 'RUPTURE0_27')
            endif
        endif
!
        do istore = 1, cgField%nb_nume
!
            call cgStudy%initialize(cgField%result_in, cgField%list_nume(istore))
!
! --------  Maillage similaire sd_fond_fissure et sd_resu
            ASSERT(cgTheta%mesh == cgStudy%mesh)
!
! --------  Récupération des champs utiles pour l'appel à calcul
            call cgStudy%getField(cgField%result_in)
!
            call cgStudy%getParameter(cgField%result_in)
!
!---------- Calcul de G(theta) pour les éléments 2D/3D option G et K
            call cgComputeGtheta(cgField, cgTheta, cgStudy)
!
!---------- Création de la table de G et des SIFS
            call cgTableG(cgField, cgTheta, cgStudy)
!
        end do
!
    end do
!
!------ Print fields
    call cgTheta%print()
    call cgField%print()
!
! --- Cleaning
    call cgField%clean()
!
!-- Création de la table container
    call cgExportTableG(cgField, cgTheta)
!
    call jedema()
!
end subroutine
