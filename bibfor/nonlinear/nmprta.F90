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
! aslint: disable=W1504
!
subroutine nmprta(model, nume_dof, numfix, ds_material, cara_elem, &
                  ds_constitutive, list_load, ds_algopara, solveu, ds_system, &
                  list_func_acti, ds_print, ds_measure, ds_algorom, sddisc, &
                  nume_inst, hval_incr, hval_algo, matass, maprec, &
                  sddyna, nlDynaDamping, &
                  ds_contact, hval_meelem, hval_measse, hval_veelem, &
                  hval_veasse, sdnume, ldccvg, faccvg, &
                  rescvg)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmassp.h"
#include "asterfort/nmforc_pred.h"
#include "asterfort/nonlinIntForce.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdep0.h"
#include "asterfort/nmfocc.h"
#include "asterfort/nmprma.h"
#include "asterfort/nmresd.h"
#include "asterfort/vtzero.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: list_func_acti(*)
    integer(kind=8) :: nume_inst, faccvg, rescvg, ldccvg
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19) :: matass, maprec
    character(len=19) :: list_load, solveu, sddisc, sdnume
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24) :: model, cara_elem
    character(len=24) :: nume_dof, numfix
    character(len=19) :: hval_algo(*), hval_incr(*)
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19) :: hval_veelem(*), hval_veasse(*)
    character(len=19) :: hval_meelem(*), hval_measse(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PREDICTION)
!
! PREDICTION PAR METHODE DE NEWTON-EULER
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  list_load        : name of datastructure for list of loads
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! In  ds_system        : datastructure for non-linear system management
! IO  ds_print         : datastructure for printing parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algorom       : datastructure for ROM parameters
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current time step
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_contact       : datastructure for contact management
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  MAPREC : NOM DE LA MATRICE DE PRECONDITIONNEMENT (GCPC)
! IN  SDNUME : SD NUMEROTATION
! OUT FACCVG : CODE RETOUR FACTORISATION MATRICE GLOBALE
!                -1 : PAS DE FACTORISATION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : MATRICE SINGULIERE
!                 2 : ERREUR LORS DE LA FACTORISATION
!                 3 : ON NE SAIT PAS SI SINGULIERE
! OUT RESCVG : CODE RETOUR RESOLUTION SYSTEME LINEAIRE
!                -1 : PAS DE RESOLUTION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : NOMBRE MAXIMUM D'ITERATIONS ATTEINT
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 2 : ERREUR SUR LA NON VERIF. DE CRITERES PHYSIQUES
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
! OUT CONDCVG : CODE RETOUR DE LA CONDANSATION STATIQUE
!                -1 : PAS DE CONDENSATION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE LA CONDENSATION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = PRED_EULER
    integer(kind=8), parameter :: iterNewtPred = 0
    integer(kind=8) :: ifm, niv
    real(kind=8) :: time_curr
    character(len=19) :: cncine, cndonn, cnpilo
    aster_logical :: leltc
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_34')
    end if
!
! - Active functionnalities
!
    leltc = isfonc(list_func_acti, 'ELT_CONTACT')
!
! --- INITIALISATIONS
!
    ldccvg = -1
    faccvg = -1
    rescvg = -1
    cndonn = '&&CNCHAR.DONN'
    cnpilo = '&&CNCHAR.PILO'
    call vtzero(cndonn)
    call vtzero(cnpilo)

! - Get time
    ASSERT(nume_inst .gt. 0)
    time_curr = diinst(sddisc, nume_inst)

!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(hval_veasse, 'VEASSE', 'CNCINE', cncine)
!
! --- INCREMENT DE DEPLACEMENT NUL EN PREDICTION
!
    call nmdep0('ON ', hval_algo)

! - Compute matrix
    call nmprma(list_func_acti, &
                model, cara_elem, &
                ds_material, ds_constitutive, &
                list_load, sddyna, nlDynaDamping, &
                sddisc, nume_inst, &
                ds_algopara, ds_contact, ds_algorom, &
                ds_print, ds_measure, &
                hval_incr, hval_algo, &
                hval_meelem, hval_measse, &
                nume_dof, numfix, &
                solveu, ds_system, &
                maprec, matass, &
                faccvg, ldccvg)

!
! --- ERREUR SANS POSSIBILITE DE CONTINUER
!
    if ((faccvg .eq. 1) .or. (faccvg .eq. 2) .or. (ldccvg .eq. 1)) then
        goto 999
    end if
!
! - Compute forces for second member at prediction
!
    call nmforc_pred(list_func_acti, &
                     model, cara_elem, list_load, &
                     nume_dof, matass, &
                     sddyna, nlDynaDamping, &
                     ds_material, ds_constitutive, &
                     ds_measure, ds_algopara, &
                     sddisc, nume_inst, &
                     hval_incr, hval_algo, &
                     hval_veelem, hval_veasse, &
                     hval_measse)
!

! --- CALCUL DU SECOND MEMBRE POUR CONTACT
!
    if (leltc) then
        call nmfocc('PREDICTION', model, ds_material, nume_dof, list_func_acti, &
                    ds_contact, ds_measure, hval_algo, hval_incr, ds_constitutive)
    end if
!
! --- INCREMENT DE DEPLACEMENT NUL EN PREDICTION
!
    call nmdep0('OFF', hval_algo)
!
! - Compute internal forces
!
    call nonlinIntForce(phaseType, &
                        model, cara_elem, &
                        list_func_acti, iterNewtPred, sdnume, &
                        ds_material, ds_constitutive, &
                        ds_system, ds_measure, &
                        hval_incr, hval_algo, &
                        ldccvg, &
                        sddyna_=sddyna, &
                        ds_algorom_=ds_algorom)
    if (ldccvg .eq. 1) then
        goto 999
    end if

! - Evaluate second member for prediction
    call nmassp(list_func_acti, &
                sddyna, nlDynaDamping, &
                ds_system, ds_contact, hval_veasse, &
                cnpilo, cndonn)
!
! --- RESOLUTION K.DU = DF
!
    call nmresd(list_func_acti, sddyna, ds_measure, solveu, nume_dof, &
                time_curr, maprec, matass, cndonn, cnpilo, &
                cncine, hval_algo, rescvg, ds_algorom)
!
999 continue
!
end subroutine
