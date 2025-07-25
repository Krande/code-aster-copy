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
! person_in_charge: mickael.abbas at edf.fr
! aslint: disable=W1504
!
subroutine ndxpre(model, nume_dof, ds_material, cara_elem, &
                  ds_constitutive, list_load, ds_algopara, solveu, ds_system, &
                  list_func_acti, sddisc, ds_measure, nume_inst, hval_incr, &
                  hval_algo, matass, maprec, &
                  sddyna, nlDynaDamping, &
                  sderro, sdnume, hval_meelem, hval_measse, hval_veelem, hval_veasse, &
                  lerrit)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/ndxprm.h"
#include "asterfort/nmassx.h"
#include "asterfort/ndxforc_pred.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmcret.h"
#include "asterfort/nmltev.h"
#include "asterfort/nmresd.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8) :: list_func_acti(*)
    integer(kind=8) :: nume_inst
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19) :: matass, maprec
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19) :: list_load, solveu
    character(len=19), intent(in) :: sdnume, sddisc
    character(len=24) :: model, cara_elem
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    character(len=24) :: nume_dof
    character(len=24) :: sderro
    character(len=19) :: hval_meelem(*), hval_veelem(*)
    character(len=19) :: hval_measse(*), hval_veasse(*)
    character(len=19) :: hval_algo(*), hval_incr(*)
    aster_logical :: lerrit
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! PHASE DE PREDICTION - CAS EXPLICITE
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
! In  ds_system        : datastructure for non-linear system management
! In  ds_inout         : datastructure for input/output management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! In  sdnume           : datastructure for dof positions
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current time step
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  hval_meelem      : hat-variable for elementary matrix
! In  hval_measse      : hat-variable for matrix
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  MAPREC : NOM DE LA MATRICE DE PRECONDITIONNEMENT (GCPC)
! OUT LERRIT  : .TRUE. SI ERREUR PENDANT PREDICTION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: instap
    character(len=19) :: cncine, cndonn, cnzero
    integer(kind=8) :: ldccvg, faccvg, rescvg
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> CALCUL DE PREDICTION'
    end if
!
! --- INITIALISATIONS
!
    instap = diinst(sddisc, nume_inst)
    cndonn = '&&CNCHAR.DONN'
    cnzero = '&&CNPART.ZERO'
    call vtzero(cndonn)
    faccvg = -1
    rescvg = -1
    ldccvg = -1
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(hval_veasse, 'VEASSE', 'CNCINE', cncine)
!
! --- CALCUL DE LA MATRICE GLOBALE
!
    call ndxprm(model, ds_material, cara_elem, ds_constitutive, ds_algopara, &
                list_load, nume_dof, solveu, ds_system, sddisc, &
                sddyna, nlDynaDamping, &
                ds_measure, nume_inst, list_func_acti, &
                hval_incr, hval_algo, hval_meelem, hval_measse, &
                maprec, matass, faccvg, ldccvg)
!
! --- ERREUR SANS POSSIBILITE DE CONTINUER
!
    if ((faccvg .eq. 1) .or. (faccvg .eq. 2) .or. (ldccvg .eq. 1)) then
        goto 99
    end if

! - Compute forces for second member at prediction
    call ndxforc_pred(list_func_acti, &
                      model, cara_elem, list_load, nume_dof, &
                      sddyna, nlDynaDamping, &
                      ds_material, ds_constitutive, ds_system, &
                      ds_measure, sdnume, &
                      sddisc, nume_inst, &
                      hval_incr, hval_algo, &
                      hval_veelem, hval_veasse, &
                      hval_measse, ldccvg)

! - Assemble second member
    call nmassx(list_func_acti, &
                sddyna, nlDynaDamping, &
                hval_veasse, ds_system, &
                cndonn)

! - Solve system
    if (ldccvg .eq. 0) then
        call nmresd(list_func_acti, sddyna, ds_measure, solveu, nume_dof, &
                    instap, maprec, matass, cndonn, cnzero, &
                    cncine, hval_algo, rescvg)
    end if
!
99  continue
!
! --- TRANSFORMATION DES CODES RETOURS EN EVENEMENTS
!
    call nmcret(sderro, 'LDC', ldccvg)
    call nmcret(sderro, 'FAC', faccvg)
    call nmcret(sderro, 'RES', rescvg)
!
! --- EVENEMENT ERREUR ACTIVE ?
!
    call nmltev(sderro, 'ERRI', 'NEWT', lerrit)
!
end subroutine
