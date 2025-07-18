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
subroutine nmpred(modele, numedd, numfix, ds_material, carele, &
                  ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                  fonact, ds_print, ds_measure, ds_algorom, sddisc, &
                  sdnume, sderro, numins, valinc, solalg, &
                  matass, maprec, ds_contact, &
                  sddyna, nlDynaDamping, &
                  meelem, measse, veelem, veasse, lerrit)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/nmcret.h"
#include "asterfort/nmltev.h"
#include "asterfort/nmprde.h"
#include "asterfort/nmprta.h"
#include "asterfort/utmess.h"
#include "asterfort/nonlinDSPrintSepLine.h"
!
    integer(kind=8) :: fonact(*)
    integer(kind=8) :: numins
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19) :: matass, maprec
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19) :: lischa, solveu, sddisc, sdnume
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24) :: modele, carele
    character(len=24) :: numedd, numfix
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(NL_DS_System), intent(in) :: ds_system
    character(len=24) :: sderro
    character(len=19) :: meelem(*), veelem(*)
    character(len=19) :: measse(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    aster_logical :: lerrit
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! PHASE DE PREDICTION
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : LISTE DES CHARGES
! IN  SOLVEU : SOLVEUR
! In  ds_system        : datastructure for non-linear system management
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! In  ds_algopara      : datastructure for algorithm parameters
! IO  ds_print         : datastructure for printing parameters
! IO  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algorom       : datastructure for ROM parameters
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  SDERRO : GESTION DES ERREURS
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IO  ds_print         : datastructure for printing parameters
! IN  SDDYNA : SD DYNAMIQUE
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  MAPREC : NOM DE LA MATRICE DE PRECONDITIONNEMENT (GCPC)
! IN  SDNUME : SD NUMEROTATION
! OUT LERRIT  : .TRUE. SI ERREUR PENDANT PREDICTION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: faccvg, rescvg, ldccvg
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
        call utmess('I', 'MECANONLINE13_33')
    end if
!
! --- INITIALISATION CODES RETOURS
!
    faccvg = -1
    rescvg = -1
    ldccvg = -1
!
! --- PREDICTION PAR LINEARISATION DU SYSTEME
!
    if ((ds_algopara%matrix_pred .eq. 'ELASTIQUE') .or. &
        (ds_algopara%matrix_pred .eq. 'TANGENTE')) then
        call nmprta(modele, numedd, numfix, ds_material, carele, &
                    ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                    fonact, ds_print, ds_measure, ds_algorom, sddisc, &
                    numins, valinc, solalg, matass, maprec, &
                    sddyna, nlDynaDamping, &
                    ds_contact, meelem, measse, veelem, &
                    veasse, sdnume, ldccvg, faccvg, &
                    rescvg)
!
! --- PREDICTION PAR EXTRAPOLATION DU PAS PRECEDENT OU PAR DEPLACEMENT
! --- CALCULE
!
    elseif ((ds_algopara%matrix_pred .eq. 'EXTRAPOLE') .or. &
            (ds_algopara%matrix_pred .eq. 'DEPL_CALCULE')) then
        call nmprde(modele, numedd, numfix, ds_material, carele, &
                    ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                    fonact, ds_print, ds_measure, ds_algorom, sddisc, numins, &
                    valinc, solalg, matass, maprec, ds_contact, &
                    sddyna, nlDynaDamping, &
                    meelem, measse, veelem, veasse, &
                    ldccvg, faccvg, rescvg)
    else
        ASSERT(ASTER_FALSE)
    end if
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
