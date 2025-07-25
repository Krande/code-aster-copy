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
subroutine nmprde(modele, numedd, numfix, ds_material, carele, &
                  ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                  fonact, ds_print, ds_measure, ds_algorom, sddisc, numins, &
                  valinc, solalg, matass, maprec, ds_contact, &
                  sddyna, nlDynaDamping, &
                  meelem, measse, veelem, veasse, &
                  ldccvg, faccvg, rescvg)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmprca.h"
#include "asterfort/nmprdc.h"
#include "asterfort/nmprex.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8) :: fonact(*)
    integer(kind=8) :: numins, ldccvg, faccvg, rescvg
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19) :: maprec, matass
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Print), intent(inout) :: ds_print
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    character(len=19) :: lischa, solveu, sddisc
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24) :: numedd, numfix
    character(len=24) :: modele, carele
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19) :: veelem(*), veasse(*)
    character(len=19) :: meelem(*), measse(*)
    character(len=19) :: solalg(*), valinc(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PREDICTION)
!
! PREDICTION PAR DEPLACEMENT DONNE (EXTRAPOL/DEPL_CALCULE)
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : LISTE DES CHARGES
! IN  MAPREC : MATRICE DE PRECONDITIONNEMENT (GCPC)
! IN  MATASS : MATRICE ASSEMBLEE
! IN  SOLVEU : SOLVEUR
! In  ds_algopara      : datastructure for algorithm parameters
! IO  ds_print         : datastructure for printing parameters
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_measure       : datastructure for measure and statistics management
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IO  ds_contact       : datastructure for contact management
! In  ds_algorom       : datastructure for ROM parameters
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
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
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: incest, depest, depmoi
    character(len=19) :: depso1, depso2
    aster_logical :: lproj
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    depest = '&&CNPART.CHP1'
    incest = '&&CNPART.CHP2'
    call vtzero(depest)
    call vtzero(incest)
    call nmchex(solalg, 'SOLALG', 'DEPSO1', depso1)
    call nmchex(solalg, 'SOLALG', 'DEPSO2', depso2)
    call vtzero(depso1)
    call vtzero(depso2)
    lproj = .true.
    faccvg = -1
    rescvg = -1
    ldccvg = -1
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPMOI', depmoi)
!
! --- VALEUR DU DEPLACEMENT -> DEPEST
! --- VALEUR DE L'INCREMENT DE DEPLACEMENT -> INCEST
!
    if (ds_algopara%matrix_pred .eq. 'EXTRAPOLE') then
        call nmprex(numedd, depmoi, solalg, sddisc, numins, &
                    incest, depest)
    else if (ds_algopara%matrix_pred .eq. 'DEPL_CALCULE') then
        call nmprdc(ds_algopara, numedd, depmoi, sddisc, numins, &
                    incest, depest)
    else
        ASSERT(.false.)
    end if

! --- RECOPIE DE LA SOLUTION
    if (numins .eq. 1) then
        call vtcopy(incest, depso1, iret)
        if (iret .ne. 0) then
            call utmess('F', 'MECANONLINE2_29')
        end if
    else
        call copisd('CHAMP_GD', 'V', incest, depso1)
    end if
!
! --- PROJECTION POUR AVOIR UN CHAMP DE DEPLACEMENT
! --- CINEMATIQUEMENT ADMISSIBLE
!
    if (lproj) then
        call nmprca(modele, numedd, numfix, ds_material, carele, &
                    ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                    fonact, ds_print, ds_measure, ds_algorom, sddisc, numins, &
                    valinc, solalg, matass, maprec, ds_contact, &
                    sddyna, nlDynaDamping, &
                    meelem, measse, veelem, veasse, &
                    depest, ldccvg, faccvg, rescvg)
    end if
!
end subroutine
