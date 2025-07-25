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
subroutine ndexpl(modele, numedd, ds_material, carele, &
                  ds_constitutive, lischa, ds_algopara, fonact, ds_system, &
                  ds_print, ds_measure, sdnume, &
                  sddyna, nlDynaDamping, &
                  sddisc, sderro, valinc, numins, solalg, solveu, &
                  matass, maprec, ds_inout, meelem, measse, &
                  veelem, veasse, nbiter)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/ndxcvg.h"
#include "asterfort/ndxdec.h"
#include "asterfort/ndxdep.h"
#include "asterfort/ndxnpa.h"
#include "asterfort/ndxpre.h"
#include "asterfort/ndxforc_step.h"
!
    integer(kind=8) :: numins
    integer(kind=8) :: fonact(*)
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=24) :: sderro
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: sdnume, sddisc
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_InOut), intent(in) :: ds_inout
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=19) :: valinc(*), solalg(*)
    character(len=19) :: meelem(*), veelem(*)
    character(len=19) :: measse(*), veasse(*)
    character(len=19) :: lischa
    character(len=19) :: solveu, maprec, matass
    character(len=24) :: modele, numedd
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24) :: carele
    integer(kind=8) :: nbiter
!
! --------------------------------------------------------------------------------------------------
!
! OPERATEUR NON-LINEAIRE MECANIQUE
!
! ALGORITHME DYNAMIQUE EXPLICITE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : L_CHARGES
! In  ds_inout         : datastructure for input/output management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  MAPREC : NOM DE LA MATRICE DE PRECONDITIONNEMENT (GCPC)
! OUT NBITER : NOMBRE D'ITERATIONS DE NEWTON
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lerrit
!
! --------------------------------------------------------------------------------------------------
!

! - Update for new time step
    call ndxnpa(modele, carele, &
                fonact, ds_print, &
                ds_material, ds_constitutive, &
                sddisc, sddyna, sdnume, numedd, numins, &
                valinc, solalg)

! - Compute forces for second member when constant in time step
    call ndxforc_step(fonact, &
                      modele, carele, numedd, &
                      lischa, sddyna, &
                      ds_material, ds_constitutive, &
                      ds_measure, ds_inout, &
                      sddisc, numins, &
                      valinc, solalg, &
                      veelem, veasse)
!
! --- PREDICTION D'UNE DIRECTION DE DESCENTE
!
    call ndxpre(modele, numedd, ds_material, carele, &
                ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                fonact, sddisc, ds_measure, numins, valinc, &
                solalg, matass, maprec, &
                sddyna, nlDynaDamping, &
                sderro, sdnume, meelem, measse, veelem, veasse, &
                lerrit)
!
! - Update displacements
!
    if (.not. lerrit) then
        call ndxdep(numedd, fonact, numins, sddisc, &
                    sddyna, nlDynaDamping, &
                    sdnume, valinc, solalg, veasse)
    end if
!
! --- ESTIMATION DE LA CONVERGENCE
!
    call ndxcvg(sddisc, sderro, valinc)
!
! --- EN L'ABSENCE DE CONVERGENCE ON CHERCHE A SUBDIVISER LE PAS
! --- DE TEMPS SI L'UTILISATEUR A FAIT LA DEMANDE
!
    call ndxdec(ds_print, sddisc, sderro, numins)
!
    nbiter = 1
!
end subroutine
