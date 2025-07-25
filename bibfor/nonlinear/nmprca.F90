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
subroutine nmprca(modele, numedd, numfix, ds_material, carele, &
                  ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                  fonact, ds_print, ds_measure, ds_algorom, sddisc, numins, &
                  valinc, solalg, matass, maprec, ds_contact, &
                  sddyna, nlDynaDamping, &
                  meelem, measse, veelem, veasse, &
                  depest, ldccvg, faccvg, rescvg)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmacin.h"
#include "asterfort/nmassd.h"
#include "asterfort/nmchex.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmprma.h"
#include "asterfort/nmreso.h"
#include "asterfort/vtzero.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/romAlgoNLSystemSolve.h"
#include "asterfort/nonlinLoadDirichletCompute.h"
#include "asterfort/romCoefComputeFromField.h"
!
    integer(kind=8) :: fonact(*)
    integer(kind=8) :: numins, ldccvg, faccvg, rescvg
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19) :: maprec, matass
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19) :: lischa, solveu, sddisc
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=24) :: modele, carele
    character(len=24) :: numedd, numfix
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19) :: veelem(*), veasse(*)
    character(len=19) :: meelem(*), measse(*)
    character(len=19) :: solalg(*), valinc(*)
    character(len=19) :: depest
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PREDICTION - DEPL. DONNE)
!
! PROJECTION DU CHAMP DONNE SUR L'ESPACE DES CONDITIONS AUX LIMITES
! CINEMATIQUEMENT ADMISSIBLES
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
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IO  ds_print         : datastructure for printing parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algorom       : datastructure for ROM parameters
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MATASS : MATRICE ASSEMBLEE
! IN  MAPREC : MATRICE DE PRECONDITIONNEMENT (GCPC)
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  DEPEST : DEPLACEMENT ESTIME
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
    integer(kind=8) :: neq, i
    character(len=19) :: depso1, depso2, cncine
    character(len=19) :: solu1, solu2, cndonn, cnpilo, cncind
    aster_logical :: l_update_redu, l_rom
    character(len=24) :: mata24, vect24, cncine24
    real(kind=8), pointer :: dep1(:) => null()
    real(kind=8), pointer :: dep2(:) => null()
    real(kind=8), pointer :: sol1(:) => null()
    real(kind=8), pointer :: sol2(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
    real(kind=8), pointer :: v_gamma(:) => null()
    real(kind=8), pointer :: v_vect(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    solu1 = '&&CNPART.CHP2'
    solu2 = '&&CNPART.CHP3'
    cndonn = '&&CNCHAR.DONN'
    cnpilo = '&&CNCHAR.PILO'
    cncind = '&&CNCHAR.CINE'
    call vtzero(solu1)
    call vtzero(solu2)
    call vtzero(cndonn)
    call vtzero(cnpilo)
    call vtzero(cncind)
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    ldccvg = -1
    faccvg = -1
    l_rom = isfonc(fonact, 'ROM')
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    call nmchex(solalg, 'SOLALG', 'DEPSO1', depso1)
    call nmchex(solalg, 'SOLALG', 'DEPSO2', depso2)
    call nmchex(veasse, 'VEASSE', 'CNCINE', cncine)
!
! --- CALCUL DE LA MATRICE GLOBALE
!
    call nmprma(fonact, &
                modele, carele, &
                ds_material, ds_constitutive, &
                lischa, sddyna, nlDynaDamping, &
                sddisc, numins, &
                ds_algopara, ds_contact, ds_algorom, &
                ds_print, ds_measure, &
                valinc, solalg, &
                meelem, measse, &
                numedd, numfix, &
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
! - Compute values of Dirichlet conditions
!
    call nonlinLoadDirichletCompute(lischa, modele, numedd, &
                                    ds_measure, matass, depest, &
                                    veelem, veasse)
!
! - Evaluate second member for Dirichlet loads (AFFE_CHAR_MECA)
!
    call nmassd(fonact, veasse, cnpilo, cndonn)
!
! --- PRISE EN COMPTE DES CL ELIMINEES
!
    call copisd('CHAMP_GD', 'V', cncine, cncind)
    call nmacin(fonact, matass, depso1, cncind)

! - Timer for solver
    call nmtime(ds_measure, 'Init', 'Solve')
    call nmtime(ds_measure, 'Launch', 'Solve')
!
! - Solving
!
    if (l_rom) then
        l_update_redu = ASTER_FALSE
        rescvg = 0
        mata24 = matass
        vect24 = cndonn
        cncine24 = cncine
        call romAlgoNLSystemSolve(mata24, vect24, cncine24, ds_algorom, solu1, l_update_redu)
    else
        call nmreso(fonact, cndonn, cnpilo, cncind, solveu, &
                    maprec, matass, solu1, solu2, rescvg)
    end if

! - End of timer for solver
    call nmtime(ds_measure, 'Stop', 'Solve')
    call nmrinc(ds_measure, 'Solve')
!
! --- ERREUR SANS POSSIBILITE DE CONTINUER
!
    if (rescvg .eq. 1) goto 999
!
! --- CORRECTION DU DEPLACEMENT DONNE POUR LE RENDRE
! --- CINEMATIQUEMENT ADMISSIBLE
!
    call jeveuo(numedd(1:14)//'.NUME.DELG', 'L', vi=delg)
    call jeveuo(solu1(1:19)//'.VALE', 'L', vr=sol1)
    call jeveuo(solu2(1:19)//'.VALE', 'L', vr=sol2)
    call jeveuo(depso1(1:19)//'.VALE', 'E', vr=dep1)
    call jeveuo(depso2(1:19)//'.VALE', 'E', vr=dep2)
!
! --- LES LAGRANGES NE SONT PAS MODIFIES
!
    do i = 1, neq
        if (delg(i) .eq. 0) then
            dep1(i) = dep1(i)+sol1(i)
            dep2(i) = sol2(i)
        end if
    end do
!
! - For HROM: update reduced coordinates
!
    if (l_rom) then
        AS_ALLOCATE(vr=v_vect, size=ds_algorom%ds_empi%nbMode)
        call jeveuo(ds_algorom%gamma, 'E', vr=v_gamma)
        call romCoefComputeFromField(ds_algorom%ds_empi, dep1, v_vect)
        v_gamma = v_gamma+v_vect
        AS_DEALLOCATE(vr=v_vect)
    end if
!
999 continue
!
end subroutine
