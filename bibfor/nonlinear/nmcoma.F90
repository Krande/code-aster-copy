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
subroutine nmcoma(listFuncActi, &
                  modelz, caraElem, &
                  ds_material, ds_constitutive, &
                  listLoad, sddyna, nlDynaDamping, &
                  sddisc, numeTime, iterNewt, &
                  ds_algopara, ds_contact, ds_algorom, &
                  ds_print, ds_measure, &
                  hval_incr, hval_algo, &
                  hval_meelem, hval_measse, &
                  numeDof, numeDofFixe, sdnume, &
                  solveu, ds_system, &
                  maprec, matrAsse, &
                  faccvg, ldccvg)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
    use NonLinear_module, only: getOption, getMatrType, isMatrUpdate, &
                                isRigiMatrCompute, isInteVectCompute, &
                                factorSystem, updateLoadBCMatrix
    use NonLinearDyna_module, only: isDampMatrUpdate, isMassMatrAssemble, &
                                    compDampMatrix, asseMassMatrix
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmari.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmelcm.h"
#include "asterfort/nmimck.h"
#include "asterfort/nmmatr.h"
#include "asterfort/nmrenu.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/nonlinIntForce.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: listLoad, sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeTime, iterNewt
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
    character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
    character(len=24), intent(inout) :: numeDof
    character(len=24), intent(in) :: numeDofFixe
    character(len=19), intent(in) :: solveu, sdnume
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: maprec
    character(len=19), intent(inout) :: matrAsse
    integer(kind=8), intent(out) :: faccvg, ldccvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE EN CORRECTION
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
! IO  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IO  ds_print         : datastructure for printing parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  listFuncActi   : list of active functionnalities
! In  ds_algorom       : datastructure for ROM parameters
! In  ds_system        : datastructure for non-linear system management
! In  numeTime        : index of current time step
! In  iterNewt        : index of current Newton iteration
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! OUT LFINT  : .TRUE. SI FORCES INTERNES CALCULEES
! OUT MATASS : MATRICE DE RESOLUTION ASSEMBLEE
! OUT MAPREC : MATRICE DE RESOLUTION ASSEMBLEE - PRECONDITIONNEMENT
! OUT FACCVG : CODE RETOUR FACTORISATION MATRICE GLOBALE
!                -1 : PAS DE FACTORISATION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : MATRICE SINGULIERE
!                 2 : ERREUR LORS DE LA FACTORISATION
!                 3 : ON NE SAIT PAS SI SINGULIERE
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 2 : ERREUR SUR LA NON VERIF. DE CRITERES PHYSIQUES
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = CORR_NEWTON
    integer(kind=8) :: ifm, niv
    aster_logical :: l_update_matr, l_renumber
    aster_logical :: lRigiCompute, lDampMatrUpdate, lRigiAssemble
    aster_logical :: lMassAssemble
    aster_logical :: lFintCompute
    aster_logical :: lContCompute
    character(len=16) :: matrType, nonLinearOption
    character(len=19) :: contElem, rigid
    integer(kind=8) :: reac_iter
    character(len=8) :: answer
    real(kind=8) :: time
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_68')
    end if

! - Initializations
    faccvg = -1
    ldccvg = -1

! - Time
    time = diinst(sddisc, numeTime-1)

! - Active functionnalites
    lContCompute = isfonc(listFuncActi, 'ELT_CONTACT')

! - Renumbering equations ?
    call nmrenu(modelz, listFuncActi, listLoad, &
                ds_measure, ds_contact, numeDof, &
                l_renumber)

! - Get type of matrix
    call getMatrType(phaseType, listFuncActi, sddisc, numeTime, ds_algopara, &
                     matrType, reac_iter_=reac_iter)

! - Update global matrix ?
    call isMatrUpdate(phaseType, matrType, listFuncActi, &
                      nlDynaDamping, ds_system, &
                      l_update_matr, &
                      iter_newt_=iterNewt, reac_iter_=reac_iter)

! - Select non-linear option for compute matrices
    call getOption(phaseType, listFuncActi, matrType, nonLinearOption, l_update_matr)

! - Do the damping matrices need to be updated ?
    call isDampMatrUpdate(nlDynaDamping, l_renumber, lDampMatrUpdate)

! - Do the mass matrix has to be assembled ?
    call isMassMatrAssemble(listFuncActi, l_update_matr, lMassAssemble)

! - Do the rigidity matrices have to be calculated/assembled ?
    call isRigiMatrCompute(phaseType, &
                           sddyna, numeTime, &
                           l_update_matr, lDampMatrUpdate, &
                           lRigiCompute, lRigiAssemble)

! - Do the internal forces vector has to be computed ?
    call isInteVectCompute(phaseType, listFuncActi, &
                           nonLinearOption, iterNewt, &
                           lRigiCompute, lFintCompute)

! - Compute internal forces elementary vectors
    if (lFintCompute) then
        call nonlinIntForce(phaseType, &
                            modelZ, caraElem, &
                            listFuncActi, iterNewt, sdnume, &
                            ds_material, ds_constitutive, &
                            ds_system, ds_measure, &
                            hval_incr, hval_algo, &
                            ldccvg, &
                            sddyna_=sddyna)
    end if

! - No error => continue
    if (ldccvg .ne. 1) then
! ----- Compute contact elementary matrices
        if (lContCompute) then
            call nmchex(hval_meelem, 'MEELEM', 'MEELTC', contElem)
            call nmelcm(modelZ, &
                        ds_material, ds_contact, &
                        ds_constitutive, ds_measure, &
                        hval_incr, hval_algo, &
                        contElem)
        end if

! ----- Update elementary matrices for loads and boundary conditions (undead cases)
        call updateLoadBCMatrix(listFuncActi, listLoad, &
                                sddisc, numeTime, &
                                modelZ, caraElem, &
                                ds_material, ds_constitutive, &
                                hval_incr, hval_algo, &
                                hval_meelem)

! ----- Assembly rigidity matrix
        if (lRigiAssemble) then
            call nmchex(hval_measse, 'MEASSE', 'MERIGI', rigid)
            call asmari(ds_system, hval_meelem, listLoad, rigid)
        end if

! ----- Assemble mass matrix
        if (lMassAssemble) then
            call asseMassMatrix(listLoad, &
                                numeDof, numeDofFixe, &
                                hval_meelem, hval_measse)
            ASSERT(l_update_matr)
        end if

! ----- Compute damping matrix
        if (lDampMatrUpdate) then
            call compDampMatrix(modelz, caraElem, &
                                ds_material, ds_constitutive, &
                                time, listLoad, numeDof, nlDynaDamping, &
                                ds_system, hval_incr, hval_meelem, sddyna)
        end if

! ----- Compute global matrix of system
        if (l_update_matr) then
            call nmmatr(phaseType, listFuncActi, listLoad, numeDof, &
                        sddyna, nlDynaDamping, &
                        numeTime, ds_contact, hval_meelem, hval_measse, &
                        matrAsse)
            call dismoi('TYPE_MATRICE', matrAsse, 'MATR_ASSE', repk=answer)
            select case (answer(1:7))
            case ('SYMETRI')
                matrType(12:16) = '(SYM)'
            case ('NON_SYM')
                matrType(10:16) = '(NOSYM)'
            case default
                ASSERT(ASTER_FALSE)
            end select
            call nmimck(ds_print, 'MATR_ASSE', matrType, ASTER_TRUE)
        else
            call nmimck(ds_print, 'MATR_ASSE', ' ', ASTER_FALSE)
        end if

! ----- Factorization of global matrix of system
        if (l_update_matr) then
            call factorSystem(listFuncActi, ds_measure, ds_algorom, &
                              numeDof, solveu, maprec, matrAsse, &
                              faccvg)
        end if
    end if
!
end subroutine
