! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine nmprma(listFuncActi, &
                  mesh, modelz, caraElem, &
                  ds_material, ds_constitutive, &
                  listLoad, sddyna, &
                  sddisc, numeTime, &
                  ds_algopara, ds_contact, ds_algorom, &
                  ds_print, ds_measure, &
                  hval_incr, hval_algo, hhoField, &
                  hval_meelem, hval_measse, &
                  numeDof, numeDofFixe, &
                  solveu, ds_system, &
                  maprec, matass, &
                  faccvg, ldccvg, condcvg)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use HHO_type
    use NonLinear_module, only: getOption, getMatrType, isMatrUpdate, &
                                isDampMatrCompute, isMassMatrCompute, isRigiMatrCompute, &
                                factorSystem, updateLoadBCMatrix
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/asmari.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/echmat.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmcmat.h"
#include "asterfort/nmelcm.h"
#include "asterfort/nmimck.h"
#include "asterfort/nmmatr.h"
#include "asterfort/nmrenu.h"
#include "asterfort/nmrigi.h"
#include "asterfort/nmxmat.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
!
    integer, intent(in) :: listFuncActi(*)
    character(len=8), intent(in) :: mesh
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: listLoad, sddyna
    character(len=19), intent(in) :: sddisc
    integer, intent(in) :: numeTime
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
    type(HHO_Field), intent(in) :: hhoField
    character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
    character(len=24), intent(inout) :: numeDof
    character(len=24), intent(in) :: numeDofFixe
    character(len=19), intent(in) :: solveu
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: maprec, matass
    integer, intent(out) :: faccvg, ldccvg, condcvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE EN PREDICTION
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
! IO  ds_print         : datastructure for printing parameters
! IN  SDDYNA : SD POUR LA DYNAMIQUE
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! IO  ds_measure       : datastructure for measure and statistics management
! In  listFuncActi   : list of active functionnalities
! In  ds_algorom       : datastructure for ROM parameters
! In  ds_system        : datastructure for non-linear system management
! In  hhoField         : datastructure for HHO
! In  numeTime        : index of current time step
! IN  SOLVEU : SOLVEUR
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
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
! OUT CONDCVG : CODE RETOUR DE LA CONDENSATION STATIQUE
!                -1 : PAS DE CONDENSATION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE LA CONDENSATION
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: phaseType = PRED_EULER
    integer, parameter :: iterNewtPred = 0
    integer :: ifm, niv
    aster_logical :: l_update_matr, l_renumber
    aster_logical :: l_comp_rigi, l_comp_damp, l_asse_rigi
    aster_logical :: l_comp_cont, lMassAssemble
    character(len=16) :: matrType, option_nonlin
    character(len=19) :: contElem, rigid
    integer :: nb_matr, reac_incr
    character(len=6) :: list_matr_type(20)
    character(len=8) :: ksym
    character(len=16) :: list_calc_opti(20), list_asse_opti(20)
    aster_logical :: list_l_asse(20), list_l_calc(20)
    character(len=3) :: mathpc
    character(len=19) :: partit
    aster_logical :: l_contact_adapt, l_cont_cont, lmhpc
    real(kind=8) :: minmat, maxmat, exponent_val
    aster_logical :: ldist
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_35')
    end if

! - Initializations
    nb_matr = 0
    list_matr_type(1:20) = ' '
    faccvg = -1
    ldccvg = -1
    condcvg = -1

! - Active functionnalites
    l_comp_cont = isfonc(listFuncActi, 'ELT_CONTACT')
    l_cont_cont = isfonc(listFuncActi, 'CONT_CONTINU')

! - Renumbering equations ?
    call nmrenu(modelz, listFuncActi, listLoad, &
                ds_measure, ds_contact, numeDof, &
                l_renumber)

! - Get type of matrix
    call getMatrType(phaseType, listFuncActi, sddisc, numeTime, ds_algopara, &
                     matrType, reac_incr_=reac_incr)

! - Update global matrix ?
    call isMatrUpdate(phaseType, matrType, listFuncActi, &
                      sddyna, ds_system, &
                      l_update_matr, &
                      nume_inst_=numeTime, reac_incr_=reac_incr)

! - Select non-linear option for compute matrices
    call getOption(phaseType, listFuncActi, matrType, option_nonlin)

! - Do the damping matrices have to be calculated ?
    call isDampMatrCompute(sddyna, l_renumber, l_comp_damp)

! - Do the mass matrices have to be calculated ?
    call isMassMatrCompute(sddyna, l_update_matr, lMassAssemble)

! - Do the rigidity matrices have to be calculated/assembled ?
    call isRigiMatrCompute(phaseType, &
                           sddyna, numeTime, &
                           l_update_matr, l_comp_damp, &
                           l_comp_rigi, l_asse_rigi)

! - Compute contact elementary matrices
    if (l_comp_cont) then
        call nmchex(hval_meelem, 'MEELEM', 'MEELTC', contElem)
        call nmelcm(mesh, modelz, &
                    ds_material, ds_contact, ds_constitutive, ds_measure, &
                    hval_incr, hval_algo, &
                    contElem)
    end if

! - Compute rigidity elementary matrices / internal forces elementary vectors
    if (l_comp_rigi) then
        call nmrigi(modelz, caraElem, &
                    ds_material, ds_constitutive, &
                    listFuncActi, iterNewtPred, sddyna, ds_measure, ds_system, &
                    hval_incr, hval_algo, hhoField, &
                    option_nonlin, ldccvg)
        if (l_asse_rigi) then
            call nmchex(hval_measse, 'MEASSE', 'MERIGI', rigid)
            call asmari(ds_system, hval_meelem, listLoad, rigid)
        end if
    end if

! - No error => continue
    if (ldccvg .ne. 1) then

! ----- Update elementary matrices for loads and boundary conditions (undead cases)
        if (matrType .ne. 'EXTRAPOLE') then
            call updateLoadBCMatrix(listFuncActi, listLoad, &
                                    sddisc, numeTime, &
                                    modelZ, caraElem, &
                                    ds_material, ds_constitutive, &
                                    hval_incr, hval_algo, &
                                    hval_meelem)
        end if

! ----- Compute damping (Rayleigh) elementary matrices
        if (l_comp_damp) then
            call nmcmat('MEAMOR', ' ', ' ', ASTER_TRUE, &
                        ASTER_TRUE, nb_matr, list_matr_type, list_calc_opti, list_asse_opti, &
                        list_l_calc, list_l_asse)
        end if

! ----- Assemble mass matrix
        if (lMassAssemble) then
            call nmcmat('MEMASS', ' ', ' ', ASTER_FALSE, &
                        ASTER_TRUE, nb_matr, list_matr_type, list_calc_opti, list_asse_opti, &
                        list_l_calc, list_l_asse)
            ASSERT(l_update_matr)
        end if

! ----- Compute and assemble matrices
        if (nb_matr .gt. 0) then
            call nmxmat(modelz, ds_material, caraElem, &
                        ds_constitutive, sddisc, numeTime, &
                        hval_incr, hval_algo, listLoad, &
                        numeDof, numeDofFixe, ds_measure, &
                        nb_matr, list_matr_type, list_calc_opti, &
                        list_asse_opti, list_l_calc, list_l_asse, &
                        hval_meelem, hval_measse, ds_system)
        end if

! ----- Compute global matrix of system
        if (l_update_matr) then
            call nmmatr(phaseType, listFuncActi, listLoad, numeDof, sddyna, &
                        numeTime, ds_contact, hval_meelem, hval_measse, matass)
            call dismoi('TYPE_MATRICE', matass, 'MATR_ASSE', repk=ksym)
            select case (ksym(1:7))
            case ('SYMETRI')
                matrType(12:16) = '(SYM)'
            case ('NON_SYM')
                matrType(10:16) = '(NOSYM)'
            case default
                ASSERT(.false.)
            end select
            call nmimck(ds_print, 'MATR_ASSE', matrType, ASTER_TRUE)
        else
            call nmimck(ds_print, 'MATR_ASSE', ' ', ASTER_FALSE)
        end if

! ----- Change matrix for contact (BEURK)
        if (l_cont_cont) then
            minmat = 0.0
            maxmat = 0.0
            exponent_val = 0.0
!   -- Avant la factorisation et pour le cas ou il y a du contact continu avec adaptation de
!      coefficient
!   -- On cherche le coefficient optimal pour eviter une possible singularite de matrice
!   -- La valeur est estimee une seule fois a la premiere prediction du premier pas de
!      temps pour l'etape de calcul
!   -- Cette valeur estimee est passee directement a mmchml_c sans passer par mmalgo car
!   -- a la premiere iteration on ne passe pas par mmalgo
            l_contact_adapt = cfdisl(ds_contact%sdcont_defi, 'EXIS_ADAP')
!            write (6,*) "l_contact_adapt", &
!                l_contact_adapt,ds_contact%update_init_coefficient
            if ((nint(ds_contact%update_init_coefficient) .eq. 0) .and. l_contact_adapt) then
                call dismoi('MATR_HPC', matass, 'MATR_ASSE', repk=mathpc)
                lmhpc = mathpc .eq. 'OUI'
                call dismoi('PARTITION', modelz, 'MODELE', repk=partit)
                ldist = partit .ne. ' '
                call echmat(matass, ldist, lmhpc, minmat, maxmat)
                ds_contact%max_coefficient = maxmat
                if (abs(log(minmat)) .ge. r8prem()) then

                    if (abs(log(maxmat))/abs(log(minmat)) .lt. 4.0d0) then
!                     Le rapport d'arete max/min est
!  un bon compromis pour initialiser le coefficient
                        ds_contact%estimated_coefficient = &
                            ((1.D3*ds_contact%arete_max)/(1.D-2*ds_contact%arete_min))
                        ds_contact%update_init_coefficient = 1.0d0
                    else
                        exponent_val = min(abs(log(minmat)), abs(log(maxmat)))/10.d0
                        ds_contact%estimated_coefficient = 10.d0**(exponent_val)
                        ds_contact%update_init_coefficient = 1.0d0
                    end if
                else
                    ds_contact%estimated_coefficient = 1.d16*ds_contact%arete_min
                    ds_contact%update_init_coefficient = 1.0d0
                end if
!             write (6,*) "min,max,coef estime,abs(log(maxmat))/abs(log(minmat))", &
!                 minmat,maxmat,ds_contact%estimated_coefficient,abs(log(maxmat))/abs(log(minmat))
            end if
        end if

! ----- Factorization of global matrix of system
        if (l_update_matr) then
            call factorSystem(listFuncActi, ds_measure, ds_algorom, &
                              numeDof, solveu, maprec, matass, &
                              faccvg)
        end if
    end if
!
end subroutine
