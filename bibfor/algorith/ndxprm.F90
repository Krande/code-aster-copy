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
subroutine ndxprm(modelz, ds_material, carele, ds_constitutive, ds_algopara, &
                  lischa, numedd, solveu, ds_system, sddisc, &
                  sddyna, nlDynaDamping, &
                  ds_measure, nume_inst, list_func_acti, &
                  valinc, solalg, meelem, measse, &
                  maprec, matass, faccvg, ldccvg)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module, only: isDampMatrUpdate, compDampMatrix
    use NonLinear_module, only: updateLoadBCMatrix
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/asmari.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/ndxmat.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmrigi.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/preres.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    integer(kind=8), intent(in) :: list_func_acti(*), nume_inst
    character(len=*) :: modelz
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24) :: carele
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24) :: numedd
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19) :: sddisc, lischa, solveu
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19) :: solalg(*), valinc(*)
    character(len=19) :: meelem(*), measse(*)
    character(len=19) :: maprec, matass
    integer(kind=8) :: faccvg, ldccvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL - UTILITAIRE)
!
! CALCUL DE LA MATRICE GLOBALE EN PREDICTION - CAS EXPLICITE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : LISTE DES CHARGES
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  SOLVEU : SOLVEUR
! In  ds_system        : datastructure for non-linear system management
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  NUMINS : NUMERO D'INSTANT
! IN  ITERAT : NUMERO D'ITERATION
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  MEELEM : VARIABLE CHAPEAU POUR NOM DES MATR_ELEM
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
    integer(kind=8), parameter :: iterNewtPred = 0
    character(len=16), parameter :: nonLinearOption = 'RIGI_MECA_TANG'
    aster_logical, parameter :: l_renumber = ASTER_FALSE
    aster_logical :: l_update_matr, l_first_step
    aster_logical :: lRigiCompute, lDampMatrUpdate, lRigiAssemble
    aster_logical :: lshima, lprmo
    character(len=16) :: explMatrType
    integer(kind=8) :: ifm, niv, ibid
    character(len=19) :: rigid
    real(kind=8) :: time
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... CALCUL MATRICE'
    end if

! - Initializations
    faccvg = -1
    ldccvg = -1

! - Time
    time = diinst(sddisc, nume_inst-1)

! - Active functionnalities
    lprmo = ndynlo(sddyna, 'PROJ_MODAL')
    lshima = ndynlo(sddyna, 'COEF_MASS_SHIFT')

! - First step ?
    l_first_step = nume_inst .le. 1

! - Get type of matrix
    explMatrType = ds_algopara%matrix_pred
    if (explMatrType .ne. 'TANGENTE') then
        call utmess('F', 'MECANONLINE5_1')
    end if

! - Update global matrix ?
    l_update_matr = ASTER_FALSE
    if (l_first_step) then
        if (lprmo) then
            l_update_matr = ASTER_FALSE
        else
            l_update_matr = ASTER_TRUE
        end if
    end if

! - Do the damping matrices need to be updated ?
    call isDampMatrUpdate(nlDynaDamping, l_renumber, lDampMatrUpdate)

! - Do the rigidity matrices have to be calculated/assembled ?
    lRigiCompute = ASTER_FALSE
    lRigiAssemble = ASTER_FALSE
    if (lDampMatrUpdate) then
        lRigiCompute = ASTER_TRUE
    end if
    if (lshima .and. l_first_step) then
        lRigiCompute = ASTER_TRUE
        lRigiAssemble = ASTER_TRUE
    end if

! - Compute rigidity elementary matrices / internal forces elementary vectors
    if (lRigiCompute) then
        call nmrigi(modelz, carele, &
                    ds_material, ds_constitutive, &
                    list_func_acti, iterNewtPred, sddyna, ds_measure, ds_system, &
                    valinc, solalg, &
                    nonLinearOption, ldccvg)
        if (lRigiAssemble) then
            call nmchex(measse, 'MEASSE', 'MERIGI', rigid)
            call asmari(ds_system, meelem, lischa, rigid)
        end if
    end if

! - Update elementary matrices for loads and boundary conditions (undead cases)
    call updateLoadBCMatrix(list_func_acti, lischa, &
                            sddisc, nume_inst, &
                            modelZ, carele, &
                            ds_material, ds_constitutive, &
                            valinc, solalg, &
                            meelem)

! - Compute damping matrix
    if (lDampMatrUpdate) then
        call compDampMatrix(modelz, carele, &
                            ds_material, ds_constitutive, &
                            time, lischa, numedd, nlDynaDamping, &
                            ds_system, valinc, meelem, sddyna)
    end if

! - No error => continue
    if (ldccvg .ne. 1) then
! ----- Compute global matrix of system
        if (l_update_matr) then
            call ndxmat(list_func_acti, lischa, numedd, sddyna, nume_inst, &
                        meelem, measse, matass)
        end if
! ----- Factorization of global matrix of system
        if (l_update_matr) then
            call nmtime(ds_measure, 'Init', 'Factor')
            call nmtime(ds_measure, 'Launch', 'Factor')
            call preres(solveu, 'V', faccvg, maprec, matass, &
                        ibid, -9999)
            call nmtime(ds_measure, 'Stop', 'Factor')
            call nmrinc(ds_measure, 'Factor')
        end if
    end if
!
end subroutine
