! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine nmconv(mesh, model, ds_material, numeDof, sdnume, listFuncActi, &
                  sddyna, nlDynaDamping, &
                  ds_conv, ds_print, ds_measure, &
                  sddisc, sdcrit, sderro, ds_algopara, ds_algorom, &
                  ds_inout, matass, solveu, ds_system, numeTime, &
                  iterNewt, eta, ds_contact, valinc, solalg, &
                  measse, veasse)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/cfmmcv.h"
#include "asterfort/dierre.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmadev.h"
#include "asterfort/nmcore.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmdivr.h"
#include "asterfort/nmresx.h"
#include "asterfort/nmeceb.h"
#include "asterfort/nmerge.h"
#include "asterfort/nmevcv.h"
#include "asterfort/nmimci.h"
#include "asterfort/nmimr0.h"
#include "asterfort/nmimrv.h"
#include "asterfort/nmlecv.h"
#include "asterfort/nmlerr.h"
#include "asterfort/nmltev.h"
#include "asterfort/nmnkft.h"
#include "asterfort/nmresi.h"
#include "asterfort/nmrvai.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: listFuncActi(*)
    integer(kind=8) :: iterNewt, numeTime
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    real(kind=8) :: eta
    character(len=19) :: sdcrit, sddisc, sdnume
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19) :: matass, solveu
    character(len=19) :: measse(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    character(len=8) :: mesh
    character(len=24) :: numeDof, model
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=24) :: sderro
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_InOut), intent(in) :: ds_inout
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Conv), intent(inout) :: ds_conv
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Global convergence of current Newton iteration
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  MODELE : NOM DU MODELE
! IO  ds_contact       : datastructure for contact management
! IO  ds_conv          : datastructure for convergence management
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_inout         : datastructure for input/output management
! IO  ds_print         : datastructure for printing parameters
! IN  NUMEDD : NUMEROTATION NUME_DDL
! IN  SDNUME : NOM DE LA SD NUMEROTATION
! IN  MATASS : MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  SOLVEU : SOLVEUR
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  ETA    : COEFFICIENT DE PILOTAGE
! In  sddisc          : datastructure for time discretization
! In  numeTime        : index of current time step
! In  iterNewt        : index of current Newton iteration
! IN  SDERRO : GESTION DES ERREURS
! IO  ds_algopara      : datastructure for algorithm parameters
! In  ds_algorom       : datastructure for ROM parameters
! In  ds_system        : datastructure for non-linear system management
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  SDCRIT : SYNTHESE DES RESULTATS DE CONVERGENCE POUR ARCHIVAGE
! In  ds_material      : datastructure for material parameters
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lLineSearch, lNewtonKrylov, lIMPLEX, lContact, lDeborst
    real(kind=8) :: pasMiniElas
    real(kind=8) :: timePrev, timeCurr
    real(kind=8) :: vresi, vchar
    aster_logical :: lEventError, lIterMaxi, dvdebo
    aster_logical :: cvnewt, cvresi
    integer(kind=8) :: nbIter, iterSupp
    integer(kind=8) :: ifm, niv
    real(kind=8) :: line_sear_coef
    integer(kind=8) :: line_sear_iter
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_64')
    end if

! - INITIALISATIONS
    lIterMaxi = ASTER_FALSE
    lEventError = ASTER_FALSE
    cvnewt = ASTER_FALSE
    pasMiniElas = ds_algopara%pas_mini_elas
    ds_algopara%l_swapToElastic = ASTER_FALSE
    line_sear_coef = r8vide()
    line_sear_iter = -1

! - Active functionnalities
    lLineSearch = isfonc(listFuncActi, 'RECH_LINE')
    lNewtonKrylov = isfonc(listFuncActi, 'NEWTON_KRYLOV')
    lIMPLEX = isfonc(listFuncActi, 'IMPLEX')
    lContact = isfonc(listFuncActi, 'CONTACT')
    lDeborst = isfonc(listFuncActi, 'DEBORST')

! - Get time
    timePrev = diinst(sddisc, numeTime-1)
    timeCurr = diinst(sddisc, numeTime)

! - Set values are not affected on rows for residuals loop
    call nmimr0(ds_print, 'RESI')

! - EVENEMENT ERREUR ACTIVE ?
    call nmltev(sderro, 'ERRI', 'NEWT', lEventError)
    if (.not. lEventError) then
!
! ----- EXAMEN DU NOMBRE D'ITERATIONS
!
        call nmlerr(sddisc, 'ITERSUP', paraValeI_=iterSupp)
        if (iterSupp .eq. 0) then
            if (abs(timeCurr-timePrev) .lt. pasMiniElas) then
                nbIter = ds_conv%iter_glob_elas
            else
                nbIter = ds_conv%iter_glob_maxi
            end if
        else
            call nmlerr(sddisc, 'NBITER', paraValeI_=nbIter)
        end if
        lIterMaxi = (iterNewt+1) .ge. nbIter
!
! ----- STATISTIQUES POUR RECHERCHE LINEAIRE
!
        if (lLineSearch) then
            line_sear_coef = ds_conv%line_sear_coef
            line_sear_iter = ds_conv%line_sear_iter
            call nmrvai(ds_measure, 'LineSearch', input_count=line_sear_iter)
        end if

! ----- Compute residuals
        call nmresi(mesh, listFuncActi, ds_material, &
                    numeDof, sdnume, &
                    sddyna, nlDynaDamping, &
                    ds_conv, ds_print, ds_contact, &
                    ds_inout, ds_algorom, ds_system, &
                    matass, numeTime, eta, &
                    valinc, solalg, &
                    veasse, measse, &
                    vresi, vchar)
!
! ----- Evaluate convergence of residuals
!
        call nmcore(sdcrit, sderro, listFuncActi, numeTime, iterNewt, &
                    line_sear_iter, eta, vresi, vchar, ds_conv)
!
! ----- METHODE IMPLEX: CONVERGENCE FORCEE
!
        if (lIMPLEX) then
            call nmeceb(sderro, 'RESI', 'CONV')
        end if
!
! ----- Evaluate convergence of contact and PRED_CONTACT
!
        if (lContact) then
            call cfmmcv(mesh, model, listFuncActi, iterNewt, numeTime, &
                        sddyna, ds_measure, sddisc, sderro, valinc, &
                        solalg, ds_print, ds_contact)
            if (ds_contact%lContStab) then
                if (ds_contact%iContStab .lt. ds_contact%sContStab) then
                    ds_algopara%l_swapToElastic = ASTER_TRUE
                end if
            end if
        end if
!
! ----- Set value of informations in convergence table (residuals are in nmimre)
!
        call nmimrv(ds_print, listFuncActi, iterNewt, line_sear_coef, line_sear_iter, &
                    eta, ds_algorom%eref_rom)
!
! ----- CAPTURE ERREUR EVENTUELLE
!
        call nmltev(sderro, 'ERRI', 'NEWT', lEventError)
        if (.not. lEventError) then
! --------- INFORMATION POUR DEBORST
            call nmlecv(sderro, 'RESI', cvresi)
            call nmerge(sderro, 'DIVE_DEBO', dvdebo)
            if (cvresi .and. dvdebo .and. lDeborst) then
                call utmess('I', 'MECANONLINE2_3')
            end if
!
! --------- EVALUATION DE LA CONVERGENCE DE L'ITERATION DE NEWTON
!
            call nmevcv(sderro, listFuncActi, 'NEWT')
            call nmlecv(sderro, 'NEWT', cvnewt)
!
! --------- ENREGISTRE LES RESIDUS A CETTE ITERATION
!
            call dierre(sddisc, sdcrit, iterNewt)
!
! --------- Check if RESI_GLOB_MAXI increase
!
            call nmdivr(sddisc, sderro, iterNewt)
!
! --------- Check if RESI_GLOB_MAXI is too large
!
            call nmresx(sddisc, sderro, iterNewt)
!
! --------- SI ON A CONVERGE: ON N'A PAS ATTEINT LE NB D'ITERATIONS MAXIMUM
!
            if (cvnewt) then
                lIterMaxi = ASTER_FALSE
            end if
!
! --------- ENREGISTREMENT EVENEMENT MAX ITERATION DE NEWTON
!
            call nmcrel(sderro, 'ITER_MAXI', lIterMaxi)
!
! --------- CALCUL CRITERE DE CONVERGENCE POUR NEWTON-KRYLOV (FORCING-TERM)
!
            if (lNewtonKrylov) then
                call nmnkft(solveu, sddisc, iterNewt)
            end if
        end if
    end if

! - Set iteration number in convergence table
    call nmimci(ds_print, 'ITER_NUME', iterNewt, ASTER_TRUE)

! - MISE A JOUR DE L'INDICATEUR DE SUCCES SUR LES ITERATIONS DE NEWTON
    call nmadev(sddisc, sderro, iterNewt)
!
end subroutine
