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
subroutine nmconv(noma, modele, ds_material, numedd, sdnume, fonact, &
                  sddyna, nlDynaDamping, &
                  ds_conv, ds_print, ds_measure, &
                  sddisc, sdcrit, sderro, ds_algopara, ds_algorom, &
                  ds_inout, matass, solveu, ds_system, numins, &
                  iterat, eta, ds_contact, valinc, solalg, &
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
    integer(kind=8) :: fonact(*)
    integer(kind=8) :: iterat, numins
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    real(kind=8) :: eta
    character(len=19) :: sdcrit, sddisc, sdnume
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19) :: matass, solveu
    character(len=19) :: measse(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    character(len=8) :: noma
    character(len=24) :: numedd, modele
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
! IN  ITERAT : NUMERO D'ITERATION
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  MEASSE : VARIABLE CHAPEAU POUR NOM DES MATR_ASSE
! IN  ETA    : COEFFICIENT DE PILOTAGE
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
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
    aster_logical :: lreli, lnkry, limpex, lcont
    real(kind=8) :: r8bid
    real(kind=8) :: pasmin
    real(kind=8) :: instam, instap
    real(kind=8) :: vresi, vchar
    aster_logical :: lerror, itemax, dvdebo
    aster_logical :: cvnewt, cvresi
    integer(kind=8) :: nbiter, itesup
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
!
! --- INITIALISATIONS
!
    itemax = .false.
    lerror = .false.
    cvnewt = .false.
    pasmin = ds_algopara%pas_mini_elas
    ds_algopara%l_swapToElastic = ASTER_FALSE
    line_sear_coef = r8vide()
    line_sear_iter = -1
!
! --- FONCTIONNALITES ACTIVEES
!
    lreli = isfonc(fonact, 'RECH_LINE')
    lnkry = isfonc(fonact, 'NEWTON_KRYLOV')
    limpex = isfonc(fonact, 'IMPLEX')
    lcont = isfonc(fonact, 'CONTACT')
!
! --- INSTANTS
!
    instam = diinst(sddisc, numins-1)
    instap = diinst(sddisc, numins)
!
! - Set values are not affected on rows for residuals loop
!
    call nmimr0(ds_print, 'RESI')
!
! --- EVENEMENT ERREUR ACTIVE ?
!
    call nmltev(sderro, 'ERRI', 'NEWT', lerror)
    if (.not. lerror) then
!
! ----- EXAMEN DU NOMBRE D'ITERATIONS
!
        call nmlerr(sddisc, 'L', 'ITERSUP', r8bid, itesup)
        if (itesup .eq. 0) then
            if (abs(instap-instam) .lt. pasmin) then
                nbiter = ds_conv%iter_glob_elas
            else
                nbiter = ds_conv%iter_glob_maxi
            end if
        else
            call nmlerr(sddisc, 'L', 'NBITER', r8bid, nbiter)
        end if
        itemax = (iterat+1) .ge. nbiter
!
! ----- STATISTIQUES POUR RECHERCHE LINEAIRE
!
        if (lreli) then
            line_sear_coef = ds_conv%line_sear_coef
            line_sear_iter = ds_conv%line_sear_iter
            call nmrvai(ds_measure, 'LineSearch', input_count=line_sear_iter)
        end if

! ----- Compute residuals
        call nmresi(noma, fonact, ds_material, &
                    numedd, sdnume, &
                    sddyna, nlDynaDamping, &
                    ds_conv, ds_print, ds_contact, &
                    ds_inout, ds_algorom, ds_system, &
                    matass, numins, eta, &
                    valinc, solalg, &
                    veasse, measse, &
                    vresi, vchar)
!
! ----- Evaluate convergence of residuals
!
        call nmcore(sdcrit, sderro, fonact, numins, iterat, &
                    line_sear_iter, eta, vresi, vchar, ds_conv)
!
! ----- METHODE IMPLEX: CONVERGENCE FORCEE
!
        if (limpex) then
            call nmeceb(sderro, 'RESI', 'CONV')
        end if
!
! ----- Evaluate convergence of contact and PRED_CONTACT
!
        if (lcont) then
            call cfmmcv(noma, modele, fonact, iterat, numins, &
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
        call nmimrv(ds_print, fonact, iterat, line_sear_coef, line_sear_iter, &
                    eta, ds_algorom%eref_rom)
!
! ----- CAPTURE ERREUR EVENTUELLE
!
        call nmltev(sderro, 'ERRI', 'NEWT', lerror)
        if (.not. lerror) then
!
! --------- INFORMATION POUR DEBORST
!
            call nmlecv(sderro, 'RESI', cvresi)
            call nmerge(sderro, 'DIVE_DEBO', dvdebo)
            if (cvresi .and. dvdebo) then
                call utmess('I', 'MECANONLINE2_3')
            end if
!
! --------- EVALUATION DE LA CONVERGENCE DE L'ITERATION DE NEWTON
!
            call nmevcv(sderro, fonact, 'NEWT')
            call nmlecv(sderro, 'NEWT', cvnewt)
!
! --------- ENREGISTRE LES RESIDUS A CETTE ITERATION
!
            call dierre(sddisc, sdcrit, iterat)
!
! --------- Check if RESI_GLOB_MAXI increase
!
            call nmdivr(sddisc, sderro, iterat)
!
! --------- Check if RESI_GLOB_MAXI is too large
!
            call nmresx(sddisc, sderro, iterat)
!
! --------- SI ON A CONVERGE: ON N'A PAS ATTEINT LE NB D'ITERATIONS MAXIMUM
!
            if (cvnewt) then
                itemax = .false.
            end if
!
! --------- ENREGISTREMENT EVENEMENT MAX ITERATION DE NEWTON
!
            call nmcrel(sderro, 'ITER_MAXI', itemax)
!
! --------- CALCUL CRITERE DE CONVERGENCE POUR NEWTON-KRYLOV (FORCING-TERM)
!
            if (lnkry) then
                call nmnkft(solveu, sddisc, iterat)
            end if
        end if
    end if
!
! - Set iteration number in convergence table
!
    call nmimci(ds_print, 'ITER_NUME', iterat, .true._1)
!
! --- MISE A JOUR DE L'INDICATEUR DE SUCCES SUR LES ITERATIONS DE NEWTON
!
    call nmadev(sddisc, sderro, iterat)
!
end subroutine
