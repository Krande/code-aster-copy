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
subroutine nmcese(modele, numedd, ds_material, carele, &
                  ds_constitutive, ds_contact, lischa, fonact, ds_measure, &
                  iterat, sdnume, sdpilo, valinc, solalg, &
                  veelem, veasse, offset, typsel, sddisc, &
                  licite, rho, eta, etaf, criter, &
                  ldccvg, pilcvg, matass, ds_system)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/exixfe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmceai.h"
#include "asterfort/nmceni.h"
#include "asterfort/nmcere.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmrcyc.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
!
    integer :: fonact(*)
    integer :: iterat
    real(kind=8) :: rho, offset, eta(2)
    character(len=19) :: lischa, sdnume, sdpilo, sddisc, matass
    character(len=24) :: modele, numedd, carele
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19) :: veelem(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    character(len=24) :: typsel
    integer :: licite(2)
    integer :: ldccvg, pilcvg
    real(kind=8) :: etaf, criter
    type(NL_DS_System), intent(in) :: ds_system
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE)
!
! Select solution
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_contact       : datastructure for contact management
! IN  LISCHA : LISTE DES CHARGES
! IN  SDPILO : SD PILOTAGE
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! IN  SDNUME : SD NUMEROTATION
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  OFFSET : DECALAGE DE ETA_PILOTAGE EN FONCTION DE RHO
! IN  TYPSEL : TYPE DE SELECTION PILOTAGE
!                'ANGL_INCR_DEPL'
!                'NORM_INCR_DEPL'
!                'RESIDU'
! In  sddisc           : datastructure for time discretization
! IN  LICITE : CODE RETOUR PILOTAGE DES DEUX PARAMETRES DE PILOTAGE
! IN  RHO    : PARAMETRE DE RECHERCHE_LINEAIRE
! IN  ETA    : LES DEUX PARAMETRES DE PILOTAGE
! OUT ETAF   : PARAMETRE DE PILOTAGE FINALEMENT CHOISI
! OUT CRITER : VALEUR DU CRITERE DE COMPARAISON
!                ANGL_INCR_DEPL
!                NORM_INCR_DEPL
!                RESIDU
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
! I/O PILCVG : CODE DE CONVERGENCE POUR LE PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
! IN  MATASS : SD MATRICE ASSEMBLEE
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ldccv(2), ierm, indic, sel
    real(kind=8) :: f(2)
    character(len=8) :: choix, txt
    character(len=19) :: depold, depdel, deppr1, deppr2
    character(len=24) :: typpil
    aster_logical :: swloun, isxfe
    aster_logical :: switch, mixte
    real(kind=8) :: miincr, miresi, fnid(2)
    real(kind=8), pointer :: plir(:) => null()
    character(len=24), pointer :: pltk(:) => null()
    real(kind=8), parameter :: contra = 0.1d0
    real(kind=8), parameter :: precyc = 5.d-2
!
! --------------------------------------------------------------------------------------------------
!

!
! --- LE CALCUL DE PILOTAGE A FORCEMENT ETE REALISE
!
    ASSERT(pilcvg .ge. 0)
!
! --- INITIALISATIONS
!
    call jeveuo(sdpilo(1:19)//'.PLTK', 'L', vk24=pltk)
    typpil = pltk(1)
    f(:) = 0.d0
    ldccv(:) = -1
    ldccvg = -1
!
! --- STRATEGIE MIXTE BASEE SUR LE CONTRASTE DES CRITERES DE CHOIX
!
    mixte = typsel .eq. 'MIXTE'
    switch = .false.
!
! --- VERIFICATION DE LA COMPATIBILITE
!
    if (mixte) then
        call utdidt('L', sddisc, 'ECHE', 'CHOIX_SOLU_PILO', valk_=choix)
        if (choix .eq. 'AUTRE') then
            call utmess('F', 'MECANONLINE_62')
        end if
    end if
!
! --- STRATEGIE BASEE SUR LES TECHNIQUES EVENT-DRIVEN 'AUTRE_PILOTAGE'
!
    swloun = .false.
!
! --- ETONNANT QUE X-FEM SE GLISSE A CE NIVEAU DU CODE (ET PLUS AVANT)
!
    call exixfe(modele, ierm)
    isxfe = (ierm .eq. 1)
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(solalg, 'SOLALG', 'DEPPR1', deppr1)
    call nmchex(solalg, 'SOLALG', 'DEPPR2', deppr2)
    call nmchex(solalg, 'SOLALG', 'DEPOLD', depold)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
!
! --- SELECTION SELON LA METHODE CHOISIE: ANGL_INCR_DEPL
!
    if (typsel .eq. 'ANGL_INCR_DEPL') then
!
        if (typpil .eq. 'LONG_ARC' .or. typpil .eq. 'SAUT_LONG_ARC') then
            call jeveuo(sdpilo(1:19)//'.PLIR', 'L', vr=plir)
            swloun = plir(1)*plir(6) .lt. 0.d0
        end if
        call nmceai(numedd, depdel, deppr1, deppr2, depold, &
                    sdpilo, rho, eta(1), isxfe, f(1), &
                    indic)
        call nmceai(numedd, depdel, deppr1, deppr2, depold, &
                    sdpilo, rho, eta(2), isxfe, f(2), &
                    indic)
        if (indic .eq. 0) then
            call nmceni(numedd, depdel, deppr1, deppr2, rho, &
                        sdpilo, eta(1), isxfe, f(1))
            call nmceni(numedd, depdel, deppr1, deppr2, rho, &
                        sdpilo, eta(2), isxfe, f(2))
        end if
        goto 500
    end if
!
! --- SELECTION SELON LA METHODE CHOISIE: NORM_INCR_DEPL OU MIXTE
!
    if (typsel .eq. 'NORM_INCR_DEPL' .or. mixte) then
        call nmceni(numedd, depdel, deppr1, deppr2, rho, &
                    sdpilo, eta(1), isxfe, f(1))
        call nmceni(numedd, depdel, deppr1, deppr2, rho, &
                    sdpilo, eta(2), isxfe, f(2))
!
! ----- SI STRATEGIE MIXTE : EXAMEN DU CONTRASTE
!
        if (mixte) then
            miincr = min(f(1), f(2))/max(f(1), f(2))
            if (miincr .le. contra) goto 600
!
! ------- ECHEC DU CONTRASTE: ON ENCHAINE PAR LA SELECTION RESIDU
!
            fnid(1) = f(1)
            fnid(2) = f(2)
        else
            goto 500
        end if
    end if
!
! - Compute residual
!
    if (typsel .eq. 'RESIDU' .or. mixte) then
        call nmcere(modele, numedd, ds_material, carele, &
                    ds_constitutive, ds_contact, lischa, fonact, ds_measure, &
                    iterat, sdnume, valinc, solalg, veelem, &
                    veasse, offset, rho, eta(1), f(1), &
                    ldccv(1), ds_system, matass)
        call nmcere(modele, numedd, ds_material, carele, &
                    ds_constitutive, ds_contact, lischa, fonact, ds_measure, &
                    iterat, sdnume, valinc, solalg, veelem, &
                    veasse, offset, rho, eta(2), f(2), &
                    ldccv(2), ds_system, matass)
    end if
!
! - Have a look on "contrast"
!
    if (typsel .eq. 'RESIDU' .or. mixte) then
        if (mixte) then
! --------- Probleme: dans nmcere, on a touhjours ldccv(1:2) .gt. 0 !
            if (ldccv(1) .eq. 0 .and. ldccv(2) .eq. 0) then
                miresi = min(f(1), f(2))/max(f(1), f(2))
                if (miresi .le. contra) goto 600
            end if
        else
            goto 500
        end if
    end if
!
! --- STRATEGIE MIXTE: LES DEUX CONTRASTES SONT INSUFFISANTS
! --- ON REVIENT SUR NORM_INCR_DEPL ET ON TESTE LES CYCLES
!
    if (mixte) then
        f(1) = fnid(1)
        f(2) = fnid(2)
        ldccv(1) = 0
        ldccv(2) = 0
        switch = nmrcyc(sddisc, iterat, precyc)
        goto 600
    end if
!
500 continue
!
! --- PERMUTATION PAR EVENT DRIVEN (HORS STRATEGIE 'MIXTE')
!
    call utdidt('L', sddisc, 'ECHE', 'CHOIX_SOLU_PILO', &
                valk_=choix)
    ASSERT(choix .eq. 'NATUREL' .or. choix .eq. 'AUTRE')
    if (choix .eq. 'AUTRE' .or. swloun) then
        switch = .true.
        txt = 'NATUREL'
        if (choix .eq. 'AUTRE') then
            call utdidt('E', sddisc, 'ECHE', 'CHOIX_SOLU_PILO', valk_=txt)
        end if
    end if
!
600 continue
!
! --- RETOUR DE LA SELECTION AVEC EVENTUELLEMENT INTERVERSION
!
    sel = 2
    if ((f(1) .le. f(2) .and. .not. switch) .or. (f(1) .gt. f(2) .and. switch)) then
        sel = 1
    end if
    etaf = eta(sel)
    pilcvg = licite(sel)
    ldccvg = ldccv(sel)
    criter = f(sel)
!
end subroutine
