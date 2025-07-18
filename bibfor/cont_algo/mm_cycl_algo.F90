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
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
! aslint: disable=W1504
!
subroutine mm_cycl_algo(ds_contact, l_frot_zone, &
                        l_glis_init, type_adap, zone_index, i_cont_poin, &
                        indi_cont_eval, indi_frot_eval, dist_cont_curr, &
                        pres_cont_curr, dist_frot_curr, pres_frot_curr, v_sdcont_cychis, &
                        v_sdcont_cyccoe, v_sdcont_cyceta, indi_cont_curr, indi_frot_curr, &
                        ctcsta, mmcvca, l_pena_frot, l_pena_cont, vale_pene, glis_maxi)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/bussetta_algorithm.h"
#include "asterfort/cfdisi.h"
#include "asterfort/mm_cycl_detect.h"
#include "asterfort/mm_cycl_trait.h"
#include "asterfort/mmstac.h"
#include "asterfort/search_opt_coef.h"
#include "Contact_type.h"
!
    type(NL_DS_Contact), intent(inout) :: ds_contact
    aster_logical, intent(in) :: l_frot_zone
    aster_logical, intent(in) :: l_glis_init
    aster_logical, intent(in) :: l_pena_frot
    aster_logical, intent(in) :: l_pena_cont
    integer(kind=8), intent(inout) :: type_adap
    integer(kind=8), intent(in) :: i_cont_poin
    integer(kind=8), intent(in) :: zone_index
    integer(kind=8), intent(inout) :: indi_cont_eval
    integer(kind=8), intent(inout) :: indi_frot_eval
    real(kind=8), intent(in) :: vale_pene
    real(kind=8), intent(in) :: glis_maxi
    real(kind=8), intent(inout) :: dist_cont_curr
    real(kind=8), intent(inout) :: pres_cont_curr
    real(kind=8), intent(inout) :: dist_frot_curr(3)
    real(kind=8), intent(in) :: pres_frot_curr(3)
    real(kind=8), pointer :: v_sdcont_cychis(:)
    real(kind=8), pointer :: v_sdcont_cyccoe(:)
    integer(kind=8), pointer :: v_sdcont_cyceta(:)
    integer(kind=8), intent(out) :: indi_cont_curr
    integer(kind=8), intent(out) :: indi_frot_curr
    integer(kind=8), intent(out) :: ctcsta
    aster_logical, intent(out) :: mmcvca
!
! ---------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE CONTINUE - CONTRAINTES ACTIVES)
! Cette routine va adapter le calcul des matrices de contact-frottement
! suivant les cas de cyclage.
! Il permet également de faire des heuristiques particulières (cas du flip-flop)
! pour aider à la convergence
! sur les statuts de contact. Cete routine permet également
! de faire une recherche optimale du coefficient de
! C'est également dans cette routine qu'on s'assure de la convergence
! de la boucle sur les statuts.
!
! TRAITEMENT DES DIFFERENTS CAS
!
! -----------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! In  l_frot_zone      : .true. if friction on zone
! In  l_coef_adap      : .true. if automatic lagrangian adaptation
! In  l_glis_init      : .true. if bilateral contact for first step
! In  i_cont_poin      : contact point index
! In  indi_cont_eval   : evaluation of new contact status
! In  indi_frot_eval   : evaluation of new friction status
! In  dist_cont_curr   : current contact gap
! In  pres_cont_curr   : current contact pressure
! In  dist_frot_curr   : current friction distance
! In  pres_frot_curr   : current friction pressure
! In  v_sdcont_cychis  : pointer to cycling history
! In  v_sdcont_cyccoe  : pointer to coefficient history
! Out indi_cont_curr   : current contact status
! Out indi_frot_curr   : current friction status
! Out mmcvca           : .true. if contact loop converged
! Out ctcsta           : number of contact points has changed their status
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: coef_cont_prev, coef_frot_prev
    real(kind=8) :: coef_cont_curr, coef_frot_curr
    aster_logical:: coef_found, treatment
    aster_logical:: l_coef_adap, mmcvca_frot
    integer(kind=8)      ::  mode_cycl
    real(kind=8) :: pres_frot_prev(3), pres_cont_prev
    real(kind=8) :: dist_frot_prev(3), dist_cont_prev
    integer(kind=8) :: indi_cont_prev, indi_frot_prev, indi(2), i_reso_cont
    real(kind=8) :: coef_frot_mini, coef_frot_maxi
    real(kind=8) :: alpha_cont_matr, alpha_cont_vect
    real(kind=8) :: alpha_frot_matr, alpha_frot_vect
    real(kind=8) :: coef_opt, pres_cont(2), dist_cont(2)
    real(kind=8) :: coef_bussetta, dist_max, coe1
    real(kind=8) ::  coef_tmp
    real(kind=8) :: bound_coef(2)
!
    bound_coef(1:2) = [1.d-8, 1.d8]
    mmcvca_frot = ASTER_TRUE
    coef_cont_prev = 0.0
    coef_frot_prev = 0.0
    coef_cont_curr = 0.0
    coef_frot_curr = 0.0
    coef_found = ASTER_FALSE
    treatment = ASTER_TRUE
    l_coef_adap = ASTER_FALSE
    mmcvca_frot = ASTER_TRUE
    mode_cycl = 0
    pres_frot_prev(:) = 0.0
    pres_cont_prev = 0.0
    dist_frot_prev(:) = 0.0
    dist_cont_prev = 0.0
    indi_cont_prev = 0
    indi_frot_prev = 0
    indi(:) = 0
    coef_frot_mini = 0.0
    coef_frot_maxi = 0.0
    alpha_cont_matr = 0.0
    alpha_cont_vect = 0.0
    alpha_frot_matr = 0.0
    alpha_frot_vect = 0.0
    coef_opt = 0.0
    pres_cont(:) = 0.0
    dist_cont(:) = 0.0
    coef_bussetta = 0.0
    dist_max = 0.0
    coe1 = 0.d0

    mmcvca = ASTER_FALSE
!
    i_reso_cont = cfdisi(ds_contact%sdcont_defi, 'ALGO_RESO_CONT')
!
! ---------------------------------------------------------------------
! - Initializations
! ---------------------------------------------------------------------------------
!
! type_adap vient de cazocc : v_sdcont_paraci(20)
    l_coef_adap = ((type_adap .eq. 1) .or. (type_adap .eq. 2) .or. &
                   (type_adap .eq. 5) .or. (type_adap .eq. 6))

!----------------------TRAITEMENT CYCLAGE -------------------------
! Quand est-ce qu'on fait de l'adaptation des matrices de contact-frottement
! MATRCF[Iteration_K]=ALPHA*MATRCF[Iteration_K-1]+(1-ALPHA)*MATRCF[Iteration_K]
! type_adap vient de cazocc : v_sdcont_paraci(20)
! CAS 1 : adaptation .eq. 'CYCLAGE' + Quelque soit ALGO_CONT/ALGO_FROT, type_adap=4
! CAS 2 : adaptation .eq. 'TOUT' + NEWT_FROT , type_adap=5
! CAS 3 : adaptation .eq. 'TOUT' + NEWT_FROT , ALGO_CONT = PENALISATION, type_adap=6
! CAS 4 : adaptation .eq. 'TOUT' + POINT_FIXE_FROT OU PAS DE FROT +
! ALGO_CONT = PENALISATION, type_adap=7
! CAS 5 : tous les autres cas du moment ou adaptation .eq. 'TOUT' actif, type_adap=4
    treatment = ((type_adap .eq. 4) .or. (type_adap .eq. 5) .or. &
                 (type_adap .eq. 6) .or. &
                 (type_adap .eq. 7) .or. (type_adap .eq. 11))

    if (nint(v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24)) .ne. &
        nint(v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+24))) then
!        Precaution :  "la maille maitre a changé : on ne fait pas de cyclage
        treatment = .false.
    end if

! - Previous informations
    indi_cont_prev = nint(v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+1))
    coef_cont_prev = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+2)
    pres_cont_prev = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+3)
    dist_cont_prev = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+4)

    indi_frot_prev = nint(v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+5))
    coef_frot_prev = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+6)
    pres_frot_prev(1) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+7)
    pres_frot_prev(2) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+8)
    pres_frot_prev(3) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+9)
    dist_frot_prev(1) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+10)
    dist_frot_prev(2) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+11)
    dist_frot_prev(3) = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+12)

! - Current max/min ratio
    coef_frot_mini = v_sdcont_cyccoe(6*(zone_index-1)+5)
    coef_frot_maxi = v_sdcont_cyccoe(6*(zone_index-1)+6)

! ---------------------------------------------------------------------------
! - Cycling detection
! -----------------------------------------------------------------
    call mm_cycl_detect(ds_contact, l_frot_zone, i_cont_poin, &
                        coef_cont_prev, coef_frot_prev, pres_cont_prev, &
                        dist_cont_prev, pres_frot_curr, pres_frot_prev, &
                        indi_frot_prev, dist_frot_prev, indi_cont_eval, &
                        indi_frot_eval, indi_cont_prev, dist_cont_curr, pres_cont_curr, &
                        dist_frot_curr, alpha_cont_matr, alpha_cont_vect, &
                        alpha_frot_matr, alpha_frot_vect)

!
! ---------------------------------------------------------------
! - Cycling treatment: automatic adaptation of augmented lagrangian ratio
! ----------------------------------------------------------------------
!

! L'ordre de traitement est le suivant :
! Step 1. Traitement de l'adaptation des coefficients
! Step 2. Traitement spécial du mot-clef GLIISIERE=OUI, label 236
! Step 3. Traitement du FLIP-FLOP : POINT_FIXE SUR LE CONTACT
! Step 4. Traitement du CYCLAGE : NEWTON SUR LE CONTACT
! Step 5 : adaptation du coefficient penalisation
! ----------------------------------------------------------------
!
! Step 1.  Traitement de l'adaptation des coefficients
!
    if (l_coef_adap) then
        call mm_cycl_trait(ds_contact, i_cont_poin, coef_cont_prev, coef_frot_prev, &
                           pres_frot_prev, dist_frot_prev, pres_frot_curr, &
                           dist_frot_curr, &
                           indi_cont_eval, indi_frot_eval, indi_cont_curr, &
                           coef_cont_curr, &
                           indi_frot_curr, coef_frot_curr)
    else
        coef_cont_curr = coef_cont_prev
        coef_frot_curr = coef_frot_prev
        indi_cont_curr = indi_cont_eval
        indi_frot_curr = indi_frot_eval
    end if

! - Saving max/min regularization coefficients
    if (coef_frot_curr .ge. coef_frot_maxi) coef_frot_maxi = coef_frot_curr
    if (coef_frot_curr .le. coef_frot_mini) coef_frot_mini = coef_frot_curr
    v_sdcont_cyccoe(6*(zone_index-1)+5) = coef_frot_mini
    v_sdcont_cyccoe(6*(zone_index-1)+6) = coef_frot_maxi

!
! Step 2. Traitement spécial du mot-clef GLIISIERE=OUI, label 236
!
    if (l_glis_init) then
        indi_cont_curr = 1
        mmcvca = .true.
        goto 236
    end if

! - Save history for automatic cycling algorithm
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+1) = indi_cont_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = coef_cont_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+3) = pres_cont_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+4) = dist_cont_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+5) = indi_frot_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6) = coef_frot_curr
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+7) = pres_frot_curr(1)
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+8) = pres_frot_curr(2)
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+9) = pres_frot_curr(3)
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+10) = dist_frot_curr(1)
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+11) = dist_frot_curr(2)
    v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+12) = dist_frot_curr(3)

!
! Step 3. Traitement du FLIP-FLOP : POINT_FIXE SUR LE CONTACT
!
!     if (i_reso_cont .eq. 0) then
!         if (v_sdcont_cyceta(4*(i_cont_poin-1)+4) .eq. -10) then
!             if (ds_contact%resi_pressure .lt. 1.d-4*ds_contact%cont_pressure) &
!             mmcvca = .true.
!         endif
!     endif
!
! Step 4. Traitement du CYCLAGE : NEWTON SUR LE CONTACT
    if ((ds_contact%iteration_newton .ge. 3) .and. &
        v_sdcont_cyceta(4*(i_cont_poin-1)+1) .ne. 0 .and. &
        treatment) then
        if (v_sdcont_cyceta(4*(i_cont_poin-1)+1) .eq. 1) then
            ! critere contact ignoré : grezillant
            mmcvca = .true.
        else if (v_sdcont_cyceta(4*(i_cont_poin-1)+4) .eq. -10 .and. .false.) then
            if (ds_contact%resi_pressure .lt. 1.d-4*ds_contact%cont_pressure) then
!           Cas 1 : Le point fait du FLIP-FLOP, meme traitement que POINT_FIXE
                mmcvca = .true.
            end if
        else
            !ADAPTATION DE MATRICES, VECTEURS ET COEFF POUR LES TE :
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+57) = 1.0
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+59) = 0.9
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+56) = 0.9
        end if
!        recherche de coefficients est-ce qu'on peut trouver coef tel que :
!        Statut_Contact(Lag_prev-coef*Gap_prev) = Statut_Contact(Lag_curr-coef*Gap_curr)
        coef_found = .false.
        indi(1) = indi_cont_curr
        indi(2) = indi_cont_prev
        pres_cont(1) = pres_cont_curr
        pres_cont(2) = pres_cont_prev
        dist_cont(1) = dist_cont_curr
        dist_cont(2) = dist_cont_prev

        call search_opt_coef(bound_coef, &
                             indi, pres_cont, dist_cont, &
                             coef_opt, coef_found)
        if (coef_found) then
            if (i_reso_cont .eq. ALGO_NEWT) then
                indi_cont_curr = indi(1)
                indi_cont_prev = indi(2)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+1) = indi_cont_curr
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+1) = indi_cont_prev
!                mmcvca = .true.
!                goto 236
            end if
            dist_cont_curr = dist_cont(1)
            dist_cont_prev = dist_cont(2)
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = coef_opt
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+2) = coef_opt
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+3) = pres_cont_curr
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+3) = pres_cont_prev
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+4) = dist_cont_curr
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+4) = dist_cont_prev
        end if

    end if

! CYCLAGE FROTTEMENT    : ADHE_GLIS

    if ((ds_contact%iteration_newton .ge. 3) .and. &
        (v_sdcont_cyceta(4*(i_cont_poin-1)+2) .ge. 10) .and. treatment) then
        ! Cyclage ADHE_GLIS purement et simplement debranchee :
        ! Se referer a la version 14.1 mmalgo, pour rebranchement si besoin
        ! test sur ssnv128p : on passe de 10 iterations a 182 itertions
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+50) = 0.0d0
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+54) = 1.0
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+55) = 1.0
!         write (6,*) "Cyclage GLI_ADHE_GLIS"
    end if

! CYCLAGE FROTTEMENT: GLIAV_AR
    if ((ds_contact%iteration_newton .ge. 3) .and. &
        (v_sdcont_cyceta(4*(i_cont_poin-1)+3) .ge. 10) .and. treatment) then
        ! Cyclage GLIAV_AR  : comme conseille dans un cours de IPSI 2004
        ! Statut de frottement mis à adherence
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+5) = 1
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+50) = 1.0d0
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+54) = 0.9
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+55) = 0.9
!           write (6,*) "Cyclage GLIAV_AR"
    end if

!
!  Algorithm of Bussetta
!

    if ((type_adap .eq. 2) .or. (type_adap .eq. 3) .or. &
        (type_adap .eq. 6) .or. (type_adap .eq. 7) .or. &
        v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+51) .eq. 4.0) then

        coef_bussetta = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2)
        coef_tmp = v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2)

        if (indi_cont_curr .eq. 1 .and. l_pena_cont) then
            if (nint(vale_pene) .eq. -1) then
                ! Mode relatif
                dist_max = 1.d-2*ds_contact%arete_min
            else
                ! Mode absolu
                dist_max = vale_pene
            end if

!           mmcvca = mmcvca .and. (ctcsta .eq. 0)
!           call bussetta_algorithm(dist_cont_curr, dist_cont_prev,dist_max, coef_bussetta)
!           v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = max(coef_bussetta,&
!                                                   v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2))
            !endif

            ! Traitement de fortes interpenetrations
            !if (indi_cont_curr .eq. 1 .and. l_pena_cont) then

            if (dist_cont_curr .gt. dist_max) then
                coef_tmp = coef_tmp*(abs(dist_cont_curr)/dist_max)
                call bussetta_algorithm(dist_cont_curr, dist_cont_prev, dist_max, coef_bussetta)
                if (coef_bussetta .lt. coef_tmp) coef_bussetta = coef_tmp
                ! On approche de la fin des iterations de Newton mais penetration pas satisfait
                ! Le calcul du coefficient n'est pas satisfaisant on l'augmente
                if (nint(ds_contact%continue_pene) .eq. 1) coef_bussetta = coef_bussetta*3
                if (coef_bussetta .gt. ds_contact%max_coefficient) then
                    coef_bussetta = coef_bussetta*0.01
                    ! critere trop severe : risque de non convergence
                    ds_contact%continue_pene = 2.0
                end if
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = coef_bussetta
                ! critere trop lache
                if (dist_max .gt. ds_contact%arete_min) &
                    ds_contact%continue_pene = 1.0
            end if
        end if
!       ! cas ALGO_CONT=PENALISATION, ALGO_FROT=STANDARD
!       ! On fixe un statut adherent en cas de fortes interpenetration
!       if ((.not. l_pena_frot .and. l_pena_cont ).and. indi_cont_curr .eq. 1) then
!           if ((dist_cont_curr .gt. dist_max) .and. (indi_frot_curr .eq. 0.)&
!                .and. (norm2(dist_frot_curr) .lt. 0.01*dist_max)) then
!               v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+5)  = 1
!               v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6) =&
!               1.d1*ds_contact%arete_max/ds_contact%arete_min*coef_frot_curr
!           endif
!       endif
!        ! cas ALGO_CONT=PENALISATION, ALGO_FROT=PENALISATION
!        ! On fixe un statut glissement en cas de fortes interpenetration
!       if ((l_pena_frot .and. l_pena_cont ).and. indi_cont_curr .eq. 1) then
!           if ((dist_cont_curr .lt. dist_max) .and. (indi_frot_curr .eq. 0.)&
!                .and. (norm2(dist_frot_curr) .lt. 1.d-6*dist_max) .and. &
!                (indi_frot_curr .eq. indi_frot_prev )) then
!               v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+5)  = 0
!               v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6) =&
!               1.d-1*ds_contact%arete_min/ds_contact%arete_max*coef_frot_curr
!           endif
!       endif

    end if

!        ! cas ALGO_CONT=STANDARD, ALGO_FROT=PENALISATION
!        ! On fixe un statut glissement en cas de fortes interpenetration
    if ((l_pena_frot .and. (.not. l_pena_cont)) .and. indi_cont_curr .eq. 1) then
        if ((norm2(dist_frot_curr) .gt. 1.d-6*glis_maxi)) then
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+6) = &
                1.d1*ds_contact%arete_max/ds_contact%arete_min*1/glis_maxi
        end if

    end if

! - Convergence ?
!
    mmcvca = (indi_cont_prev .eq. indi_cont_curr) .or. mmcvca
    if (.not. mmcvca .and. treatment) then
!       On Bascule en mode penalise automatiquement jusqu'à convergence puis
!       On revient en standard quand le statut de contact se stabilise : ssnv128z
        mode_cycl = 1
        if (mode_cycl .eq. 1 .and. &
            ds_contact%iteration_newton .gt. ds_contact%it_cycl_maxi+3) then
            ! On fait la projection sur le cône négatif des valeurs admissibles
            if (dist_cont_curr .gt. 1.d-6) dist_cont_curr = 0.0
            if (pres_cont_curr .gt. 1.d-6) pres_cont_curr = -1.d-15
            if (dist_cont_prev .gt. 1.d-6) dist_cont_prev = 0.0
            if (pres_cont_prev .gt. 1.d-6) pres_cont_prev = -1.d-15
            if (i_reso_cont .eq. ALGO_NEWT) then
                call mmstac(dist_cont_curr, pres_cont_curr, coe1, indi_cont_curr)
                call mmstac(dist_cont_prev, pres_cont_prev, coe1, indi_cont_prev)
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+1) = indi_cont_curr
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+1) = indi_cont_prev
            end if
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+3) = pres_cont_curr
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+3) = pres_cont_prev
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+4) = dist_cont_curr
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+4) = dist_cont_prev
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+57) = 1.0
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+59) = 0.99
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+56) = 1.0
            v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+51) = 4.0
            v_sdcont_cyceta(4*(i_cont_poin-1)+1) = 10
            if (.not. l_pena_cont) then
                !coef_tmp =max(1.d0/ds_contact%arete_min,&
                !              v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2))
                coef_tmp = 1.d0/ds_contact%arete_min

                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+2) = coef_tmp
                v_sdcont_cychis(NB_DATA_CYCL*(i_cont_poin-1)+24+2) = coef_tmp
            end if
            mmcvca = indi_cont_prev .eq. indi_cont_curr
        end if
    end if

    if (.not. mmcvca) ctcsta = ctcsta+1
    mmcvca = (mmcvca .and. (ctcsta .eq. 0))
! Ne pas effracer ce comment &  .or. (ds_contact%resi_pressure .lt. 1.d-6*ds_contact%cont_pressure))
236 continue
end subroutine
