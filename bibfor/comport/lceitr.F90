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

subroutine lceitr(fami, kpg, ksp, mat, option, &
                  mu, su, de, ddedt, vim, &
                  vip, r)
!
! person_in_charge: jerome.laverne at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    character(len=16) :: option
    integer(kind=8) :: mat, kpg, ksp
    real(kind=8) :: mu(3), su(3), de(6), ddedt(6, 6), vim(*), vip(*), r
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
!            LOI DE COMPORTEMENT COHESIVE CZM_TRA_MIX
!            POUR LES ELEMENTS D'INTERFACE 2D ET 3D.
!
! IN : FAMI,KPG,KSP,MAT,OPTION
!      MU  : LAGRANGE
!      SU  : SAUT DE U
!      VIM : VARIABLES INTERNES
!
! OUT : DE    : DELTA, SOLUTION DE LA MINIMISATION
!       DDEDT : D(DELTA)/DT
!       VIP   : VARIABLES INTERNES MISES A JOUR
!       R     : PENALISATION DU LAGRANGE
!-----------------------------------------------------------------------
    integer(kind=8) :: nbpar
    parameter(nbpar=6)
    aster_logical :: resi, rigi, elas
    integer(kind=8) :: regime
    real(kind=8) :: sc, gc, c, h, ka, sk, st, val(nbpar), tmp, kap, skp, gap
    real(kind=8) :: dn, tn, t(3), ddndtn, dele, delp, delc, coee, coep
    integer(kind=8) :: cod(nbpar)
    character(len=16) :: nom(nbpar)
    character(len=1) :: poum
!
    data nom/'GC', 'SIGM_C', 'COEF_EXTR', 'COEF_PLAS',&
     &          'PENA_LAGR', 'RIGI_GLIS'/
!-----------------------------------------------------------------------
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
    resi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH'
    rigi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RIGI'
    elas = option(11:14) .eq. 'ELAS'
!
! RECUPERATION DES PARAMETRES PHYSIQUES
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mat, &
                ' ', 'RUPT_DUCT', 0, ' ', [0.d0], &
                nbpar, nom, val, cod, 2)
!
    gc = val(1)
    sc = val(2)
    coee = val(3)
    coep = val(4)
!
    if (coee .gt. coep) then
        call utmess('F', 'COMPOR1_67')
    end if
!
! EVALUATION DES SAUTS CRITIQUE, EXTRINSEQUE ET PLASTIQUE
    delc = 2*gc/(sc*(1-coee+coep))
    dele = coee*delc
    delp = coep*delc
!
    h = sc/(delc-delp)
    r = h*val(5)
    c = h*val(6)
!
! SEUIL COURANT ET CONTRAINTE CRITIQUE COURANTE
    ka = max(dele, vim(1))
    sk = max(0.d0, min(sc, sc*(ka-delc)/(delp-delc)))
!
!    FORCES COHESIVES AUGMENTEES
    t(1) = mu(1)+r*su(1)
    t(2) = mu(2)+r*su(2)
    t(3) = mu(3)+r*su(3)
    tn = t(1)
!
! -- CALCUL DE DELTA
! ------------------
!
!     ON VA TESTER ST : VALEUR DE LA DROITE EN KA
    st = -r*ka+tn
!
!    SI RIGI_MECA_*
    if (.not. resi) then
        regime = nint(vim(2))
        goto 5000
    end if
!
!    CONTACT
    if (st .le. -r*dele*sk/sc) then
        regime = -1
        dn = ka-dele*sk/sc
!
!    DECHARGE
    else if ((-r*dele*sk/sc .lt. st) .and. (st .le. sk)) then
        regime = 0
        dn = (dele*(tn-sk)+sc*ka)/(r*dele+sc)
!
!    PLATEAU
    elseif ((sk .lt. st) .and. (st .le. (max(0.d0, r*(delp-ka))+sk))) &
        then
        regime = 3
        dn = (tn-sc)/r
!
!    ENDOMMAGEMENT
    elseif (((max(0.d0, r*(delp-ka))+sk) .lt. st) .and. (st .le. r*( &
                                                         delc-ka))) then
        regime = 1
        dn = (tn-h*delc)/(r-h)
!
!    RUPTURE (SURFACE LIBRE)
    else
        regime = 2
        dn = tn/r
    end if
!
    call r8inir(6, 0.d0, de, 1)
!
!    COMPOSANTE DE L'OUVERTURE :
    de(1) = dn
!    COMPOSANTES DE CISAILLEMENT : ELASTIQUE
    de(2) = t(2)/(c+r)
    de(3) = t(3)/(c+r)
!
!
! -- ACTUALISATION DES VARIABLES INTERNES
! ---------------------------------------
!
!   V1 :  PLUS GRANDE NORME DU SAUT (SEUIL EN SAUT)
!   V2 :  REGIME DE LA LOI
!        -1 : CONTACT
!         0 : ADHERENCE OU DECHARGE
!         3 : PLATEAU
!         1 : ENDOMMAGEMENT
!         2 : RUPTURE (SURFACE LIBRE)
!   V3 :  INDICATEUR D'ENDOMMAGEMENT
!         0 : SAIN
!         1 : ENDOMMAGE
!         2 : CASSE
!   V4 :  POURCENTAGE D'ENERGIE DISSIPEE
!   V5 :  VALEUR DE L'ENERGIE DISSIPEE (V4*GC)
!   V6 :  ENERGIE RESIDUELLE COURANTE
!        (NULLE POUR CE TYPE D'IRREVERSIBILITE)
!   V7 A V9 : VALEURS DE DELTA
!
    kap = min(max(ka, dn), delc)
    skp = max(0.d0, min(sc, sc*(kap-delc)/(delp-delc)))
    tmp = kap+dele*(1.d0-skp/sc)
    gap = sc*(kap-dele*skp/sc)-(tmp-delp)*(sc-skp)/2.d0
    gap = gap/gc
!
    vip(1) = kap
    vip(2) = regime
!
    if (kap .eq. dele) then
        vip(3) = 0.d0
    else if (kap .eq. delc) then
        vip(3) = 2.d0
    else
        vip(3) = 1.d0
    end if
!
    vip(4) = gap
    vip(5) = gc*vip(4)
    vip(6) = 0.d0
    vip(7) = de(1)
    vip(8) = de(2)
    vip(9) = de(3)
!
!
! -- MATRICE TANGENTE
!--------------------
!
5000 continue
    if (.not. rigi) goto 9999
!
!    AJUSTEMENT POUR PRENDRE EN COMPTE *_MECA_ELAS
    if (elas) then
        if (regime .eq. 1) regime = 0
    end if
!
    call r8inir(36, 0.d0, ddedt, 1)
!
    ddedt(2, 2) = 1/(c+r)
    ddedt(3, 3) = 1/(c+r)
!
    if (regime .eq. 0) then
        ddndtn = dele/(dele*r+sc)
    else if (regime .eq. 3) then
        ddndtn = 1/r
    else if (regime .eq. 1) then
        ddndtn = 1/(r-h)
    else if (regime .eq. 2) then
        ddndtn = 1/r
    else if (regime .eq. -1) then
        ddndtn = 0
    end if
    ddedt(1, 1) = ddndtn
!
9999 continue
!
end subroutine
