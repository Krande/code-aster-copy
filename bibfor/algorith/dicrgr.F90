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

subroutine dicrgr(DD, icodma, varim, klv, varip, fono, sip)
!
! --------------------------------------------------------------------------------------------------
!
!  COMPORTEMENT DIS_GRICRA : LIAISON GRILLE-CRAYON
!
!       ÉLASTIQUE PARTOUT SAUF SUIVANT Y LOCAL : FROTTEMENT DE COULOMB
!
!       ÉLÉMENTS MECA_DIS_TR_L
!
! --------------------------------------------------------------------------------------------------
!
! IN :
!   DD      : Cf appel
!   icodma  : Adresse du matériau code
!   varim   : Variables internes à l'instant précédent
!
! OUT :
!   klv     : Matrice tangente dans le repère local
!   varip   : Variables internes réactualisées
!   fono    : Forces nodales
!   sip     : Efforts internes
!
! --------------------------------------------------------------------------------------------------
!
! Variables Internes
!   V1 : déplacement plastique cumulé (direction axiale)
!   V2 : indicateur de contact/frottement (1 si glissement, 0 si non glissement)
!   V3 : indicateur de décollement en rotation
!   V4 : angle plastique (glissement)
!   V5 : angle plastique cumulé
!   V6 : angle plastique (glissement) nr 2
!   V7 : effort dans la direction du glissement (Y)
!   V8 : mémorisation de l’historique d’irradiation (fluence)
! --------------------------------------------------------------------------------------------------
!
!   Les relations sont définies dans le repère local.
!
    use te0047_type
    implicit none
#include "jeveux.h"
#include "asterfort/diklvraid.h"
#include "asterfort/moytem.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "asterfort/assert.h"
!
    type(te0047_dscr), intent(in) :: DD
    integer(kind=8) :: icodma
    real(kind=8) :: varim(*)
    real(kind=8) :: klv(78), varip(*), fono(*), sip(*)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: option
    integer(kind=8) :: neq
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    character(len=8) :: nompar(2)
!
    integer(kind=8) :: npg, nno, nbpar, iretp, iretm, iret2
    real(kind=8) :: valpar(2), irram, irrap
    real(kind=8) :: fser, h1
    real(kind=8) :: dux, duy, duz, dph, dth
    real(kind=8) :: uxm, uym, uzm, phm, thm
    real(kind=8) :: kxx, kyy, kzz, kpp, ktt, kkk
    real(kind=8) :: temp, prec
    real(kind=8) :: kt_ax, coul_ax, et_ax, kn_ax, et_rot
    real(kind=8) :: ang1, ang2, ang3, pen1, pen2, pen3, pen4
    real(kind=8) :: ktrig, dp, pm, rtm, rtp, reac_y
    real(kind=8) :: ecrotra, forceseuil
    real(kind=8) :: phipl, ppm, fphi, mophi, motheta, dkh, khm, fphi2
    real(kind=8) :: phitan, thetan
    real(kind=8) :: phiseuil, kphi, thetac, ktheta, kphi2, phiseuil2, phipl2
    real(kind=8) :: dpp, dphipl, dphipl2, ktheta2
    real(kind=8) :: tempp, tempm
    real(kind=8) :: fl(12), raideldc(6)
!
    integer(kind=8)             :: codre1(6), codre2(1), codre3(1), codre4(7), codre5(7)
    real(kind=8)        :: valp_fix(6), valp_tr(1), valp_rot(7)
    character(len=8)    :: nomre1(6), nomre2, nomre3, nomre4(7), nomre5(7)
!
    data nomre1/'KN_AX', 'KT_AX', 'ET_AX', 'ET_ROT', 'COUL_AX', 'KP'/
    data nomre2/'F_SER'/
    data nomre3/'F_SER_FO'/
    data nomre4/'ANG1', 'ANG2', 'PEN1', 'PEN2', 'PEN3', 'ANG3', 'PEN4'/
    data nomre5/'ANG1_FO', 'ANG2_FO', 'PEN1_FO', 'PEN2_FO', 'PEN3_FO', 'ANG3_FO', 'PEN4_FO'/
!
! --------------------------------------------------------------------------------------------------
!
    option = DD%option
    neq = DD%nno*DD%nc
!
!   Type d’élément : SEG2
    npg = 2
!
    prec = 1.0d-06
    fl = 0.0d0
!
!   Incrément d'irradiation sur le 1er pg (Sur les discrets il n'y en a qu'un)
    call rcvarc(' ', 'IRRA', '-', 'RIGI', 1, 1, irram, iret2)
    if (iret2 .gt. 0) irram = 0.d0
    call rcvarc(' ', 'IRRA', '+', 'RIGI', 1, 1, irrap, iret2)
    if (iret2 .gt. 0) irrap = 0.d0
!   L'irradiation ne peut que croître
    ASSERT((irrap-irram) .gt. -prec)
    irrap = irrap-irram+varim(8)
!
!   Récupération des données matériau
!
!   'KN_AX', 'KT_AX', 'ET_AX', 'ET_ROT', 'COUL_AX', 'KP'
!
    call rcvalb(fami, 1, 1, '+', icodma, ' ', 'DIS_GRICRA', 0, ' ', [0.d0], &
                6, nomre1, valp_fix, codre1, 0)
!
    ASSERT(ALL(codre1 .eq. 0))
!
!   F_SER
!   ou F_SER_FO
!
    call rcvalb(fami, 1, 1, '+', icodma, ' ', 'DIS_GRICRA', 0, ' ', [0.d0], &
                1, nomre2, valp_tr, codre2, 0)

    if (.not. ALL(codre2 .eq. 0)) then
        call moytem(fami, npg, 1, '+', tempp, iretp)
        call moytem(fami, npg, 1, '-', tempm, iretm)
        if ((iretp+iretm) .ge. 1) then
            call utmess('F', 'COMPOR5_40', sk='DIS_GRICRA')
        end if
        temp = (tempp+tempm)/2.d0
        !
        nbpar = 2
        nompar(1) = 'TEMP'
        valpar(1) = temp
        nompar(2) = 'IRRA'
        valpar(2) = irrap
        !
        call rcvalb(fami, 1, 1, '+', icodma, ' ', 'DIS_GRICRA', nbpar, nompar, valpar, &
                    1, nomre3, valp_tr, codre3, 0)
        ASSERT(ALL(codre3 .eq. 0))
    end if
!
!   'ANG1', 'ANG2', 'PEN1', 'PEN2', 'PEN3', 'ANG3', 'PEN4'
!   ou 'ANG1_FO', 'ANG2_FO', 'PEN1_FO', 'PEN2_FO', 'PEN3_FO', 'ANG3_FO', 'PEN4_FO'
!
    call rcvalb(fami, 1, 1, '+', icodma, ' ', 'DIS_GRICRA', 0, ' ', [0.d0], &
                7, nomre4, valp_rot, codre4, 0)
    if (.not. ALL(codre4 .eq. 0)) then
        call moytem(fami, npg, 1, '+', tempp, iretp)
        call moytem(fami, npg, 1, '-', tempm, iretm)
        if ((iretp+iretm) .ge. 1) then
            call utmess('F', 'COMPOR5_40', sk='DIS_GRICRA')
        end if
        temp = (tempp+tempm)/2.d0
        nbpar = 2
        nompar(1) = 'TEMP'
        valpar(1) = temp
        nompar(2) = 'IRRA'
        valpar(2) = irrap
        call rcvalb(fami, 1, 1, '+', icodma, ' ', 'DIS_GRICRA', nbpar, nompar, valpar, &
                    7, nomre5, valp_rot, codre5, 0)
        ASSERT(ALL(codre5 .eq. 0))
    end if
!
!   Extraction
!
    kn_ax = valp_fix(1)
    kt_ax = valp_fix(2)/4.d0
    et_ax = valp_fix(3)
    et_rot = valp_fix(4)
    coul_ax = valp_fix(5)
!
    fser = valp_tr(1)/4.d0
!
    ang1 = valp_rot(1)
    ang2 = valp_rot(2)
    pen1 = valp_rot(3)
    pen2 = valp_rot(4)
    pen3 = valp_rot(5)
    ang3 = valp_rot(6)
    pen4 = valp_rot(7)
!
    phiseuil = ang1
    thetac = ang2
    ktheta = (pen2-pen3+pen4)/2.d0
    kphi = (pen1-pen2)/2.d0
    ktheta2 = pen4/2.d0
    phiseuil2 = ang3
    kphi2 = (pen3-pen4)/2.d0
!
!   Variables internes de contact au temps moins
    h1 = 1.d0-varim(3)
!
!   Calcul de l'évolution des variables internes et des forces pour FULL_MECA et RAPH_MECA
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')) then
        ! extension de l'élément au pas de temps précédent
        uxm = DD%ulm(1+DD%nc)-DD%ulm(1)
        uym = DD%ulm(2+DD%nc)-DD%ulm(2)
        uzm = DD%ulm(3+DD%nc)-DD%ulm(3)
        phm = DD%ulm(4+DD%nc)-DD%ulm(4)
        khm = DD%ulm(5+DD%nc)-DD%ulm(5)
        thm = DD%ulm(6+DD%nc)-DD%ulm(6)
        ! variation d'extension
        dux = DD%dul(1+DD%nc)-DD%dul(1)
        duy = DD%dul(2+DD%nc)-DD%dul(2)
        duz = DD%dul(3+DD%nc)-DD%dul(3)
        dph = DD%dul(4+DD%nc)-DD%dul(4)
        dkh = DD%dul(5+DD%nc)-DD%dul(5)
        dth = DD%dul(6+DD%nc)-DD%dul(6)
        !
        ! Calcul des forces en translation
        !   on garde -FN0 dans la direction du discret
        !   possibilité de frottement dans la direction 2 (verticalement)
        !   force nulle dans la 3e direction

        ! Réaction pour le calcul du glissement à l'instant précédent
        rtm = varim(7)
        ! Plasticité à l'instant précédent
        pm = varim(1)
        ! Réaction actualisée
        reac_y = rtm+kt_ax*duy
        ! Seuil glissement
        forceseuil = coul_ax*fser
        !
        ecrotra = kt_ax*et_ax/(1-et_ax)
        !
        if (abs(reac_y) .le. (forceseuil+ecrotra*pm)) then
            ! Absence de glissement
            rtp = reac_y
            ktrig = kt_ax
            varip(1) = pm
            varip(2) = 0.d0
        else
            ! Glissement
            dp = (abs(reac_y)-(forceseuil+ecrotra*pm))/(ecrotra+kt_ax)
            rtp = (forceseuil+ecrotra*(pm+dp))*reac_y/abs(reac_y)
            ktrig = kt_ax*et_ax
            varip(1) = pm+dp
            varip(2) = 1.d0
        end if
        !
        varip(7) = rtp
        !
        ! Mise en // de KP sur (kx,ky,kz)
        sip(1) = -fser+kn_ax*(uxm+dux)
        sip(2) = rtp
        sip(3) = 0.d0
        sip(1+DD%nc) = sip(1)
        sip(2+DD%nc) = sip(2)
        sip(3+DD%nc) = sip(3)
        !
        fl(1) = -sip(1)
        fl(2) = -sip(2)
        fl(3) = -sip(3)
        fl(1+DD%nc) = -fl(1)
        fl(2+DD%nc) = -fl(2)
        fl(3+DD%nc) = -fl(3)
        !
        ! Calcul des forces en rotation
        !
        ! Calcul des forces et évolution des variables internes
        ! sur chacun des sous éléments du système de liaison
        !
        !  angle theta (DRZ): possibilité de décollement cela équivaut à une loi bilinéaire
        if ((thetac-thm-dth) .lt. 0.d0) then
            h1 = 0.d0
            motheta = ktheta*thetac+ktheta2*(thm+dth-thetac)
        else if ((thetac+thm+dth) .lt. 0.d0) then
            h1 = 0.d0
            motheta = -ktheta*thetac+ktheta2*(thm+dth+thetac)
        else
            h1 = 1.d0
            motheta = ktheta*(thm+dth)
        end if
        !
        thetan = h1*ktheta+(1.d0-h1)*ktheta2
        varip(3) = 1.d0-h1
        !
        sip(6) = motheta
        sip(6+DD%nc) = sip(6)
        fl(6) = -sip(6)
        fl(6+DD%nc) = -fl(6)
        !
        ! angle phi (DRX): loi de frottement
        phipl = varim(4)
        ! Plasticité à l'instant précédent
        ppm = varim(5)
        ! Angle PHI actualisé
        fphi = phm+dph-phipl
        !
        if (abs(fphi) .lt. (phiseuil+et_rot*ppm)) then
            mophi = kphi*fphi
            phitan = kphi
            varip(4) = phipl
            varip(5) = ppm
        else
            dpp = (abs(fphi)-(phiseuil+et_rot*ppm))/(1.d0+et_rot*(fphi)/abs(fphi))
            dphipl = dpp*(fphi)/abs(fphi)

            mophi = kphi*(phm+dph-phipl-dphipl)
            phitan = kphi*et_rot
            varip(4) = phipl+dphipl
            varip(5) = ppm+dpp
        end if
        !
        phipl2 = varim(6)
        fphi2 = phm+dph-phipl2
        if (abs(fphi2) .lt. phiseuil2) then
            mophi = mophi+kphi2*fphi2
            phitan = phitan+kphi2
            varip(6) = phipl2
        else
            dphipl2 = (abs(fphi2)-phiseuil2)*(fphi2)/abs(fphi2)
            mophi = mophi+kphi2*(phm+dph-phipl2-dphipl2)
            varip(6) = phipl2+dphipl2
        end if
        !
        sip(4) = mophi
        sip(4+DD%nc) = sip(4)
        fl(4) = -sip(4)
        fl(4+DD%nc) = -fl(4)
        !
        ! Pour le dernier angle, on met quand même une rigidité, même si
        !   les conditions limites doivent imposer que ça ne tourne pas
        sip(5) = kphi*(khm+dkh)
        sip(5+DD%nc) = sip(5)
        fl(5) = -sip(5)
        fl(5+DD%nc) = -fl(5)
        !
        varip(8) = irrap
        !
    end if
!
!   Matrice tangente pour FULL_MECA(_ELAS) et RIGI_MECA(_ELAS)
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RIGI_MECA')) then
        if ((option(1:9) .eq. 'RIGI_MECA') .or. (option(10:14) .eq. '_ELAS')) then
            kxx = kn_ax
            kyy = kt_ax
            kzz = 0.d0
            kpp = kphi+kphi2
            kkk = kphi+kphi2
            ktt = ktheta
        else
            kxx = kn_ax
            kyy = ktrig
            kzz = 0.d0
            kpp = phitan
            kkk = kphi
            ktt = thetan
        end if
        ! Mise en // de KP sur (kx,ky,kz)
        raideldc(1) = kxx
        raideldc(2) = kyy
        raideldc(3) = kzz
        raideldc(4) = kpp
        raideldc(5) = kkk
        raideldc(6) = ktt
        !
        klv(1:78) = 0.0d0
        call diklvraid(DD%nomte, klv, raideldc)
    end if
!
!   Calcul des forces nodales dans le repère global
    if (DD%lVect) then
        nno = 2
        call utpvlg(nno, DD%nc, DD%pgl, fl, fono)
    end if
!
end subroutine
