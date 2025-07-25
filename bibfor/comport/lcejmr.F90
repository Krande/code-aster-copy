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
! aslint: disable=W0413,W1306
!
subroutine lcejmr(BEHinteg, fami, kpg, ksp, ndim, &
                  mate, option, epsm, deps, &
                  sigma, dsidep, vim, vip, typmod, &
                  instam, instap)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8pi.h"
#include "asterfort/matinv.h"
#include "asterfort/pmavec.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: mate, ndim, kpg, ksp
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigma(6), dsidep(6, 6)
    real(kind=8) :: vim(*), vip(*), instam, instap
    character(len=8) :: typmod(*)
    character(len=16) :: option
    character(len=*) :: fami
!-----------------------------------------------------------------------
!     LOI DE COMPORTEMENT DES JOINTS DE BARRAGE : JOINT_MECA_RUPT
!     POUR LES ELEMENTS DE JOINT ET JOINT_HYME 2D ET 3D
!
! DEUX TYPES DE CALCUL SONT POSSIBLES:
!
! 1) MODELISATIONS *_JOINT_HYME, AVEC PARAMETRE HYDRO PRESENTS
!    CALCUL COUPLE HYDRO-MECANIQUE SUR UN MAILLAGE QUADRATIQUE
!
! 2) MODELISATIONS *_JOINT, PAS DE PARAMETRES HYDRO POSSIBLE
!    CALCUL MECA (AVEC ENVENTUELLEMENT PRES_CLAVAGE OU PRES_FLUIDE)
!    SUR UN MAILLAGE LINEAIRE OU QUADRATIQUE
!
! IN : EPSM - SAUT INSTANT MOINS ET GRAD PRESSION ET PRES FLUIDE SI HYME
! IN : DEPS - INC DE SAUT  ET INC GRAD PRESSION ET INC PRES FLUIDE SI HYME
! IN : MATE, OPTION, VIM, COOROT, INSTAM, INSTAP
! IN : SIGMO - SIGMA INSTANT MOINS ET FLUX HYDRO SI HYME
! OUT : SIGMA , DSIDEP , VIP
!-----------------------------------------------------------------------
    integer(kind=8) :: nbpa
    parameter(nbpa=12)
    integer(kind=8) :: cod(nbpa)
    integer(kind=8) :: i, n, diss, cass
    real(kind=8) :: sc, lc, lct, k0, val(nbpa), presfl, presg, prescl, sciage, tmp
    real(kind=8) :: gp(ndim-1), gploc(ndim), gpglo(ndim), fhloc(ndim), fhglo(ndim)
    real(kind=8) :: delta(ndim), ddelta(ndim), ka, r0, rc, alpha, beta, rk, ra, rt
    real(kind=8) :: rt0, r8bid
    real(kind=8) :: offset(ndim), doffset(ndim), inst, valpar(ndim+1), rhof, visf, amin
    real(kind=8) :: coorot(ndim+ndim*ndim), invrot(ndim, ndim), rigart
    character(len=8) :: nompar(ndim+1)
    character(len=16) :: nom(nbpa)
    character(len=1) :: poum
    aster_logical :: resi, rigi, elas, ifpahm, ifhyme
    blas_int :: b_incx, b_incy, b_n
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
!
! INDICATEUR AVEC/SANS HYDRO
    if (typmod(2) .eq. 'EJ_HYME') ifhyme = .true.
    if (typmod(2) .eq. 'ELEMJOIN') ifhyme = .false.
!
! SAUT DE DEPLACEMENT EN T- OU T+
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, delta, b_incy)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, deps, b_incx, ddelta, b_incy)
    if (resi) then
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, delta, &
                   b_incy)
    end if
!
! GRADIENT DE PRESSION ET PRESSION EN T- OU T+
    if (ifhyme) then
!
        do n = 1, ndim-1
            gp(n) = epsm(ndim+n)
            if (resi) gp(n) = gp(n)+deps(ndim+n)
        end do
!
        presg = epsm(2*ndim)
        if (resi) presg = presg+deps(2*ndim)
!
    end if
!
! INSTANT DE CALCUL T- OU T+
    inst = instam
    if (resi) inst = instap
!
! RECUPERATION DES PARAMETRES PHYSIQUES
!--------------------------------------
    nom(1) = 'K_N'
    nom(2) = 'SIGM_MAX'
    nom(3) = 'PENA_RUPTURE'
    nom(4) = 'PENA_CONTACT'
    nom(5) = 'ALPHA'
    nom(6) = 'K_T'
    nom(7) = 'PRES_FLUIDE'
    nom(8) = 'PRES_CLAVAGE'
    nom(9) = 'RHO_FLUIDE'
    nom(10) = 'VISC_FLUIDE'
    nom(11) = 'OUV_MIN'
    nom(12) = 'SCIAGE'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', 0, ' ', [0.d0], &
                5, nom, val, cod, 2)
!
! CONTRAINTE CRITIQUE SANS PENALISATION
    sc = val(2)*(1.d0+val(3))/val(3)
! LONGUEUR CRITIQUE AVANT LA RUPTURE COMPLETE DU JOINT
    lc = (1.d0+val(3))*val(2)/val(1)
! LONGUEUR AVANT L'ADOUCISSEMENT
    k0 = val(2)/val(1)
! PENTE NORMALE INITIAL
    r0 = val(1)
    beta = val(4)
! PARAMETRE QUI DEFINI LA LONGUEUR CRITIQUE TANGENTIELLE (0<=ALPHA<=2)
    alpha = val(5)
! LONGUEUR CRITIQUE TANGENTIELLE
! ALPHA=0: LCT=0; ALPHA=1: LCT=LC; ALPHA=2;LCT=INFTY
    if ((alpha .ge. 0.d0) .and. (alpha .lt. 2.d0)) then
        lct = lc*tan(alpha*r8pi()/4.d0)
    else
! PRESENTATION D'UNE INFINITE NUMERIQUE
        lct = (1.d0+lc)*1.d8
    end if
! PENTE TANGENTIELLE INITIAL (SI ELLE N'EST PAS DEFINI ALORS K_T=K_N)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', 0, ' ', [0.d0], &
                1, nom(6), val(6), cod(6), 0)
    if (cod(6) .eq. 0) then
        rt0 = val(6)
    else
        rt0 = r0
    end if
!
! DEFINITION DES PARAMETRES POUR LA RECUPERATION DES FONCTIONS
    coorot = 0.d0
    do i = 1, ndim
        coorot(i) = BEHinteg%behavESVA%behavESVAGeom%coorElga(kpg, i)
    end do
    do i = 1, ndim*ndim
        coorot(ndim+i) = BEHinteg%behavESVA%behavESVAOther%rotpg(i)
    end do
!
    nompar(1) = 'INST'
    nompar(2) = 'X'
    nompar(3) = 'Y'
    valpar(1) = inst
    valpar(2) = coorot(1)
    valpar(3) = coorot(2)
    if (ndim .eq. 3) then
        nompar(4) = 'Z'
        valpar(4) = coorot(3)
    end if
!
! RECUPERATION DE LA PRESS FLUIDE, CLAVAGE ET SCIAGE (MODELISATION MECA PURE)
!-----------------------------------------------------------------------
! RECUPERATION DE LA PRESS FLUIDE (FONCTION DE L'ESPACE ET DU TEMPS)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', ndim+1, nompar, valpar, &
                1, nom(7), val(7), cod(7), 0)
!
    if (cod(7) .eq. 0) then
        presfl = val(7)
    else
        presfl = 0.d0
    end if
!
! RECUPERATION DE LA PRESS CLAVAGE (FONCTION DE L'ESPACE ET DU TEMPS)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', ndim+1, nompar, valpar, &
                1, nom(8), val(8), cod(8), 0)
!
    if (cod(8) .eq. 0) then
        prescl = val(8)
    else
        prescl = -1.d0
    end if
!
! RECUPERATION DE LA TAILLE DE SCIE = SCIAGE (FONCTION DE L'ESPACE ET DU TEMPS)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', ndim+1, nompar, valpar, &
                1, nom(12), val(12), cod(12), 0)
!
    if (cod(12) .eq. 0) then
        sciage = val(12)
    else
        sciage = 0.d0
    end if
!
!
! RECUPERATION DE LA MASSE VOL ET DE LA VISCO (MODELISATION JOINT HM)
!--------------------------------------------------------------------
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', 0, ' ', [0.d0], &
                1, nom(9), val(9), cod(9), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', 0, ' ', [0.d0], &
                1, nom(10), val(10), cod(10), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_RUPT', 0, ' ', [0.d0], &
                1, nom(11), val(11), cod(11), 0)
!
    if (cod(9) .eq. 0) rhof = val(9)
    if (cod(10) .eq. 0) visf = val(10)
    if (cod(11) .eq. 0) amin = val(11)
!
! INDICATEUR SI LES PARAMETRES HYDRO SONT RENSEIGNES
    ifpahm = (cod(9) .eq. 0) .and. (cod(10) .eq. 0) .and. (cod(11) .eq. 0)
!
! VERIFICATION DE LA PRESENCE/ABSENCE DE PARAMETRES
! EN FONCTION DE LA MODELISATION MECA PUR OU HYDRO MECA
!
    if (ifhyme) then
!       POUR LE CALCUL HYDRO => PAS DE PRES_CLAVAGE
        if (cod(8) .eq. 0) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PAS DE SCIAGE
        if (cod(12) .eq. 0) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PAS DE PRES_FLUIDE
        if (cod(7) .eq. 0) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PRESENCE DE PARA_HM
        if (.not. ifpahm) then
            call utmess('F', 'ALGORITH17_15')
        end if
    else
!       POUR LE CALCUL MECA => PAS DE PARAMETRE HYDRO
        if (ifpahm) then
            call utmess('F', 'ALGORITH17_16')
        end if
    end if
!
! INITIALISATION DU SEUIL D'ENDOMMAGEMENT ACTUEL
    ka = max(k0, vim(1))
!
! DANS LE CAS DU CLAVAGE/SCIAGE
! INITIALISATION DU POINT D'EQUILIBRE POUR LA LDC (OFFSET)
!-----------------------
! CLAVAGE
! EPAISSEUR DE JOINT NE PEUT QUE AUGMENTER; DOFFSET(1) > 0
    doffset(1) = 0.d0
    if (prescl .ge. 0.d0) doffset(1) = max(0.d0, delta(1)+prescl/(beta*r0))
! SCIAGE
! LE JOINT EST ENDOMMAGE PAR LE SCIAGE
    if (sciage .gt. 0.d0) ka = lc
! L'EPASSEUR SCIEE EST DIMINUEE DE L'OUVERTURE INITALE DE JOINT
    sciage = sciage-max(0., epsm(1))
    if (sciage .gt. 0.d0) doffset(1) = doffset(1)-sciage
    offset(1) = vim(10)+doffset(1)
! OFFSET TANGENTIEL
    offset(2) = vim(19)
    if (ndim .eq. 3) offset(3) = vim(20)
!
! LA LDC EST DEFINIE PAR RAPPORT A NOUVEAU POINT D'EQUILIBRE
    do i = 1, ndim
        delta(i) = delta(i)-offset(i)
    end do
!
!
! CALCUL DES PENTES
!------------------
!
! PENTE ACTUEL EN DECHARGE DANS LA ZONE DE TRACTION
    if (lc*ka .ne. 0.d0) then
        rk = max(0.d0, sc*(1.d0-ka/lc)/ka)
    else
        rk = 0.d0
    end if
!
! DANS LE DOMAINE DE COMPRESSION
    rc = beta*r0
!
! PENTE TANGENTIELLE ACTUELLE
! SI ALPHA=2 RT CSTE, SI ALPHA=0 ALORS LCT=0 ET DONC RT=0)
    rt = rt0
    if (delta(1) .gt. 0.d0) then
        if (lct .ne. 0.d0) rt = max(0.d0, rt0*(1.d0-delta(1)/lct))
        if (lct .eq. 0.d0) rt = 0.d0
    end if
    if (alpha .eq. 2.d0) rt = rt0
!
! INITIALISATION COMPLEMENTAIRE POUR RIGI_MECA_TANG (SECANTE PENALISEE)
    if (.not. resi) then
        if (elas) then
            diss = 0
        else
            diss = nint(vim(2))
        end if
        cass = nint(vim(3))
        goto 5000
    end if
!
!     INITIALISATION DE LA CONTRAINTE
    call r8inir(6, 0.d0, sigma, 1)
!
! CALCUL DE LA CONTRAINTE HYDRO : DEBIT (LOI CUBIQUE)
! RECUP DE LA PRESSION AU PG
!----------------------------------------------------
!
    if (ifhyme) then
        do n = 1, ndim-1
            sigma(ndim+n) = -rhof*gp(n)*(max(amin, delta(1)+amin))**3/(12*visf)
        end do
    end if
!
! CALCUL DE LA CONTRAINTE MECANIQUE
!----------------------------------
!    CONTRAINTE DE CONTACT PENALISE
!    ET PRISE EN COMPTE DE LA PRESSION DE FLUIDE EVENTUELLE
!    PRESFL : IMPOSEE, PRESG : CALCULEE (MODELISATION HYME)
!
    if (ifhyme) then
        sigma(1) = rc*min(0.d0, delta(1))-presg
    else
        sigma(1) = rc*min(0.d0, delta(1))-presfl
    end if
!
!    PARTIE TANGENTIELLE
    do i = 2, ndim
        if (rt .gt. 0.d0) then
!           sigma(i) = sigmo(i) + rt*ddelta(i) ! version incrementale
            sigma(i) = rt*delta(i)
        else
            sigma(i) = 0.d0
!          glissement de joint de la valeur de delta
            offset(i) = delta(i)+offset(i)
        end if
    end do
!
!    CONTRAINTE DE FISSURATION NORMALE
    if ((delta(1) .ge. lc) .or. (ka .ge. lc)) then
        diss = 0
        cass = 2
    else
        if (delta(1) .le. ka) then
!
            diss = 0
            if (ka .gt. k0) then
                cass = 1
            else
                cass = 0
            end if
            sigma(1) = sigma(1)+rk*max(0.d0, delta(1))
!
        else
!
            diss = 1
            cass = 1
            if (lc .ne. 0.d0) then
                ra = max(0.d0, sc*(1.d0-delta(1)/lc)/delta(1))
            else
                ra = 0.d0
            end if
            sigma(1) = sigma(1)+ra*max(0.d0, delta(1))
!
        end if
    end if
!
! ACTUALISATION DES VARIABLES INTERNES
!-------------------------------------
! V1 : SEUIL, PLUS GRANDE NORME DU SAUT
! V2 : INDICATEUR DE DISSIPATION (0 : NON, 1 : OUI)
! V3 : INDICATEUR D'ENDOMMAGEMENT NORMAL (0 : SAIN, 1: ENDOM, 2: CASSE)
! V4 : POURCENTAGE D'ENDOMMAGEMENT NORMAL (DANS LA ZONE ADOUCISSANTE)
! V5 : INDICATEUR D'ENDOMMAGEMENT TANGENTIEL (0:SAIN, 1:ENDOM, 2:CASSE)
! V6 : POURCENTAGE D'ENDOMMAGEMENT TANGENTIEL
! V7 A V9 : VALEURS DU SAUT DANS LE REPERE LOCAL
! V10: EPAISSEUR DU JOINT
! V11 : CONTRAINTE MECANIQUE NORMALE (SANS INFLUENCE PRESSION DE FLUIDE)
! V12 A V14 : COMPOSANTES DU GRADIENT DE PRESSION DANS LE REPERE GLOBAL
! V15 A V17 : COMPOSANTES DU FLUX HYDRO DANS LE REPERE GLOBAL
! V18 : PRESSION DE FLUIDE IMPOSEE OU CALCULEE ET INTERPOLEE (EN HYME)
! V19-V20 : GLISSEMENT TANGENTIELS
    vip(1) = max(ka, delta(1))
    vip(2) = diss
    vip(3) = cass
    if (lc .ne. 0.d0) then
        tmp = max(0.d0, (vip(1)-val(2)/val(1))/(lc-val(2)/val(1)))
        vip(4) = min(1.d0, tmp)
    else
        vip(4) = 1.d0
    end if
    vip(5) = 0.d0
    if (rt .lt. rt0) vip(5) = 1.d0
    if (rt .eq. 0.d0) vip(5) = 2.d0
    vip(6) = 1.d0-rt/rt0
    vip(7) = delta(1)+offset(1)
    vip(8) = delta(2)+offset(2)
    if (ndim .eq. 3) then
        vip(9) = delta(3)+offset(3)
    else
        vip(9) = 0.d0
    end if
!
!     CALCUL DU NOUVEAU POINT D'EQUILIBRE V10 EN CAS DE CLAVAGE/SCIAGE
!     LE CLAVAGE FAIT AUGMENTER L'EPAISSEUR DU JOINT
!     => OFFSET(1) EST CROISSANT (CLAVAGE)
!     LE SCIAGE FAIT DIMINUER L'EPAISSEUR DU JOINT
!     => OFFSET(1) EST DECROISSANT (SCIAGE)
    vip(10) = offset(1)
!
!     FLUX, GRAD DE PRESSION ET PRESSION DANS LE REPERE GLOBAL
    if (ifhyme) then
        gploc(1) = 0.d0
        gploc(2) = gp(1)
        if (ndim .eq. 3) then
            gploc(3) = gp(2)
        end if
!
        fhloc(1) = 0.d0
        fhloc(2) = sigma(ndim+1)
        if (ndim .eq. 3) then
            fhloc(3) = sigma(2*ndim-1)
        end if
!
        call matinv('S', ndim, coorot(ndim+1), invrot, r8bid)
        call pmavec('ZERO', ndim, invrot, gploc, gpglo)
        call pmavec('ZERO', ndim, invrot, fhloc, fhglo)
!
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE CALCULEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presg
        vip(12) = gpglo(1)
        vip(13) = gpglo(2)
        vip(15) = fhglo(1)
        vip(16) = fhglo(2)
        if (ndim .eq. 3) then
            vip(14) = gpglo(3)
            vip(17) = fhglo(3)
        else
            vip(14) = 0.d0
            vip(17) = 0.d0
        end if
!
!       PRESSION DE FLUIDE CALCULEE AUX NOEUDS (DDL) ET INTERPOL AU PG
        vip(18) = presg
    else
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE IMPOSEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presfl
!       VI PAS UTILISEES EN MODELISATION NON HYME
        vip(12) = 0.d0
        vip(13) = 0.d0
        vip(14) = 0.d0
        vip(15) = 0.d0
        vip(16) = 0.d0
        vip(17) = 0.d0
!       PRESSION DE FLUIDE IMPOSEE AU PG :
        vip(18) = presfl
    end if
    vip(19) = offset(2)
    if (ndim .eq. 3) then
        vip(20) = offset(3)
    else
        vip(20) = 0.d0
    end if
!
5000 continue
!
! INITIALISATION DE LA MATRICE TANGENTE
    call r8inir(6*6, 0.d0, dsidep, 1)
!
    if (.not. rigi) goto 999
!
! CALCUL DE LA MATRICE TANGENTE HYDRO
!------------------------------------
    if (ifhyme) then
!
!       TERME : DW/DGP  (POUR KTAN P P)
        do n = 1, ndim-1
            dsidep(ndim+n, ndim+n) = -rhof*(max(amin, delta(1)+amin))**3/(12*visf)
        end do
!
!       TERME : DW/DDELTA_N  (POUR KTAN P U)
        do n = 1, ndim-1
            if (delta(1) .lt. 0.d0) then
                dsidep(ndim+n, 1) = 0.d0
            else
                dsidep(ndim+n, 1) = -3*rhof*gp(n)*(delta(1)+amin)**2/(12*visf)
            end if
        end do
!
    end if
!
! CALCUL DE LA MATRICE TANGENTE MECA (POUR KTAN U U)
!-----------------------------------
!
    rigart = 1.d-8
!   MATRICE TANGENTE DE CONTACT FERME
    if (delta(1) .le. 0.d0) then
        dsidep(1, 1) = rc
! POUR LE JOINT CLAVE LA MATRICE DE RIGIDITE NORMALE EST ZERO
        if ((prescl .ge. 0.d0) .and. (doffset(1) .gt. 0.d0)) dsidep(1, 1) = rigart*r0
        do i = 2, ndim
            dsidep(i, i) = rt0
        end do
    else
!   MATRICE TANGENTE DE CONTACT OUVERT
!   (NOTER QUE DANS LA SUITE DELTA(1)>0)
!
!       MATRICE TANGENTE DE FISSURATION
        if ((diss .eq. 0) .or. elas) then
            dsidep(1, 1) = rk
        else
            if (lc .ne. 0.d0) dsidep(1, 1) = -sc/lc
        end if
!
        do i = 2, ndim
            dsidep(i, i) = rt
            if ((lct .ne. 0.d0) .and. (delta(1) .lt. lct)) then
                dsidep(i, 1) = -delta(i)*rt0/lct
            end if
        end do
!
!       DANS LE CAS OU L'ELEMENT EST TOTALEMENT CASSE ON INTRODUIT UNE
!       RIGIDITE ARTIFICIELLE DANS LA MATRICE TANGENTE POUR ASSURER
!       LA CONVERGENCE
        if (cass .eq. 2) dsidep(1, 1) = rigart*r0
!
        if (abs(rt) .lt. rigart*rt0) then
            do i = 2, ndim
                dsidep(i, i) = rigart*rt0
            end do
        end if
    end if
!
999 continue
end subroutine
