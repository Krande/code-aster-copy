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
!
subroutine lcumfp(fami, kpg, ksp, ndim, typmod, &
                  imate, compor, tinstm, tinstp, epsm, &
                  deps, sigm, vim, option, rela_plas, &
                  sigp, vip, dsidep)
!
    implicit none
!
#include "asterfort/lceibt.h"
#include "asterfort/lcldsb.h"
#include "asterfort/lcmaza.h"
#include "asterfort/lcumef.h"
#include "asterfort/lcummd.h"
#include "asterfort/lcumme.h"
#include "asterfort/lcumsf.h"
#include "asterfort/lcumvi.h"
#include "asterfort/mgauss.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/sigela.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/get_varc.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
!
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imate
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    character(len=8), intent(in) :: typmod(*)
    character(len=16), intent(in) :: compor(*)
    character(len=16), intent(in) :: rela_plas
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: fami
    real(kind=8) :: tinstm, tinstp
    real(kind=8) :: epsm(*), deps(*), sigm(*), sigp(*), vim(*), vip(*)
    real(kind=8) :: dsidep(6, 6), tbid(36)
!
!---&s---1---------2---------3---------4---------5---------6---------7--
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION MECANIQUE   (--> IFOU SOUS CASTEM)
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  TINSTM  : INSTANT AU CALCUL PRECEDENT
! IN  TINSTP  : INSTANT DU CALCUL
! IN  DEPS    : INCREMENT DE DEFORMATION TOTALE
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : (1) OPTION DEMANDEE:RIGI_MECA_TANG, FULL_MECA ,RAPH_MECA
!               (2) MODELE MECA DE COUPLAGE EVENTUEL (MAZARS OU EIB)
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE
!
! ANCIENNE ROUTINE CASTEM FLU     SOURCE    BENBOU   02/02/7
! UMLVFP SOURCE YLP septembre 2002
!_______________________________________________________________________
!
! ROUTINE CALCULANT :
!
!    - CONTRAINTES FINALES          : SIGP (NSTRS)
!    - VARIABLES INTERNES FINALES   : VIP (NVARI)
!_______________________________________________________________________
!
! STRUCTURE DES PARAMETRES MATERIAU ET AUTRE
! PARAMETRES ELASTIQUES
!     CMAT(1)     = YOUN : MODULE D YOUNG
!     CMAT(2)     = XNU  : COEFFICIENT DE POISSON ELASTIQUE
! PARAMETRES DU FLUAGE PROPRE
!     CMAT(3)     = KRS   : RIGIDITE SPHERIQUE APPARENTE SQUELETTE
!     CMAT(4)     = ETARS : VISCOSITE SPHERIQUE APPARENTE EAU LIBRE
!     CMAT(5)     = KIS   : RIGIDITE SPHERIQUE APPARENTE HYDRATES
!     CMAT(6)     = ETAIS : VISCOSITE SPHERIQUE APPARENTE EAU LIEE
!     CMAT(7)     = KRD   : RIGIDITE DEVIATORIQUE APPARENTE
!     CMAT(8)     = ETARD : VISCOSITE DEVIATORIQUE APPARENTE EAU LIBRE
!     CMAT(9)     = ETAID : VISCOSITE DEVIATORIQUE APPARENTE EAU LIEE
! LES DEUX PARAMETRES SUIVANTS CONCERNENT UNIQUEMENT
! LE FLUAGE DE DESSICCATION
!     CMAT(10)    = RDES : COEFFICIENT UMLV FLUAGE DESSICCATION
!     CMAT(11)    = VDES : COEFFICIENT UMLV/BAZANT FLUAGE DESSICCATION
! OBSOLETE (CONCERNE UNIQUEMENT CASTEM) --> TYPMOD (CODE_ASTER)
!     CMAT(12)    = IFOU : HYPOTHESE DE CALCUL AUX ELEMENTS FINIS
!                     -2 : CONTRAINTES PLANES
!                     -1 : DEFORMATION PLANE
!                      0 : AXISYMETRIQUE
!                      2 : TRIDIMENSIONEL
! PAR PRINCIPE CMAT(13) = 2
!     CMAT(13)    = IFPO : OPTION SUR LE MODELE FLUAGE PROPRE
!                      0 : PAS DE FLUAGE PROPRE
!                      1 : PAS D INFUENCE DE L HUMIDITE REALTIVE
!                      2 : INFLUENCE DE L HUMIDITE RELATIVE
! OBSOLETE --> DEBRANCHE
!     CMAT(14)    = IDES : OPTION SUR LE MODELE FLUAGE DESSICCATION
!                      0 : PAS PRIS EN COMPTE
!                      1 : MODELE BAZANT
!                      2 : MODELE UMLV
!     CMAT(15)    = ICOU : OPTION COUPLAGE MECANIQUE/FLUAGE
!                      0 : PAS DE COUPLAGE (CALCUL CHAINE)
!                      1 : COUPLAGE FORT
!_______________________________________________________________________
!
!  STRUCTURE DES CONTRAINTES : SIGM,SIGP ( X = M ou P )
!    IFOU = -2 : CONTRAINTES PLANES
!      - SIGX(1) = SIGMA_XX
!      - SIGX(2) = SIGMA_YY
!      - SIGX(3) = SIGMA_XY
!      - (SIGX(4) = SIGMA_ZZ = 0)
!
!    IFOU = -1 : DEFORMATION PLANE
!      - SIGX(1) = SIGMA_XX
!      - SIGX(2) = SIGMA_YY
!      - SIGX(3) = SIGMA_ZZ
!      - SIGX(4) = SIGMA_XY
!
!    IFOU = 0  : AXISYMETRIQUE
!      - SIGX(1) = SIGMA_RR
!      - SIGX(2) = SIGMA_ZZ
!      - SIGX(3) = SIGMA_TT
!      - SIGX(4) = SIGMA_RZ
!
!    IFOU = 2  : TRIDIMENSIONEL
!      - SIGX(1) = SIGMA_XX
!      - SIGX(2) = SIGMA_YY
!      - SIGX(3) = SIGMA_ZZ
!      - SIGX(4) = SIGMA_XY
!      - SIGX(5) = SIGMA_ZX
!      - SIGX(6) = SIGMA_YZ
!_______________________________________________________________________
!
!  STRUCTURE DE L'INCREMENT DEFORMATION TOTALE : DEPS (NSTRS)
!
!    IFOU = -2 : CONTRAINTES PLANES
!      - DEPS(1) = DEPSILON_XX
!      - DEPS(2) = DEPSILON_YY
!      - DEPS(3) = DEPSILON_XY
!      - DEPS(4) = DEPSILON_ZZ
!
!    IFOU = -1 : DEFORMATION PLANE
!      - DEPS(1) = DEPSILON_XX
!      - DEPS(2) = DEPSILON_YY
!      - DEPS(3) = DEPSILON_XY
!      - (DEPS(4) = DEPSILON_ZZ = 0)
!
!    IFOU = 0  : AXISYMETRIQUE
!      - DEPS(1) = DEPSILON_RR
!      - DEPS(2) = DEPSILON_ZZ
!      - DEPS(3) = DEPSILON_TT
!      - DEPS(4) = DEPSILON_RZ
!
!    IFOU = 2  : TRIDIMENSIONEL
!      - DEPS(1) = SIGMA_XX
!      - DEPS(2) = SIGMA_YY
!      - DEPS(3) = SIGMA_ZZ
!      - DEPS(4) = SIGMA_XY
!      - DEPS(5) = SIGMA_ZX
!      - DEPS(6) = SIGMA_YZ
!_______________________________________________________________________
!
!  STRUCTURE DES VARIABLES INTERNES : VIM,VIP ( X = I ou F )
!
!     VIX(1)     = ERSP  : DEFORMATION DE FLUAGE REV SPHERIQUE
!     VIX(2)     = EISP  : DEFORMATION DE FLUAGE IRR SPHERIQUE
!     VIX(3)     = ERD11 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 11
!     VIX(4)     = EID11 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 11
!     VIX(5)     = ERD22 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 22
!     VIX(6)     = EID22 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 22
!     VIX(7)     = ERD33 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 33
!     VIX(8)     = EID33 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 33
!     VIX(9)     = EFD11 : DEFORMATION DE FLUAGE DE DESSICCATION  11
!     VIX(10)    = EFD22 : DEFORMATION DE FLUAGE DE DESSICCATION  22
!     VIX(11)    = EFD33 : DEFORMATION DE FLUAGE DE DESSICCATION  33
!     VIX(12)    = ERD12 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 12
!     VIX(13)    = EID12 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 12
!     VIX(14)    = ERD23 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 23
!     VIX(15)    = EID23 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 23
!     VIX(16)    = ERD31 : DEFORMATION DE FLUAGE REV DEVIATORIQUE 31
!     VIX(17)    = EID31 : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 31
!     VIX(18)    = EFD12 : DEFORMATION DE FLUAGE DE DESSICCATION  12
!     VIX(19)    = EFD23 : DEFORMATION DE FLUAGE DE DESSICCATION  23
!     VIX(20)    = EFD31 : DEFORMATION DE FLUAGE DE DESSICCATION  31
!     VIX(21)    = INDICATEUR DU FLUAGE SPHERIQUE (0 ou 1)
!     VARIABLES INTERNES UTILISEES SI COUPLAGE AVEC ENDO_ISOT_BETON
!     VIX(22)    = ENDOMMAGEMENT D DONNE PAR EIB
!     VIX(23)    = INDICATEUR D'ENDOMMAGEMENT DE EIB
!     VARIABLES INTERNES UTILISEES SI COUPLAGE AVEC MAZARS
!     VIX(22)    = ENDOMMAGEMENT D DONNE PAR MAZARS
!     VIX(23)    = INDICATEUR D'ENDOMMAGEMENT DE EIB
!     VIX(24)    = TEMPERATURE MAXIMALE ATTEINTE PAR LE MATERIAU
!     VIX(25)    = VALEUR DE EPSEQ (UTILE POUR POSTTRAITER)
!_______________________________________________________________________
!
    character(len=16) :: option2
    real(kind=8) :: det
    integer(kind=8) :: iret
    character(len=16) :: nomres(16), phenbid
    integer(kind=8) :: icodre(16)
    real(kind=8) :: cfps, cfpd
    integer(kind=8) :: i, j, k, l, nstrs, ifou, isph
    real(kind=8) :: tdt
    real(kind=8) :: youn, xnu
    real(kind=8) :: bendo, kdess
    real(kind=8) :: krs, etars, kis, etais, krd, etard, etaid
    real(kind=8) :: etafd
    real(kind=8) :: cmat(15), dep(6, 12), depm(6, 6)
    real(kind=8) :: an(6), bn(6, 6), cn(6, 6), valres(16)
    real(kind=8) :: hygrm, hygrp, rbid
    real(kind=8) :: matn(6, 6), invn(6, 6), eps(6), epsf(6)
    real(kind=8) :: epsrm, epsrp, epsfm(6)
    real(kind=8) :: epsme(6), epse(6)
    real(kind=8) :: hydrm, hydrp, sechm, sechp, sref, tm, tp, tref
    real(kind=8) :: epsthp, epsthm
!
    real(kind=8) :: tmaxp, tmaxm, younm, xnum
    real(kind=8) :: sigelm(6), sigelp(6), epsel(6)
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    blas_int :: b_incx, b_incy, b_n
!
!
!   CALCUL DE L'INTERVALLE DE TEMPS
!
    tdt = tinstp-tinstm
!
!   DIMENSION
!
    nstrs = 2*ndim
!
!   TYPE DE CALCUL
!
    if (typmod(1) .eq. 'C_PLAN') then
        ifou = -2
    else if (typmod(1) .eq. 'D_PLAN') then
        ifou = -1
    else if (typmod(1) .eq. 'AXIS') then
        ifou = 0
    else
        ifou = 2
    end if
!
!   INITIALISATION DU FLUAGE SPHERIQUE PROPRE
!
    isph = 1
!
! - Get temperatures
!
    call get_varc(fami, kpg, ksp, 'T', tm, &
                  tp, tref)
!
!
!  ------- LECTURE DES CARACTERISTIQUES ELASTIQUES
!  MB: LA DEPENDENCE DES PARAMETRES PAR RAPPORT A LA TEMPERATURE
!  CHANGE PAR RAPPORT A LA LOI D ENDOMMAGEMENT (COUPLAGE)
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    nomres(4) = 'ALPHA'
!
    if (rela_plas .eq. 'ENDO_ISOT_BETON') then
!
        call rcvalb(fami, 1, 1, '+', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [0.d0], &
                    2, nomres, valres, icodre, 1)
!
        call rcvalb(fami, 1, 1, '+', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [0.d0], &
                    1, nomres(3), valres(3), icodre(3), 0)
        valres(4) = valres(3)
        icodre(4) = icodre(3)
!
    else if (rela_plas .eq. 'MAZARS') then
        tmaxm = vim(24)
        tmaxp = max(tmaxm, tp)
!
        call rcvalb(fami, 1, 1, '-', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [tmaxm], &
                    2, nomres, valres, icodre, 1)
        younm = valres(1)
        xnum = valres(2)
!
        call rcvalb(fami, 1, 1, '+', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [tmaxp], &
                    2, nomres, valres, icodre, 1)
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [tmaxm], &
                    1, nomres(3), valres(3), icodre(3), 0)
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ELAS', 1, 'TEMP', [tmaxp], &
                    1, nomres(4), valres(4), icodre(4), 0)
!
    else
        call rcvalb(fami, kpg, ksp, '-', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
        younm = valres(1)
        xnum = valres(2)
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
    end if
!
    youn = valres(1)
    xnu = valres(2)
!
!
!  -------CALCUL DES DEFORMATIONS THERMIQUES
!
    if ((rela_plas .eq. 'MAZARS') .or. (rela_plas .eq. 'ENDO_ISOT_BETON')) then
        if ((isnan(tref)) .or. (icodre(3) .ne. 0) .or. (icodre(4) .ne. 0)) then
            call utmess('F', 'CALCULEL_15')
        else
            if (.not. isnan(tm)) then
                epsthm = valres(3)*(tm-tref)
            else
                epsthm = 0.d0
            end if
            if (.not. isnan(tp)) then
                epsthp = valres(4)*(tp-tref)
            else
                epsthp = 0.d0
            end if
        end if
    else
!
        call verift(fami, kpg, ksp, '+', imate, &
                    epsth_=epsthp)
        call verift(fami, kpg, ksp, '-', imate, &
                    epsth_=epsthm)
    end if
!
!  ------- CARACTERISTIQUES DE RETRAIT ENDOGENE ET DE DESSICCATION
!
    nomres(1) = 'B_ENDOGE'
    nomres(2) = 'K_DESSIC'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres, valres, icodre, 1)
    bendo = valres(1)
    kdess = valres(2)
!
!  ------- CARACTERISTIQUES FLUAGE PROPRE UMLV
!
    nomres(1) = 'K_RS'
    nomres(2) = 'ETA_RS'
    nomres(3) = 'K_IS'
    nomres(4) = 'ETA_IS'
    nomres(5) = 'K_RD'
    nomres(6) = 'ETA_RD'
    nomres(7) = 'ETA_ID'
!
    rbid = 0.0d0
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'BETON_UMLV', 0, ' ', [rbid], &
                7, nomres, valres, icodre, 2)
    krs = valres(1)
    etars = valres(2)
    kis = valres(3)
    etais = valres(4)
    krd = valres(5)
    etard = valres(6)
    etaid = valres(7)
!
! ------- CARACTERISTIQUE FLUAGE DE DESSICATION DE BAZANT
!
    nomres(8) = 'ETA_FD'
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'BETON_UMLV', 0, ' ', [rbid], &
                8, nomres, valres, icodre, 0)
!     FLUAGE DE DESSICCATION NON ACTIVE
    if (icodre(8) .ne. 0) then
        cmat(14) = 0
        etafd = -1.0d0
!     FLUAGE DE DESSICCATION ACTIVE
    else
        cmat(14) = 1
        etafd = valres(8)
    end if
!
!  ------- CARACTERISTIQUES HYGROMETRIE H
!
    call rccoma(imate, 'BETON_DESORP', 0, phenbid, icodre(1))
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH4_93')
    end if

    nomres(1) = 'FONC_DESORP'
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'BETON_DESORP', 0, ' ', [rbid], &
                1, nomres(1), valres(1), icodre(1), 0)
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH4_94')
    end if
    hygrm = valres(1)
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'BETON_DESORP', 0, ' ', [rbid], &
                1, nomres(1), valres(1), icodre(1), 0)
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ALGORITH4_94')
    end if
    hygrp = valres(1)
!
! CONSTRUCTION DU VECTEUR CMAT CONTENANT LES CARACTERISTIQUES MECANIQUES
!
!     CMAT(1)     = YOUN   : MODULE D YOUNG
!     CMAT(2)     = XNU    : COEFFICIENT DE POISSON ELASTIQUE
!     CMAT(3)     = KRS   : RIGIDITE SPHERIQUE APPARENTE SQUELETTE
!     CMAT(4)     = ETARS : VISCOSITE SPHERIQUE APPARENTE EAU LIBRE
!     CMAT(5)     = KIS   : RIGIDITE SPHERIQUE APPARENTE HYDRATES
!     CMAT(6)     = ETAIS : VISCOSITE SPHERIQUE APPARENTE EAU LIEE
!     CMAT(7)     = KRD   : RIGIDITE DEVIATORIQUE APPARENTE
!     CMAT(8)     = ETARD : VISCOSITE DEVIATORIQUE APPARENTE EAU LIBRE
!     CMAT(9)     = ETAID : VISCOSITE DEVIATORIQUE APPARENTE EAU LIEE
!
    cmat(1) = youn
    cmat(2) = xnu
    cmat(3) = krs
    cmat(4) = etars
    cmat(5) = kis
    cmat(6) = etais
    cmat(7) = krd
    cmat(8) = etard
    cmat(9) = etaid
! MODIFI FD BAZANT
    cmat(11) = etafd
    cmat(12) = ifou
    cmat(13) = 2
! MODIFI 25/08/04 YLP - ACTIVATION DU FLUAGE DE DESSICATION DE BAZANT
! CMAT(14)=IDES 0 --> 1
!      CMAT(14)    = 1
    cmat(15) = 1
!
!   DANS LE CAS OU LE TEST DE DEFORMATION DE FLUAGE PROPRE
!        IRREVE A ECHOUE : ISPH = 0
!
10  continue
!
! INITIALISATION DES VARIABLES
!
    cfps = 0.d0
    cfpd = 0.d0
!
    do i = 1, 6
        an(i) = 0.d0
        do j = 1, 6
            dep(i, j) = 0.d0
            bn(i, j) = 0.d0
            cn(i, j) = 0.d0
        end do
    end do
!_______________________________________________________________________
!
! CALCUL DES MATRICES DES DEFORMATIONS DE FLUAGE TOTAL
!   DFLUT(N+1) = AN + BN * SIGMA(N) + CN * SIGMA(N+1)
!_______________________________________________________________________
!
    if (tdt .ne. 0.d0) then
        if (option(1:9) .eq. 'RIGI_MECA') then
            isph = nint(vim(21))
        end if
        call lcummd(vim, 20, cmat, 15, sigm, &
                    nstrs, isph, tdt, hygrm, hygrp, &
                    an, bn, cn, cfps, cfpd)
    end if
!
!_______________________________________________________________________
!
! RECUPERATION DE L HYDRATATION E DU SECHAGE
! CALCUL DE LA SIGMA ELASTIQUE AU TEMP M POUR COUPLAGE AVEC MAZARS
!_______________________________________________________________________
!
!
    call lcumvi('FT', vim, epsfm)
!
!
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')) then
        call rcvarc(' ', 'HYDR', '+', fami, kpg, &
                    ksp, hydrp, iret)
        if (iret .ne. 0) hydrp = 0.d0
        call rcvarc(' ', 'HYDR', '-', fami, kpg, &
                    ksp, hydrm, iret)
        if (iret .ne. 0) hydrm = 0.d0
        call rcvarc(' ', 'SECH', '+', fami, kpg, &
                    ksp, sechp, iret)
        if (iret .ne. 0) sechp = 0.d0
        call rcvarc(' ', 'SECH', '-', fami, kpg, &
                    ksp, sechm, iret)
        if (iret .ne. 0) sechm = 0.d0
        call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                    ksp, sref, iret)
        if (iret .ne. 0) sref = 0.d0
!
        epsrm = kdess*(sechm-sref)-bendo*hydrm+epsthm
        epsrp = kdess*(sechp-sref)-bendo*hydrp+epsthp
!
!
! - MB: CALCUL DE LA DEFORMATION ELASTIQUE AU TEMP M
!    (LA SEULE QUI CONTRIBUE A FAIRE EVOLUER L'ENDOMMAGEMENT)
!    POUR LE COUPLAGE AVEC MAZARS
!
        if (rela_plas .eq. 'MAZARS') then
            call r8inir(6, 0.d0, epsel, 1)
            do k = 1, nstrs
                epsel(k) = epsm(k)-epsrm*kron(k)-epsfm(k)
            end do
!
!  -  ON CALCUL LES CONTRAINTES ELASTIQUES AU TEMP M
!          CALL SIGELA (NDIM,'LAMBD',LAMBDA,DEUXMU,EPSEL,SIGELM)
!
            call sigela(typmod, ndim, younm, xnum, epsel, &
                        sigelm)
        end if
!
! ________________________________________________________________
!
!  1. CONSTRUCTION DE LA MATRICE D ELASTICITE DE HOOKE POUR MAZARS
!     OU UMLV SANS COUPLAGE, OU DE LA MATRICE ELASTO-ENDOMMAGEE POUR EIB
!  2. MISE A JOUR DE L ENDOMMAGEMENT ET DES SIGMA POUR EIB
! ________________________________________________________________
!
        if (rela_plas .eq. 'ENDO_ISOT_BETON') then
!    MATRICE ELASTO-ENDOMMAGEE ET MISE A JOUR DE L ENDOMMAGEMENT
            call lcldsb(fami, kpg, ksp, ndim, imate, &
                        epsm, deps, vim(22), 'RAPH_COUP       ', tbid, &
                        vip(22), dep)
        else
!    MATRICE D ELASTICITE DE HOOKE POUR MAZARS ET UMLV SANS COUPLAGE
            if (rela_plas .eq. 'MAZARS') then
                call lcumme(youn, xnu, ifou, dep)
            else
                call lcumme(youn, xnu, ifou, dep)
                call lcumme(younm, xnum, ifou, depm)
            end if
        end if
! ________________________________________________________________
!   MODIFIE MB 20 OCT 2008
!
!  1. MISE A JOUR DES SIGMA POUR EIB ET UMLV SANS COUPLAGE
!     CALCUL DES SIGMA ELASTIQUES POUR MAZARS
!      (LCUMEF)
!  2. MISE A JOUR DES VARIABLES INTERNES FINALES DE FLUAGE
!      (LCUMSF)
! ________________________________________________________________
!
!  PRISE EN COMPTE DU FLUAGE PROPRE ET DE DESSICCATION
!   MODIFI DU 18 AOUT 2004 YLP - CORRECTION DE LA DEFORMATION DE FLUAGE
!   PAR LES DEFORMATIONS DE RETRAIT
!
        if (rela_plas .eq. 'MAZARS') then
            call lcumef(rela_plas, dep, dep, an, bn, &
                        cn, epsm, epsrm, epsrp, deps, &
                        epsfm, sigelm, nstrs, sigelp)
            call lcumsf(sigelm, sigelp, nstrs, vim, 20, &
                        cmat, 15, isph, tdt, hygrm, &
                        hygrp, vip)
        else if (rela_plas .eq. 'ENDO_ISOT_BETON') then
            call lcumef(rela_plas, dep, dep, an, bn, &
                        cn, epsm, epsrm, epsrp, deps, &
                        epsfm, sigm, nstrs, sigp)
            call lcumsf(sigm, sigp, nstrs, vim, 20, &
                        cmat, 15, isph, tdt, hygrm, &
                        hygrp, vip)
        else
            call lcumef(rela_plas, dep, depm, an, bn, &
                        cn, epsm, epsrm, epsrp, deps, &
                        epsfm, sigm, nstrs, sigp)
            call lcumsf(sigm, sigp, nstrs, vim, 20, &
                        cmat, 15, isph, tdt, hygrm, &
                        hygrp, vip)
        end if
        vip(21) = 1
!
!  TEST DE LA CROISSANCE SUR LA DEFORMATION DE FLUAGE PROPRE SPHERIQUE
!
        if (isph .eq. 2) then
            isph = 0
            goto 10
        end if
!
!___________________________________________________________
!
!  MB: MISE A JOUR DE L ENDOMMAGEMENT ET DES SIGMA POUR MAZARS
!_________________________________________________________
!
!
        if (rela_plas .eq. 'MAZARS') then
            option2 = 'RAPH_COUP'
            call lcmaza(fami, kpg, ksp, ndim, typmod, &
                        imate, compor, epsm, deps, vim(22), &
                        option2, sigp, vip, tbid)
        end if
    end if
!
!_______________________________________________________________________
!
! CONSTRUCTION DE LA MATRICE TANGENTE
!_______________________________________________________________________
!
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RIGI_MECA')) then
!
! - MB: SI COUPLAGE AVEC MAZARS, ON UTILISE POUR LE COUPLAGE
!       LA MATRICE TANGENTE DE CETTE LOI
        if (rela_plas .eq. 'MAZARS') then
            option2 = option
            if (option(1:9) .eq. 'FULL_MECA') then
                option2 = 'RIGI_COUP'
            end if
            call lcmaza(fami, kpg, ksp, ndim, typmod, &
                        imate, compor, epsm, deps, vim(22), &
                        option2, tbid, vip, dsidep)
        else
            option2 = 'RIGI_COUP'
            if (option(1:9) .eq. 'RIGI_MECA') then
                if (rela_plas .eq. 'ENDO_ISOT_BETON') then
                    call lcldsb(fami, kpg, ksp, ndim, imate, &
                                epsm, tbid, vim(22), option2, tbid, &
                                tbid, dep)
                else
                    call lcumme(youn, xnu, ifou, dep)
                end if
            end if
            call r8inir(36, 0.d0, matn, 1)
            do i = 1, nstrs
                matn(i, i) = 1.d0
            end do
!
            do i = 1, nstrs
                do j = 1, nstrs
                    do k = 1, nstrs
                        matn(i, j) = matn(i, j)+cn(i, k)*dep(k, j)
                    end do
                end do
            end do
!
            call r8inir(36, 0.d0, invn, 1)
!
            do i = 1, nstrs
                invn(i, i) = 1.d0
            end do
!
            call mgauss('NFVP', matn, invn, 6, nstrs, &
                        nstrs, det, iret)
!
            call r8inir(36, 0.d0, dsidep, 1)
!
            do i = 1, nstrs
                do j = 1, nstrs
                    do k = 1, nstrs
                        dsidep(i, j) = dsidep(i, j)+invn(k, j)*dep(i, k)
                    end do
                end do
            end do
!
            if (rela_plas .eq. 'ENDO_ISOT_BETON') then
                if (option .eq. 'RIGI_MECA_TANG') then
                    call rcvarc(' ', 'HYDR', '+', fami, kpg, &
                                ksp, hydrp, iret)
                    if (iret .ne. 0) hydrp = 0.d0
                    call rcvarc(' ', 'HYDR', '-', fami, kpg, &
                                ksp, hydrm, iret)
                    if (iret .ne. 0) hydrm = 0.d0
                    call rcvarc(' ', 'SECH', '+', fami, kpg, &
                                ksp, sechp, iret)
                    if (iret .ne. 0) sechp = 0.d0
                    call rcvarc(' ', 'SECH', '-', fami, kpg, &
                                ksp, sechm, iret)
                    if (iret .ne. 0) sechm = 0.d0
                    call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                                ksp, sref, iret)
                    if (iret .ne. 0) sref = 0.d0
                    if (nint(vim(23)) .eq. 1) then
                        epsrm = kdess*(sechm-sref)-bendo*hydrm+epsthm
                        do i = 1, nstrs
                            epsme(i) = epsm(i)-epsrm*kron(i)
                        end do
                        call lceibt(nstrs, epsme, epsfm, dep, invn, &
                                    cn, dsidep)
                    end if
                else if ((option .eq. 'RAPH_MECA') .or. (option .eq. 'FULL_MECA')) then
                    if (nint(vip(23)) .eq. 1) then
                        b_n = to_blas_int(nstrs)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        call dcopy(b_n, epsm, b_incx, eps, b_incy)
                        b_n = to_blas_int(nstrs)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        call daxpy(b_n, 1.d0, deps, b_incx, eps, &
                                   b_incy)
                        do i = 1, nstrs
                            epse(i) = epsm(i)+deps(i)-epsrp*kron(i)
                        end do
                        call lcumvi('FT', vip, epsf)
                        call lceibt(nstrs, epse, epsf, dep, invn, &
                                    cn, dsidep)
                    end if
                end if
            end if
!----------- CORRECTION POUR LES CONTRAINTES PLANES :
            if (ifou .eq. -2) then
                do k = 1, nstrs
                    if (k .ne. 3) then
                        do l = 1, nstrs
                            if (l .ne. 3) then
                                dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep&
                                               &(3, l)
                            end if
                        end do
                    end if
                end do
            end if
        end if
    end if
!
end subroutine
