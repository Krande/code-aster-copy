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

subroutine lkcomp(fami, kpg, ksp, typmod, imate, instam, instap, &
                  deps, sigm, vinm, &
                  option, sigp, vinp, dside, retcom, &
                  invi)
!
    implicit none
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lkcrip.h"
#include "asterfort/lkcriv.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lkdgde.h"
#include "asterfort/lkelas.h"
#include "asterfort/lkgamp.h"
#include "asterfort/lklmat.h"
#include "asterfort/lkoptg.h"
#include "asterfort/r8inir.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
#include "asterfort/get_varc.h"
    integer(kind=8) :: retcom, imate, invi
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    character(len=8) :: typmod(*)
    character(len=16) :: option
    real(kind=8) :: instam, instap
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6), vinm(invi)
    real(kind=8) :: sigp(6), vinp(invi)
    real(kind=8) :: dside(6, 6)
! --- MODELE LETK : LAIGLE ET KLEINE (CIH)  MODELE VISCOPLASTIQUE--
! =================================================================
! --- BUT : ROUTINE PRINCIPALE---------------------------
! =================================================================
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  MOD     : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VINM    : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VINP    : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDE   : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
! OUT RETCOM  : CODE RETOUR POUR LE REDECOUPAGE DU PAS DE TEMPS
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!=======================================================================
!=======================================================================
! --- ATTENTION : CHANGEMENT DE SIGNES DES CHAMPS DE CONTRAINTES ET DES
! ----DEFORMATIONS - DANS CE MODELE CONVENTION MECANIQUE DES SOLS A L
! ----OPPPOSE DE CELLES DE LA MECANIQUE DES MILIEUX CONTINUS - EN
! ----COMPRESSION LA CONTRAINTE EST POSITIVE ET EN CONTRACTANCE :
! ----DEFORMATION VOLUMIQUE POSITIVE
!=======================================================================
!=======================================================================
    integer(kind=8) :: nbmat, ndt, ndi, nvi, val, varv, i, k, matr
    integer(kind=8) :: iret
    integer(kind=8) :: indal
    aster_logical ::  l_temp, lVari
    real(kind=8) :: mun, un, zero, deux, trois
!      REAL*8        LGLEPS
    parameter(nbmat=90)
    real(kind=8) :: materd(nbmat, 2), materf(nbmat, 2)
    real(kind=8) :: dt, alpha, coef
    real(kind=8) :: sigml(6), sigpl(6), depml(6), depsth(6)
    real(kind=8) :: i1ml, sml(6), siim
    real(kind=8) :: iel, i1el, sel1(6)
    real(kind=8) :: dvml, devml(6), tm, tp, tref
    real(kind=8) :: dvml1, devml1(6)
    real(kind=8) :: sel(6), sigel(6)
    real(kind=8) :: seuilv, seuilp
    real(kind=8) :: ucrvm, seuvm, ucrpm, seupm
    real(kind=8) :: depsv(6), dgamv, dxivm, xipic
    real(kind=8) :: depsp(6), dgamp, xivm, dxip, dxiv
    real(kind=8) :: seuivm, ucrivm, ucrip, ucriv, irrev(6)
    real(kind=8) :: dsig(6), vecd(6)
    real(kind=8) :: de(6, 6), kk, mu
    real(kind=8) :: kron(6), vintr
    real(kind=8) :: somme, variTmp(9)
    character(len=3) :: matcst
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(mun=-1.d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
!      PARAMETER       (LGLEPS =  1.0D-8 )
! =================================================================
    common/tdim/ndt, ndi
! =================================================================
    data kron/un, un, un, zero, zero, zero/
!
    ASSERT(invi .eq. 9)
! - Flag to modify internal state variable
    lVari = L_VARI(option)

! =================================================================
! --- INITIALISATION LOCALES
! =================================================================

    dt = 0.d0
    retcom = 0
    variTmp = 0.d0

    materd = 0.d0
    materf = 0.d0
    alpha = 0.d0
    coef = 0.d0

    sigml = 0.d0
    sigpl = 0.d0
    depml = 0.d0
    depsth = 0.d0

    i1ml = 0.d0
    siim = 0.d0
    sml = 0.d0

    iel = 0.d0
    i1el = 0.d0
    sel1 = 0.d0

    dvml = 0.d0
    devml = 0.d0
    tm = 0.d0
    tp = 0.d0
    tref = 0.d0

    dvml1 = 0.d0
    devml1 = 0.d0
    sel = 0.d0
    sigel = 0.d0
    seuilv = 0.d0
    seuilp = 0.d0
    ucrvm = 0.d0
    seuvm = 0.d0
    ucrpm = 0.d0
    seupm = 0.d0

    depsv = 0.d0
    dgamv = 0.d0
    dxivm = 0.d0
    xipic = 0.d0
    depsp = 0.d0
    dgamp = 0.d0
    xivm = 0.d0
    dxip = 0.d0
    dxiv = 0.d0

    seuivm = 0.d0
    ucrivm = 0.d0
    ucrip = 0.d0
    ucriv = 0.d0
    irrev = 0.d0
    dsig = 0.d0
    vecd = 0.d0
    de = 0.d0
    kk = 0.d0
    mu = 0.d0
    vintr = 0.d0
    somme = 0.d0

!
!
    dt = instap-instam

!
! - Get temperatures
!
    call get_varc(fami, kpg, ksp, 'T', &
                  tm, tp, tref, l_temp)

! =================================================================
! --- RECUPERATION DES PARAMETRES DU MODELE -----------------------
! --- LES COEFFICIENTS MATERIAU N EVOLUENT PAS AVEC LE TEMPS-------
! =================================================================
!
    matcst = 'OUI'
    call lklmat(typmod(1), imate, nbmat, tm, materd, &
                materf, matcst, ndt, ndi, nvi, &
                indal)
    ASSERT(invi .eq. nvi)
!      SIGC   = MATERD(3,2)
    xipic = materd(18, 2)
    xivm = materd(20, 2)
! =================================================================
! --- CONVENTIONS DE SIGNE DU MODELE LAIGLE VISCOPLASTIQUE --------
! =================================================================
!
    do i = 1, ndt
        sigml(i) = mun*sigm(i)
        depml(i) = mun*deps(i)
    end do
! =================================================================
! --- DEFINITION DES INVARIANTS ET DU DEVIATEUR A L'INSTANT MOINS--
! =================================================================
!
    i1ml = trace(ndi, sigml)
!
    call lcdevi(sigml, sml)
!
    siim = norm2(sml(1:ndt))
!
! =================================================================
! ---PRISE EN COMPTE DE LA DILATATION THERMIQUE--------------------
! =================================================================
!
    alpha = materd(3, 1)
!
    if (l_temp) then
        coef = alpha*(tp-tref)-alpha*(tm-tref)
    else
        coef = zero
    end if
!
! =================================================================
! --- DEFINITION DES DEFORMATIONS VOLUMIQUES ET DEVIATORIQUES -----
! =================================================================
    dvml = 0.d0
!
    do k = 1, ndt
        depsth(k) = depml(k)
    end do
    do k = 1, 3
        depsth(k) = depsth(k)-coef
        dvml = dvml+depsth(k)
    end do
    do k = 1, ndt
        devml(k) = depsth(k)-dvml/3.d0*kron(k)
    end do
!
! =================================================================
! --- VERIFICATION D'UN ETAT INITIAL PLASTIQUEMENT ADMISSIBLE -----
! =================================================================
    somme = sum(vinm(1:nvi))
    if (abs(somme) .lt. r8prem()) then
        call lkcrip(i1ml, sml, vinm, nbmat, materd, &
                    ucrpm, seupm)
        if (seupm/materd(4, 1) .gt. 1.0d-6) then
            call utmess('F', 'ALGORITH2_2')
        end if
    end if
! =================================================================
! --- PREDICTION ELASTIQUE ----------------------------------------
! =================================================================
    call lkelas(ndi, ndt, nbmat, materd, depsth, &
                sigml, de, kk, mu)
!
    iel = i1ml+trois*kk*dvml
!
    do i = 1, ndt
        sel(i) = sml(i)+deux*mu*devml(i)
    end do
!
    do i = 1, ndt
        sigel(i) = sel(i)+iel/trois*kron(i)
    end do
!
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
! =================================================================
! --- CRITERE VISQUEUX --------------------------------------------
! =================================================================
! =================================================================
! --- CALCUL DE fv(SIGE, XIVM) ---CRITERE VISQUEUX MAX-------------
! =================================================================
        call lkcriv(xivm, iel, sel, vinm, nbmat, &
                    materd, ucrivm, seuivm)
!
!           IF (UCRIVM.LT.ZERO)  CALL UTMESS('F','COMPOR1_27')
!
!---- VARV : EN DESSOUS DU CRITERE VISQUEUX MAX : CONTRACTANCE: VARV=0
!---- VARV : AU DESSUS DU CRITERE VISQUEUX MAX  : DILATANCE:    VARV=1
!
!---- VAL  : INDICATEUR POUR LES LOIS DE DILALANCE
!----      : EN DESSOUS DU PIC ET POUR LA VISCOSITE : VAL = 0
!----      : AU DESSUS DU PIC  : VAL = 1
!
        if (seuivm .lt. zero) then
            varv = 0
        else
            varv = 1
        end if
!
        vintr = vinm(3)
!
! =================================================================
! --- CALCUL DE fv(SIGE, XIVM) ---CRITERE VISCOPLASTIQUE ---------
! =================================================================
        call lkcriv(vintr, iel, sel, vinm, nbmat, &
                    materd, ucriv, seuilv)
!
! --- VERIFICATION DU SIGNE DE U A L INSTANT MOINS AVANT ENTREE
! --- DANS LKDGDE
!
        call lkcriv(vintr, i1ml, sml, vinm, nbmat, &
                    materd, ucrvm, seuvm)
!
! =================================================================
! --- PAS DE VISCOSITE  -------------------------------------------
! =================================================================
        if (seuilv .lt. zero) then
            val = 0
            dgamv = zero
            dxiv = zero
            dvml1 = zero
!
            do i = 1, ndt
                depsv(i) = zero
                devml1(i) = zero
            end do
!
!---- XIV A T + DT ------------------------------------------------
!
            variTmp(3) = vinm(3)
!
!---- GAMMAV A T + DT ---------------------------------------------
!
            variTmp(4) = vinm(4)
!
! --  INDICATEUR DE VISCOSITE
            variTmp(6) = 0.d0
!
        else
! =================================================================
! --- VISCOSITE  --------------------------------------------------
! =================================================================
            val = 0
!
! -------------CALCUL DE DEPSV ET DE GAMMAV ----CRITERE VISQUEUX---
            call lkdgde(val, vintr, dt, seuilv, ucrvm, &
                        i1ml, sml, vinm, nbmat, materd, &
                        depsv, dgamv, iret)
            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if
!
            dvml1 = trace(ndi, depsv)
            call lcdevi(depsv, devml1)
!
! -------------DELTA XIV
!
            dxivm = xivm-vinm(3)
            dxiv = min(dgamv, dxivm)
!
!---- XIV A T + DT ------------------------------------------------
!
            variTmp(3) = vinm(3)+dxiv
!
!---- GAMMAV A T + DT ---------------------------------------------
!
            variTmp(4) = vinm(4)+dgamv
!
! --  INDICATEUR DE VISCOSITE
            variTmp(6) = 1.d0
!
        end if
!
! --- MISE A JOUR DE LA PREDICTION DE LA CONTRAINTE ---------------
!
        i1el = iel-trois*kk*dvml1
!
        do i = 1, ndt
            sel1(i) = sel(i)-deux*mu*devml1(i)
        end do
! =================================================================
! --- CRITERE ELASTOPLASTIQUE  ------------------------------------
! =================================================================
! --- VERIFICATION DU SIGNE DE U A L INSTANT MOINS AVANT ENTREE
! --- DANS LKGAMP et LKOPTG
!
        call lkcrip(i1ml, sml, vinm, nbmat, materd, &
                    ucrpm, seupm)
!
! =================================================================
! --- CALCUL DE fp(SIGE, XIPM) ---CRITERE ELASTOPLASTIQUE ---------
! =================================================================
        call lkcrip(i1el, sel1, vinm, nbmat, materd, &
                    ucrip, seuilp)
!
        if ((ucrip .lt. zero) .or. (ucrpm .lt. zero)) then
            retcom = 1
            goto 999
        end if
!
!==================================================================
!--------- ELASTICITE ---------------------------------------------
!==================================================================
        if (seuilp .lt. zero) then
            dgamp = zero
!
            do i = 1, ndt
                depsp(i) = zero
            end do
!
!---- REACTUALISATION DES CONTRAINTES -----------------------------
!
            do i = 1, ndt
                sigel(i) = sel1(i)+i1el/trois*kron(i)
                sigpl(i) = sigel(i)
            end do
!
! -------- DELTA XIP
!
            if (varv .eq. 0) then
!
!--------- CONTRACTANCE
!---------- ELASTICITE EN DESSOUS DU CRITERE VISQUEUX MAX
                dxip = zero
                variTmp(5) = 0.0d0
!
            else if (varv .eq. 1) then
!
! -------- DILATANCE
!---------- ELASTICITE EN DESSUS DU CRITERE VISQUEUX MAX
!
                dxip = dgamv
                variTmp(5) = 1.0d0
!
            end if
!
!---- XIP A T + DT ------------------------------------------------
!
            variTmp(1) = vinm(1)+dxip
!
!---- GAMMAP A T + DT ---------------------------------------------
!
            variTmp(2) = vinm(2)
!
! --  INDICATEUR DE PLASTICITE
            variTmp(7) = 0.d0
!
        else
! =================================================================
! -------- PLASTIFICATION -----------------------------------------
! =================================================================
            if (vinm(1) .lt. xipic) then
                val = 0
            else
                val = 1
            end if
!
! ------- CALCUL DE  GAMMAP -------------CRITERE ELASTOPLASTIQUE--
!
            call lkgamp(val, varv, i1ml, sml, ucrpm, &
                        seupm, vinm, nbmat, materd, de, &
                        depsth, depsv, dgamv, depsp, dgamp, &
                        iret)
!
            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if
!
! -------- DELTA XIP
!
            if (varv .eq. 0) then
!
!--------- CONTRACTANCE
!--------- PLASTIFICATION ET EN DESSOUS DU CRITERE VISQUEUX MAX
!
                dxip = dgamp
                variTmp(5) = 0.0d0
!
            else if (varv .eq. 1) then
!
! -------- DILATANCE
!--------- PLASTIFICATION ET EN DESSUS DU CRITERE VISQUEUX MAX
!
                dxip = dgamp+dgamv
                variTmp(5) = 1.0d0
!
            end if
! =================================================================
! --- REACTUALISATION DES CONTRAINTES  ----------------------------
! =================================================================
! --- DEFORMATIONS IRREVERSIBLES ----------------------------------
!
            irrev(1:ndt) = depsv(1:ndt)+depsp(1:ndt)
!
            vecd(1:ndt) = depsth(1:ndt)-irrev(1:ndt)
!
            dsig(1:ndt) = matmul(de(1:ndt, 1:ndt), vecd(1:ndt))
!
            do i = 1, ndt
                sigpl(i) = sigml(i)+dsig(i)
            end do
!
!==================================================================
!--------- REACTUALISATION DES VARIABLES INTERNES PLASTIQUES ------
!==================================================================
!---- XIP A T + DT ------------------------------------------------
!
            variTmp(1) = vinm(1)+dxip
!
!---- GAMMAP A T + DT ---------------------------------------------
!
            variTmp(2) = vinm(2)+dgamp
!
! --  INDICATEUR DE PLASTICITE
!
            variTmp(7) = 1.d0
        end if
    end if
!
! =================================================================
! --- TERMES DE L OPERATEUR TANGENT -------------------------------
! =================================================================
    if (option(11:14) .eq. 'ELAS') then
        call lkelas(ndi, ndt, nbmat, materd, depsth, &
                    sigml, de, kk, mu)
        dside(1:ndt, 1:ndt) = de(1:ndt, 1:ndt)
    end if
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
            if ((vinm(7) .eq. 0.d0) .and. (vinm(6) .eq. 0.d0)) then
                matr = 0
            else if ((vinm(7) .eq. 1.d0) .or. (vinm(6) .eq. 1.d0)) then
                matr = 1
            end if
        end if
        if (option(1:9) .eq. 'FULL_MECA') then
            if ((variTmp(7) .eq. 0.d0) .and. (variTmp(6) .eq. 0.d0)) then
                matr = 0
            else if ((variTmp(7) .eq. 1.d0) .or. (variTmp(6) .eq. 1.d0)) then
                matr = 1
            end if
        end if
        call r8inir(6*6, 0.d0, dside, 1)
        call lkelas(ndi, ndt, nbmat, materd, depsth, &
                    sigml, de, kk, mu)
!
        if (matr .eq. 0) then
!
!
            do i = 1, ndt
                do k = 1, ndt
                    dside(i, k) = de(i, k)
                end do
            end do
!
        else
!
            if (vinm(1) .lt. xipic) then
                val = 0
            else
                val = 1
            end if
!
            if (seuivm .lt. zero) then
                varv = 0
            else
                varv = 1
            end if
!
            vintr = vinm(3)
!
            call lkcrip(i1ml, sml, vinm, nbmat, materd, &
                        ucrpm, seupm)
!
            call lkcriv(vintr, i1ml, sml, vinm, nbmat, &
                        materd, ucrvm, seuvm)
!
            call lkcriv(vintr, iel, sel, vinm, nbmat, &
                        materd, ucriv, seuilv)
!
            call lkoptg(val, varv, dt, nbmat, materd, &
                        i1ml, sml, iel, sel, ucrpm, &
                        ucrvm, ucriv, seuilv, vinm, de, &
                        depsv, dside, iret)
!
!
            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if
!
        end if
!
    end if
!==================================================================
!--------- CONTRAINTES DE SORTIE:
! -------- RETABLISSEMENT DES SIGNES POUR ASTER --
!==================================================================
    do i = 1, ndt
        sigp(i) = mun*sigpl(i)
        deps(i) = mun*depsth(i)
    end do
! =================================================================
999 continue
! - Copy internal state variables
    if (lVari) then
        vinp = variTmp
    end if
end subroutine
