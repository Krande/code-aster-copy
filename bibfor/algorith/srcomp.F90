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

subroutine srcomp(mod, imate, instam, instap, &
                  deps, sigm, vinm, &
                  option, sigp, vinp, dside, retcom, &
                  invi)

!

!!!
!!! MODELE LKR : ROUTINE PRINCIPALE D'INTEGRATION EXPLICITE
!!!

! ===================================================================================
! IN  : NDIM    : DIMENSION DE L'ESPACE
!     : MOD     : TYPE DE MODELISATION
!     : IMATE   : ADRESSE DU MATERIAU CODE
!     : COMPOR  : COMPORTEMENT
!     : CRIT    : CRITERES DE CONVERGENCE LOCAUX
!     : INSTAM  : INSTANT DU CALCUL PRECEDENT
!     : INSTAP  : INSTANT DU CALCUL
!     : TM      : TEMPERATURE A L'INSTANT PRECEDENT
!     : TP      : TEMPERATURE A L'INSTANT DU CALCUL
!     : TREF    : TEMPERATURE DE REFERENCE
!     : DEPS    : INCREMENT DE DEFORMATION
!     : SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
!     : VINM    : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
!     : OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT : SIGP    : CONTRAINTES A L'INSTANT ACTUEL
!     : VINP    : VARIABLES INTERNES A L'INSTANT ACTUEL
!     : DSIDE   : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!     : RETCOM  : CODE RETOUR POUR LE REDECOUPAGE DU PAS DE TEMPS
!                 ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!                 L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!====================================================================================
! !!! ATTENTION !!! : CHANGEMENT DE SIGNES DES CHAMPS DE CONTRAINTES ET DES
!                     DEFORMATIONS - DANS CE MODELE CONVENTION MECANIQUE DES
!                     SOLS A L OPPPOSE DE CELLES DE LA MECANIQUE DES MILIEUX
!                     CONTINUS - EN COMPRESSION LA CONTRAINTE EST POSITIVE
!                     ET EN CONTRACTANCE : DEFORMATION VOLUMIQUE POSITIVE
!====================================================================================

    implicit none

#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/lcdevi.h"
#include "asterfort/srcrip.h"
#include "asterfort/srcriv.h"
#include "asterfort/srdgde.h"
#include "asterfort/srelas.h"
#include "asterfort/srgamp.h"
#include "asterfort/srlmat.h"
#include "asterfort/sroptg.h"
#include "asterfort/r8inir.h"
#include "asterfort/trace.h"
#include "asterfort/utmess.h"
#include "asterfort/get_varc.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: retcom, imate, invi
    character(len=8) :: mod(*)
    character(len=16) :: option
    real(kind=8) :: instam, instap, tm, tp, tref
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6), vinm(invi)
    real(kind=8) :: sigp(6), vinp(invi)
    real(kind=8) :: dside(6, 6)
    aster_logical :: l_temp

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: nbmat, ndt, ndi, nvi, val, varv, i, k, matr, iret, indal

    parameter(nbmat=90)

    real(kind=8) :: mun, un, zero, deux, trois, materd(nbmat, 2), materf(nbmat, 2)
    real(kind=8) :: dt, alpha, coef, coupl, sigml(6), sigpl(6), depml(6), depsth(6)
    real(kind=8) :: i1ml, sml(6), siim, iel, i1el, sel1(6), dvml, devml(6)
    real(kind=8) :: dvml1, devml1(6), sel(6), sigel(6), seuilv, seuilp
    real(kind=8) :: ucrvm, seuvm, ucrpm, seupm, depsv(6), dgamv, dxi5, xi1
    real(kind=8) :: depsp(6), dgamp, xi5, dxip, dxiv, seuip2, ucrip2, ucrip, ucriv
    real(kind=8) :: dsig(6), vecd(6), irrev(6), de(6, 6), kk, mu, kron(6), vintr
    real(kind=8) :: xi10, xi50, rx1, rx5, tmm, tpp, trr, somme, dtempm, dtempp, dtemp
    real(kind=8) :: xi20, xi2, rx2

    character(len=3) :: matcst

    common/tdim/ndt, ndi

    parameter(mun=-1.d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    data kron/un, un, un, zero, zero, zero/

    dt = instap-instam
    retcom = 0
    dgamp = zero
    dgamv = zero
    dxip = zero
    dxiv = zero
    seuip2 = zero
    call r8inir(6, zero, depsp, 1)
    call r8inir(6, zero, depsv, 1)
!
! - Get temperatures
!
    call get_varc('RIGI', 1, 1, 'T', &
                  tm, tp, tref, l_temp)
!
! - Glute pour LKR
!
    if (.not. l_temp) then
        tm = 0.d0
        tp = 0.d0
        tref = 0.d0
    end if

    !!!
    !!! Recuperation des parametres du modele : les coefficients materiau n'evoluent
    !!! pas avec le temps ; les temperatures tm, tp et tref sont stockees dans les
    !!! tableaux materd et materf et les increments (tm-tref), (tp-tref) et (tp-tm)
    !!! aussi
    !!!

    call srlmat(mod(1), imate, nbmat, tm, tp, tref, materd, materf, matcst, ndt, ndi, nvi, indal)

    ASSERT(invi .eq. nvi)

    !!!
    !!! T-   : materd(6,1)
    !!! T+   : materd(7,1)
    !!! Tref : materd(8,1)
    !!! dT-  : materd(9,1)
    !!! dT+  : materd(10,1)
    !!! dT   : materd(11,1)
    !!!

    tmm = materd(6, 1)
    tpp = materd(7, 1)
    trr = materd(8, 1)
    dtempm = materd(9, 1)
    dtempp = materd(10, 1)
    dtemp = materd(11, 1)

    !!! Definition de xi_1, xi_2 et xi_5 a t+dt

    xi10 = materd(12, 2)
    xi20 = materd(13, 2)
    xi50 = materd(14, 2)
    rx1 = materd(24, 2)
    rx2 = materd(25, 2)
    rx5 = materd(26, 2)
    xi1 = xi10*exp(rx1*dtempp)
    xi2 = xi20*exp(rx2*dtempp)
    xi5 = xi50*exp(rx5*dtempp)

    !!! Couplage entre les deux mecanismes

    coupl = materd(28, 2)

    !!! Changement de signe des contraintes et deformations

    do i = 1, ndt
        sigml(i) = mun*sigm(i)
        depml(i) = mun*deps(i)
    end do

    !!! Definition de i1, s et sii à t-

    i1ml = trace(ndi, sigml)
    call lcdevi(sigml, sml)
    siim = norm2(sml(1:ndt))

    !!!
    !!! Prise en compte de la dilatation thermique
    !!!

    alpha = materd(3, 1)
    coef = alpha*(dtempp-dtempm)

    !!!
    !!! Definition des deformations volumiques et deviatoriques
    !!!

    do k = 1, ndt
        depsth(k) = depml(k)
    end do

    dvml = zero
    do k = 1, ndi
        depsth(k) = depsth(k)+coef
        dvml = dvml+depsth(k)
    end do

    do k = 1, ndt
        devml(k) = depsth(k)-dvml*kron(k)/trois
    end do

    !!! actualisation de eps_v^te --vin(9)
    if (option .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        vinp(9) = vinm(9)-trois*coef
    end if

    !!!
    !!! Verification d'un etat initial plastiquement admissible
    !!!

    somme = zero
    do i = 1, nvi
        somme = somme+vinm(i)
    end do
    if (abs(somme) .lt. r8prem()) then
        call srcrip(i1ml, sml, vinm, nvi, nbmat, materd, tmm, ucrpm, seupm)
        if (seupm/materd(4, 1) .gt. 1.0d-6) then
            call utmess('F', 'ALGORITH2_81')
        end if
    end if

    !!!
    !!! Prediction elastique
    !!!

    call srelas(ndi, ndt, nbmat, materd, sigml, de, kk, mu)

    iel = i1ml+trois*kk*dvml

    do i = 1, ndt
        sel(i) = sml(i)+deux*mu*devml(i)
    end do

    do i = 1, ndt
        sigel(i) = sel(i)+iel/trois*kron(i)
    end do

    !!!
    !!! Test sur les criteres visco. et plast.
    !!!

    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        !!! criter visco.
        vintr = vinm(3)

        call srcriv(vintr, iel, sel, nbmat, materd, tpp, ucriv, seuilv)
        call srcriv(vintr, i1ml, sml, nbmat, materd, tmm, ucrvm, seuvm)

        if (seuilv .lt. zero) then
            !!! elasticite -- maj des variables internes
            val = 0
            dgamv = zero
            dxiv = zero
            dvml1 = zero
            do i = 1, ndt
                depsv(i) = zero
                devml1(i) = zero
            end do
            vinp(3) = vinm(3)
            vinp(4) = vinm(4)
            vinp(6) = 0
            vinp(11) = vinm(11)
        else
            !!! viscoplasticite
            val = 0

            !!! calcul de depsv et dgamv

            call srdgde(val, vintr, dt, seuilv, ucrvm, &
                        i1ml, sml, vinm, nvi, nbmat, materd, &
                        tmm, depsv, dgamv, iret)
            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if
            dvml1 = trace(ndi, depsv)
            call lcdevi(depsv, devml1)

            !!! maj des variables internes
            dxi5 = xi5-vinm(3)
            dxiv = min(dgamv, dxi5)
            vinp(3) = vinm(3)+dxiv
            vinp(4) = vinm(4)+dgamv
            vinp(6) = 1
            vinp(11) = vinm(11)+dvml1
        end if

        !!! pas de maj de la prediction elastique -- difference / letk
        !!! avec maj : i1el = iel - trois*kk*dvml1
        !!! sel1(i) = sel(i) - deux*mu*devml1(i)

        i1el = iel-trois*kk*dvml1
        do i = 1, ndt
            sel1(i) = sel(i)-deux*mu*devml1(i)
        end do

        !!!
        !!! Critere plastique
        !!!

        call srcrip(i1ml, sml, vinm, nvi, nbmat, materd, tmm, ucrpm, seupm)

        !!! verification couplage entre les deux mecanismes
        call srcriv(xi5, i1el, sel1, nbmat, materd, tpp, ucrip2, seuip2)

        if (seuip2 .lt. zero) then
            !!! contractance
            varv = 0
            vinp(5) = zero
        else
            !!! dilatance
            varv = 1
            vinp(5) = un
        end if

        !!! calcul de fp
        call srcrip(i1el, sel1, vinm, nvi, nbmat, materd, tpp, ucrip, seuilp)
        if ((ucrip .lt. zero) .or. (ucrpm .lt. zero)) then
            retcom = 1
            goto 999
        end if

        if (seuilp .lt. zero) then
            !!! elasticite
            dgamp = zero
            do i = 1, ndt
                depsp(i) = zero
            end do

            !!! maj des contraintes
            irrev(1:ndt) = depsv(1:ndt)+depsp(1:ndt)
            vecd(1:ndt) = depsth(1:ndt)-irrev(1:ndt)
            dsig(1:ndt) = matmul(de(1:ndt, 1:ndt), vecd(1:ndt))

            do i = 1, ndt
                sigpl(i) = sigml(i)+dsig(i)
            end do

            !!! maj des variables internes
            if ((varv .eq. 1) .and. (coupl .ge. un/deux)) then
                vinp(1) = vinm(1)+dgamv
            else
                vinp(1) = vinm(1)
            end if

            vinp(2) = vinm(2)
            vinp(7) = 0
            vinp(10) = vinm(10)
            vinp(8) = vinm(8)+depml(1)+depml(2)+depml(3)+trois*coef- &
                      (vinp(10)-vinm(10))-(vinp(11)-vinm(11))

        else

            !!! plasticite
            if (vinm(1) .lt. xi1) then
                !!! loi de dilatance pre-pic
                val = 0
            else
                !!! loi de dilatance post-pic
                val = 1
            end if

            !!! calcul de gammap
            call srgamp(val, varv, i1ml, sml, ucrpm, &
                        seupm, vinm, nvi, nbmat, materd, de, &
                        depsth, depsv, dgamv, depsp, dgamp, iret)

            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if

            !!!maj des contraintes
            irrev(1:ndt) = depsv(1:ndt)+depsp(1:ndt)
            vecd(1:ndt) = depsth(1:ndt)-irrev(1:ndt)
            dsig(1:ndt) = matmul(de(1:ndt, 1:ndt), vecd(1:ndt))

            do i = 1, ndt
                sigpl(i) = sigml(i)+dsig(i)
            end do

            !!! maj des variables internes
            if ((varv .eq. 1) .and. (coupl .ge. un/deux)) then
                vinp(1) = vinm(1)+dgamp+dgamv
            else
                vinp(1) = vinm(1)+dgamp
            end if

            vinp(2) = vinm(2)+dgamp
            vinp(7) = 1
            vinp(10) = vinm(10)+depsp(1)+depsp(2)+depsp(3)
            vinp(8) = vinm(8)+depml(1)+depml(2)+depml(3)+trois*coef- &
                      (vinp(10)-vinm(10))-(vinp(11)-vinm(11))

        end if

        !!! post-traitement domaine
        if (vinp(1) .le. zero) then
            vinp(12) = 0
        else if ((vinp(1) .gt. zero) .and. (vinp(1) .lt. xi1)) then
            vinp(12) = 1
        else if ((vinp(1) .ge. xi1) .and. (vinp(1) .lt. xi2)) then
            vinp(12) = 2
        else if (vinp(1) .ge. xi2) then
            vinp(12) = 3
        end if

    end if

!!!
!!! Operateur tangent
!!!

    if (option(11:14) .eq. 'ELAS') then

        call srelas(ndi, ndt, nbmat, materd, sigml, de, kk, mu)
        dside(1:ndt, 1:ndt) = de(1:ndt, 1:ndt)

    end if

    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then

        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
            if ((vinm(7) .le. 0) .and. (vinm(6) .le. 0)) then
                matr = 0
            else if ((vinm(7) .gt. 0) .or. (vinm(6) .gt. 0)) then
                matr = 1
            end if
        end if

        if (option(1:9) .eq. 'FULL_MECA') then
            if ((vinp(7) .le. 0) .and. (vinp(6) .le. 0)) then
                matr = 0
            else if ((vinp(7) .gt. 0) .or. (vinp(6) .gt. 0)) then
                matr = 1
            end if
        end if

        call r8inir(6*6, 0.d0, dside, 1)
        call srelas(ndi, ndt, nbmat, materd, sigml, de, kk, mu)

        if (matr .eq. 0) then

            do i = 1, ndt
                do k = 1, ndt
                    dside(i, k) = de(i, k)
                end do
            end do

        else

            if (vinm(1) .lt. xi1) then
                val = 0
            else
                val = 1
            end if

            if (seuip2 .lt. zero) then
                varv = 0
            else
                varv = 1
            end if

            vintr = vinm(3)

            call srcrip(i1ml, sml, vinm, nvi, nbmat, materd, tmm, ucrpm, seupm)
            call srcriv(vintr, i1ml, sml, nbmat, materd, tmm, ucrvm, seuvm)
            call srcriv(vintr, iel, sel, nbmat, materd, tpp, ucriv, seuilv)
            call sroptg(val, varv, dt, nbmat, materd, i1ml, sml, sel, ucrpm, &
                        ucrvm, ucriv, seuilv, vinm, nvi, de, depsv, dside, iret)

            if (iret .eq. 1) then
                retcom = 1
                goto 999
            end if

        end if

    end if

!!!
!!! retablissement des contraintes pour aster
!!!
    if (option .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        do i = 1, ndt
            sigp(i) = mun*sigpl(i)
        end do
    end if

    do i = 1, ndt
        deps(i) = mun*depsth(i)
    end do

999 continue

end subroutine
