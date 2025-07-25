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
subroutine nzcizi(fami, kpg, ksp, ndim, imat, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, vim, option, sigp, &
                  vip, dsidep, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nzcalc.h"
#include "asterfort/verift.h"
#include "asterfort/metaGetMechanism.h"
#include "asterfort/metaGetType.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetParaVisc.h"
#include "asterfort/metaGetParaHardLine.h"
#include "asterfort/metaGetParaMixture.h"
#include "asterfort/metaGetParaPlasTransf.h"
#include "asterfort/metaGetParaAnneal.h"
#include "asterfort/metaGetParaElas.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imat
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    real(kind=8), intent(in) :: epsm(*)
    real(kind=8), intent(in) :: deps(*)
    real(kind=8), intent(in) :: sigm(*)
    real(kind=8), intent(in) :: vim(25)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(*)
    real(kind=8), intent(out) :: vip(25)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Comportment
!
! META_P_C* / META_V_C* for small strains and zirconium metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMAT    : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE
!     IRET    : CODE RETOUR DE LA RESOLUTION DE L'EQUATION SCALAIRE
!               (NZCALC)
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX YY ZZ XY XZ YZ
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_phase, meta_type
    integer(kind=8) :: ndimsi, i, j, k, l, mode
    real(kind=8) :: phase(5), phasm(5), zalpha, deltaz(5)
    real(kind=8) :: dt, coef_hard
    real(kind=8) :: epsth, e, deuxmu, deumum, troisk
    real(kind=8) :: fmel, sy(3), symoy, h(3), hmoy, rprim
    real(kind=8) :: theta(4)
    real(kind=8) :: eta(5), n(3), unsurn(5), c(3), m(3), cmoy, mmoy, cr
    real(kind=8) :: dz(2), dz1(2), dz2(2), vi(18), dvin, vimt(18)
    real(kind=8) :: xmoy(6), ds(6), xmoyeq
    real(kind=8) :: trans, kpt(2), fpt(2)
    real(kind=8) :: trepsm, trdeps, trsigm, trsigp
    real(kind=8) :: dvdeps(6), dvsigm(6), dvsigp(6)
    real(kind=8) :: sigel(6), sigel2(6), sig0(6), sieleq, sigeps
    real(kind=8) :: plasti, dp, seuil
    real(kind=8) :: coef1, coef2, coef3, dv, n0(3), b
    character(len=1) :: poum
    aster_logical :: resi, rigi
    aster_logical :: l_visc, l_plas, l_anneal, l_plas_tran
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    character(len=16) :: metaRela, metaGlob
!
! --------------------------------------------------------------------------------------------------
!

    iret = 0
    ndimsi = 2*ndim
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    dt = instap-instam

! - Behaviour in kit
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)

!
! - Get metallurgy type
!
    call metaGetType(meta_type, nb_phase)
    ASSERT(meta_type .eq. META_ZIRC)
    ASSERT(nb_phase .eq. 3)
!
! - Get phasis
!
    if (resi) then
        poum = '+'
        call metaGetPhase(fami, '+', kpg, ksp, meta_type, &
                          nb_phase, phase, zcold_=zalpha)
        call metaGetPhase(fami, '-', kpg, ksp, meta_type, &
                          nb_phase, phasm)
    else
        poum = '-'
        call metaGetPhase(fami, '-', kpg, ksp, meta_type, &
                          nb_phase, phase, zcold_=zalpha)
    end if
    do k = 1, nb_phase-1
        deltaz(k) = phase(k)-phasm(k)
    end do
!
! - Compute thermic strain
!
    call verift(fami, kpg, ksp, poum, imat, &
                epsth_meta_=epsth)
!
! - Mechanisms of comportment law
!
    call metaGetMechanism(metaRela, metaGlob, &
                          l_plas=l_plas, &
                          l_visc=l_visc, &
                          l_anneal=l_anneal, &
                          l_plas_tran=l_plas_tran)
!
! - Get elastic parameters
!
    call metaGetParaElas(poum, fami, kpg, ksp, imat, &
                         e_=e, deuxmu_=deuxmu, troisk_=troisk, &
                         deuxmum_=deumum)
    plasti = vim(25)
!
! - Mixture law (yield limit)
!
    call metaGetParaMixture(poum, fami, kpg, ksp, imat, &
                            l_visc, meta_type, nb_phase, zalpha, fmel, &
                            sy)
!
! - Get hardening slope (linear)
!
    coef_hard = (2.d0/3.d0)
    call metaGetParaHardLine(poum, fami, kpg, ksp, imat, &
                             meta_type, nb_phase, &
                             e, coef_hard, h)
    hmoy = 0.d0
    do k = 1, nb_phase
        hmoy = hmoy+phase(k)*h(k)
    end do
!
    if (resi) then
        sigp(1:2*ndim) = 0.d0
        vip(1:25) = 0.d0
! ----- Parameters for annealing
        if (l_anneal) then
            call metaGetParaAnneal(poum, fami, kpg, ksp, imat, &
                                   meta_type, nb_phase, &
                                   theta)
        else
            do i = 1, 4
                theta(i) = 1.d0
            end do
        end if
! ----- Parameters for viscosity
        if (l_visc) then
            call metaGetParaVisc(poum, fami, kpg, ksp, imat, &
                                 meta_type, nb_phase, eta, n, unsurn, &
                                 c, m)
        else
            eta(:) = 0.d0
            n(:) = 20.d0
            unsurn(:) = 1.d0
            c(:) = 0.d0
            m(:) = 20.d0
        end if
!
! 2.7 - CALCUL DE VIM+DG
!
        do k = 1, nb_phase-1
            dz(k) = phase(k)-phasm(k)
            if (dz(k) .ge. 0.d0) then
                dz1(k) = dz(k)
                dz2(k) = 0.d0
            else
                dz1(k) = 0.d0
                dz2(k) = -dz(k)
            end if
        end do
        if (phase(nb_phase) .gt. 0.d0) then
            do i = 1, ndimsi
                dvin = 0.d0
                do k = 1, nb_phase-1
                    l = i+(k-1)*6
                    dvin = dvin+dz2(k)*(theta(2+k)*vim(l)-vim(12+i))/phase(nb_phase)
                end do
                vi(12+i) = vim(12+i)+dvin
                if ((vi(12+i)*vim(12+i)) .lt. 0.d0) vi(12+i) = 0.d0
            end do
        else
            do i = 1, ndimsi
                vi(12+i) = 0.d0
            end do
        end if
        do k = 1, nb_phase-1
            do i = 1, ndimsi
                l = i+(k-1)*6
                if (phase(k) .gt. 0.d0) then
                    dvin = dz1(k)*(theta(k)*vim(12+i)-vim(l))/phase(k)
                    vi(l) = vim(l)+dvin
                    if ((vi(l)*vim(l)) .lt. 0.d0) vi(l) = 0.d0
                else
                    vi(l) = 0.d0
                end if
            end do
        end do
!
!    - MISE AU FORMAT DES CONTRAINTES DE RAPPEL
!
        do i = 4, ndimsi
            do k = 1, nb_phase
                l = i+(k-1)*6
                vi(l) = vi(l)*rac2
            end do
        end do
!
! 2.8 - RESTAURATION D ORIGINE VISQUEUSE
!
        do i = 1, ndimsi
            xmoy(i) = 0.d0
            do k = 1, nb_phase
                l = i+(k-1)*6
                xmoy(i) = xmoy(i)+phase(k)*h(k)*vi(l)
            end do
        end do
        xmoyeq = 0.d0
        do i = 1, ndimsi
            xmoyeq = xmoyeq+xmoy(i)**2.d0
        end do
        xmoyeq = sqrt(1.5d0*xmoyeq)
        cmoy = 0.d0
        mmoy = 0.d0
        do k = 1, nb_phase
            cmoy = cmoy+phase(k)*c(k)
            mmoy = mmoy+phase(k)*m(k)
        end do
        cr = cmoy*xmoyeq
        if (xmoyeq .gt. 0.d0) then
            do i = 1, ndimsi
                ds(i) = 3.d0*dt*(cr**mmoy)*xmoy(i)/(2.d0*xmoyeq)
            end do
        else
            do i = 1, ndimsi
                ds(i) = 0.d0
            end do
        end if
        do k = 1, nb_phase
            do i = 1, ndimsi
                l = i+(k-1)*6
                if (phase(k) .gt. 0.d0) then
                    vimt(l) = vi(l)
                    vi(l) = vi(l)-ds(i)
                    if ((vi(l)*vimt(l)) .lt. 0.d0) vi(l) = 0.d0
                end if
            end do
        end do
! ----- Parameters for plasticity of tranformation
        trans = 0.d0
        if (l_plas_tran) then
            call metaGetParaPlasTransf('+', fami, 1, 1, imat, &
                                       meta_type, nb_phase, deltaz, zalpha, &
                                       kpt, fpt)
            do k = 1, nb_phase-1
                if (deltaz(k) .gt. 0.d0) then
                    trans = trans+kpt(k)*fpt(k)*deltaz(k)
                end if
            end do
        end if
    else
        do k = 1, nb_phase
            do i = 1, ndimsi
                l = i+(k-1)*6
                vi(l) = vim(l)
                if (i .gt. 3) then
                    vi(l) = vi(l)*rac2
                end if
            end do
        end do
        trans = 0.d0
        do i = 1, ndimsi
            xmoy(i) = 0.d0
            do k = 1, nb_phase
                l = i+(k-1)*6
                xmoy(i) = xmoy(i)+phase(k)*h(k)*vi(l)
            end do
        end do
    end if
!
! 2.10 - CALCUL DE SYMOY
!
    if (zalpha .gt. 0.d0) then
        symoy = phase(1)*sy(1)+phase(2)*sy(2)
        symoy = symoy/zalpha
    else
        symoy = 0.d0
    end if
    symoy = (1.d0-fmel)*sy(nb_phase)+fmel*symoy
!
! ********************************
! 3 - DEBUT DE L ALGORITHME
! ********************************
!
    trdeps = (deps(1)+deps(2)+deps(3))/3.d0
    trepsm = (epsm(1)+epsm(2)+epsm(3))/3.d0
    trsigm = (sigm(1)+sigm(2)+sigm(3))/3.d0
    trsigp = troisk*(trepsm+trdeps)-troisk*epsth
    do i = 1, ndimsi
        dvdeps(i) = deps(i)-trdeps*kron(i)
        dvsigm(i) = sigm(i)-trsigm*kron(i)
    end do
    sieleq = 0.d0
    do i = 1, ndimsi
        sigel(i) = deuxmu*dvsigm(i)/deumum+deuxmu*dvdeps(i)
        sigel2(i) = sigel(i)-(1.5d0*deuxmu*trans+1.d0)*xmoy(i)
        sieleq = sieleq+sigel2(i)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    if (sieleq .gt. 0.d0) then
        do i = 1, ndimsi
            sig0(i) = sigel2(i)/sieleq
        end do
    else
        do i = 1, ndimsi
            sig0(i) = 0.d0
        end do
    end if
!
! ************************
! 4 - RESOLUTION
! ************************
!
    if (resi) then
!
! 4.2.1 - CALCUL DE DP
!
        seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*symoy
        if (seuil .lt. 0.d0) then
            vip(25) = 0.d0
            dp = 0.d0
        else
            vip(25) = 1.d0
            rprim = 3.d0*hmoy/2.d0
            if (l_plas) then
                dp = seuil/(1.5d0*deuxmu+(1.5d0*deuxmu*trans+1.d0)* &
                            rprim)
            else
                call nzcalc(carcri, nb_phase, phase, zalpha, &
                            fmel, seuil, dt, trans, &
                            rprim, deuxmu, eta, unsurn, &
                            dp, iret)
                if (iret .eq. 1) goto 999
            end if
        end if
!
! 4.2.2 - CALCUL DE SIGMA
!
        plasti = vip(25)
        do i = 1, ndimsi
            dvsigp(i) = sigel(i)-1.5d0*deuxmu*dp*sig0(i)
            dvsigp(i) = dvsigp(i)/(1.5d0*deuxmu*trans+1.d0)
            sigp(i) = dvsigp(i)+trsigp*kron(i)
        end do
!
! 4.2.3 - CALCUL DE VIP ET XMOY
!
        do k = 1, nb_phase
            do i = 1, ndimsi
                l = i+(k-1)*6
                if (phase(k) .gt. 0.d0) then
                    vip(l) = vi(l)+3.d0*dp*sig0(i)/2.d0
                    if (i .gt. 3) then
                        vip(l) = vip(l)/rac2
                    end if
                else
                    vip(l) = 0.d0
                end if
            end do
        end do
        do i = 1, ndimsi
            vip(18+i) = xmoy(i)+3.d0*hmoy*dp*sig0(i)/2.d0
            if (i .gt. 3) then
                vip(18+i) = vip(18+i)/rac2
            end if
        end do
    end if
!
! *******************************
! 5 - MATRICE TANGENTE DSIGDF
! *******************************
!
    if (rigi) then
        mode = 2
        if (l_visc) mode = 1
        dsidep(1:ndimsi, 1:ndimsi) = 0.d0
        do i = 1, ndimsi
            dsidep(i, i) = 1.d0
        end do
        do i = 1, 3
            do j = 1, 3
                dsidep(i, j) = dsidep(i, j)-1.d0/3.d0
            end do
        end do
        if (option(1:9) .eq. 'FULL_MECA') then
            coef1 = (1.5d0*deuxmu*trans+1.d0)
        else
            coef1 = 1.d0
        end if
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*deuxmu/coef1
            end do
        end do
!
! 5.2 - PARTIE PLASTIQUE
!
        b = 1.d0
        coef2 = 0.d0
        coef3 = 0.d0

        if (plasti .ge. 0.5d0) then
            if (option(1:9) .eq. 'FULL_MECA') then
                sigeps = 0.d0
                do i = 1, ndimsi
                    dvsigp(i) = dvsigp(i)-xmoy(i)
                    sigeps = sigeps+dvsigp(i)*dvdeps(i)
                end do
                if ((mode .eq. 1) .or. ((mode .eq. 2) .and. (sigeps .ge. 0.d0))) then
                    b = 1.d0-(1.5d0*deuxmu*dp/sieleq)
                    dv = 0.d0
                    if (mode .eq. 1) then
                        do k = 1, nb_phase
                            n0(k) = (1-n(k))/n(k)
                        end do
                        dv = (1-fmel)*phase(nb_phase)*(eta(nb_phase)/n(nb_phase)/dt)* &
                             ((dp/dt)**n0(nb_phase))
                        if (zalpha .gt. 0.d0) then
                            do k = 1, nb_phase-1
                                if (phase(k) .gt. 0.d0) dv = dv+fmel*(phase(k)/zalpha)*&
                                                             & (eta(k)/n(k)/dt)*((dp/dt)**n0&
                                                             &(k))
                            end do
                        end if
                    end if
                    coef2 = 3.d0*hmoy/2.d0+dv
                    coef2 = (1.5d0*deuxmu*trans+1.d0)*coef2
                    coef2 = (1.5d0*deuxmu)+coef2
                    coef2 = 1/coef2-dp/sieleq
                    coef2 = ((1.5d0*deuxmu)**2)*coef2
                end if
            end if
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                if (mode .eq. 2) coef2 = ((1.5d0*deuxmu)**2)/(1.5d0*deuxmu+1.5d0*hmoy)
            end if
            coef3 = coef2/coef1
        end if
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*b
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                dsidep(i, j) = dsidep(i, j)+troisk/3.d0
            end do
        end do
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)-coef3*sig0(i)*sig0(j)
            end do
        end do
    end if
!
999 continue
!
end subroutine
