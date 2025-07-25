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
subroutine nzisfw(fami, kpg, ksp, ndim, imat, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, vim, option, sigp, &
                  vip, dsidep, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/nzcalc.h"
#include "asterfort/rcvarc.h"
#include "asterfort/verift.h"
#include "asterfort/metaGetMechanism.h"
#include "asterfort/metaGetType.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetParaVisc.h"
#include "asterfort/metaGetParaHardLine.h"
#include "asterfort/metaGetParaHardTrac.h"
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
    real(kind=8), intent(in) :: vim(*)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(*)
    real(kind=8), intent(out) :: vip(*)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Comportment
!
! META_P_I* / META_V_I* for small strains and steel metallurgy
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
    integer(kind=8) :: maxval, nb_phase, meta_type
    integer(kind=8) :: ndimsi, i, j, k, mode, iret2
    real(kind=8) :: phase(5), phasm(5), zalpha, deltaz(5)
    real(kind=8) :: temp, dt, coef_hard
    real(kind=8) :: epsth, e, deuxmu, deumum, troisk
    real(kind=8) :: fmel, sy(5), h(5), hmoy, hplus(5), r(5), rmoy
    real(kind=8) :: theta(8)
    real(kind=8) :: eta(5), n(5), unsurn(5), c(5), m(5), cmoy, mmoy, cr
    real(kind=8) :: dz(4), dz1(4), dz2(4), vi(5), dvin, vimoy, ds
    real(kind=8) :: trans, kpt(4), fpt(4)
    real(kind=8) :: trepsm, trdeps, trsigm, trsigp
    real(kind=8) :: dvdeps(6), dvsigm(6), dvsigp(6)
    real(kind=8) :: sigel(6), sig0(6), sieleq, sigeps
    real(kind=8) :: plasti, dp, seuil
    real(kind=8) :: coef1, coef2, coef3, dv, n0(5), b
    real(kind=8) :: precr
    character(len=1) :: poum
    integer(kind=8) :: test
    aster_logical :: resi, rigi, l_temp
    aster_logical :: l_visc, l_plas, l_anneal, l_plas_tran, l_hard_isotline, l_hard_isotnlin
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    character(len=16) :: metaRela, metaGlob
!
! --------------------------------------------------------------------------------------------------
!
    ndimsi = 2*ndim
    iret = 0
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    dt = instap-instam
    precr = r8prem()

! - Behaviour in kit
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)
!
! - Get metallurgy type
!
    call metaGetType(meta_type, nb_phase)
    ASSERT(meta_type .eq. META_STEEL)
    ASSERT(nb_phase .eq. 5)
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
    call rcvarc(' ', 'TEMP', poum, fami, kpg, &
                ksp, temp, iret2)
    l_temp = iret2 .eq. 0
!
! - Mechanisms of comportment law
!
    call metaGetMechanism(metaRela, metaGlob, &
                          l_plas=l_plas, &
                          l_visc=l_visc, &
                          l_anneal=l_anneal, &
                          l_plas_tran=l_plas_tran, &
                          l_hard_isotline=l_hard_isotline, &
                          l_hard_isotnlin=l_hard_isotnlin)
!
! - Get elastic parameters
!
    call metaGetParaElas(poum, fami, kpg, ksp, imat, &
                         e_=e, deuxmu_=deuxmu, troisk_=troisk, &
                         deuxmum_=deumum)
    plasti = vim(IDX_I_IPLAS)
    trans = 0.d0
!
! - Mixture law (yield limit)
!
    call metaGetParaMixture(poum, fami, kpg, ksp, imat, &
                            l_visc, meta_type, nb_phase, zalpha, fmel, &
                            sy)
!
    if (resi) then
! ----- Parameters for annealing
        if (l_anneal) then
            call metaGetParaAnneal(poum, fami, kpg, ksp, imat, &
                                   meta_type, nb_phase, &
                                   theta)
        else
            do i = 1, 8
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
! 2.6 - CALCUL DE VIM+DG-DS ET DE RMOY
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
            dvin = 0.d0
            do k = 1, nb_phase-1
                dvin = dvin+dz2(k)*(theta(4+k)*vim(k)-vim(nb_phase))/phase(nb_phase)
            end do
            vi(nb_phase) = vim(nb_phase)+dvin
            vimoy = phase(nb_phase)*vi(nb_phase)
        else
            vi(nb_phase) = 0.d0
            vimoy = 0.d0
        end if
        do k = 1, nb_phase-1
            if (phase(k) .gt. 0.d0) then
                dvin = dz1(k)*(theta(k)*vim(nb_phase)-vim(k))/phase(k)
                vi(k) = vim(k)+dvin
                vimoy = vimoy+phase(k)*vi(k)
            else
                vi(k) = 0.d0
            end if
        end do
!
! 2.7 - RESTAURATION D ORIGINE VISQUEUSE
!
        cmoy = 0.d0
        mmoy = 0.d0
        do k = 1, nb_phase
            cmoy = cmoy+phase(k)*c(k)
            mmoy = mmoy+phase(k)*m(k)
        end do
        cr = cmoy*vimoy
        if (cr .le. 0.d0) then
            ds = 0.d0
        else
            ds = dt*(cr**mmoy)
        end if
        do k = 1, nb_phase
            if (phase(k) .gt. 0.d0) then
                vi(k) = vi(k)-ds
                if (vi(k) .le. 0.d0) vi(k) = 0.d0
            end if
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
            vi(k) = vim(k)
        end do
    end if
!
! 2.9 - CALCUL DE HMOY ET RMOY (ON INCLUE LE SIGY)
!
    if (l_hard_isotline) then
! ----- Get hardening slope (linear)
        coef_hard = (1.d0)
        call metaGetParaHardLine(poum, fami, kpg, ksp, imat, &
                                 meta_type, nb_phase, &
                                 e, coef_hard, h)
        do k = 1, nb_phase
            r(k) = h(k)*vi(k)+sy(k)
        end do
    end if
    if (l_hard_isotnlin) then
! ----- Get hardening slope (non-linear)
        call metaGetParaHardTrac(imat, meta_type, nb_phase, &
                                 l_temp, temp, &
                                 vi, h, r, maxval)
        do k = 1, nb_phase
            r(k) = r(k)+sy(k)
        end do
    end if
    if (zalpha .gt. 0.d0) then
        rmoy = phase(1)*r(1)+phase(2)*r(2)+phase(3)*r(3)+phase(4)*r(4)
        rmoy = rmoy/zalpha
        hmoy = phase(1)*h(1)+phase(2)*h(2)+phase(3)*h(3)+phase(4)*h(4)
        hmoy = hmoy/zalpha
    else
        rmoy = 0.d0
        hmoy = 0.d0
    end if
    rmoy = (1.d0-fmel)*r(nb_phase)+fmel*rmoy
    hmoy = (1.d0-fmel)*h(nb_phase)+fmel*hmoy
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
        sieleq = sieleq+sigel(i)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    if (sieleq .gt. 0.d0) then
        do i = 1, ndimsi
            sig0(i) = sigel(i)/sieleq
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
        vip(1:IDX_I_IPLAS) = 0.d0
        seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*rmoy
        if (seuil .lt. 0.d0) then
            vip(IDX_I_IPLAS) = 0.d0
            dp = 0.d0
        else
            vip(IDX_I_IPLAS) = 1.d0
            call nzcalc(carcri, nb_phase, phase, zalpha, &
                        fmel, seuil, dt, trans, &
                        hmoy, deuxmu, eta, unsurn, &
                        dp, iret)
            if (iret .eq. 1) goto 999
!
! DANS LE CAS NON LINEAIRE
! VERIFICATION QU ON EST DANS LE BON INTERVALLE
!
            if (l_hard_isotnlin) then
!
                do j = 1, maxval
                    test = 0
                    vip(1:nb_phase) = vi(1:nb_phase)+dp
                    hplus(1:nb_phase) = h(1:nb_phase)
                    call metaGetParaHardTrac(imat, meta_type, nb_phase, &
                                             l_temp, temp, &
                                             vip, h, r)
                    do k = 1, nb_phase
                        if (phase(k) .gt. 0.d0) then
                            r(k) = r(k)+sy(k)
                            if (abs(h(k)-hplus(k)) .gt. precr) test = 1
                        end if
                    end do
                    if (test .eq. 0) goto 600
                    hmoy = 0.d0
                    rmoy = 0.d0
                    if (zalpha .gt. 0.d0) then
                        do k = 1, nb_phase-1
                            if (phase(k) .gt. 0.d0) then
                                rmoy = rmoy+phase(k)*(r(k)-h(k)*dp)
                                hmoy = hmoy+phase(k)*h(k)
                            end if
                        end do
                        rmoy = fmel*rmoy/zalpha
                        hmoy = fmel*hmoy/zalpha
                    end if
                    if (phase(nb_phase) .gt. 0.d0) then
                        rmoy = (1.d0-fmel)*(r(nb_phase)-h(nb_phase)*dp)+rmoy
                        hmoy = (1.d0-fmel)*h(nb_phase)+hmoy
                    end if
                    seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*rmoy
                    call nzcalc(carcri, nb_phase, phase, zalpha, &
                                fmel, seuil, dt, trans, &
                                hmoy, deuxmu, eta, unsurn, &
                                dp, iret)
                    if (iret .eq. 1) goto 999
                end do
                ASSERT((test .ne. 1) .or. (j .ne. maxval))
600             continue
            end if
        end if
!
! 4.2.2 - CALCUL DE SIGMA
!
        plasti = vip(IDX_I_IPLAS)
        sigp(1:2*ndim) = 0.d0
        do i = 1, ndimsi
            dvsigp(i) = sigel(i)-1.5d0*deuxmu*dp*sig0(i)
            dvsigp(i) = dvsigp(i)/(1.5d0*deuxmu*trans+1.d0)
            sigp(i) = dvsigp(i)+trsigp*kron(i)
        end do
!
! 4.2.3 - CALCUL DE VIP ET RMOY
!
        do k = 1, nb_phase
            if (phase(k) .gt. 0.d0) then
                vip(k) = vi(k)+dp
            else
                vip(k) = 0.d0
            end if
        end do
        vip(IDX_I_EPSEQ) = 0.d0
        if (phase(nb_phase) .gt. 0.d0) then
            if (l_hard_isotline) then
                vip(IDX_I_EPSEQ) = vip(IDX_I_EPSEQ)+(1-fmel)*h(nb_phase)*vip(nb_phase)
            end if
            if (l_hard_isotnlin) then
                vip(IDX_I_EPSEQ) = vip(IDX_I_EPSEQ)+(1-fmel)*(r(nb_phase)-sy(nb_phase))
            end if
        end if
!
        if (zalpha .gt. 0.d0) then
            do k = 1, nb_phase-1
                if (l_hard_isotline) then
                    vip(IDX_I_EPSEQ) = vip(IDX_I_EPSEQ)+fmel*phase(k)*h(k)*vip(k)/zalpha
                end if
                if (l_hard_isotnlin) then
                    vip(IDX_I_EPSEQ) = vip(IDX_I_EPSEQ)+fmel*phase(k)*(r(k)-sy(k))/zalpha
                end if
            end do
        end if
    end if
!
! *******************************
! 5 - MATRICE TANGENTE DSIGDF
! *******************************
!
    if (rigi) then
        mode = 2
        if (l_visc) mode = 1
        dsidep(1:6, 1:6) = 0.d0
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
                                if (phase(k) .gt. 0.d0) then
                                    dv = dv+fmel*(phase(k)/zalpha)* &
                                         (eta(k)/n(k)/dt)*((dp/dt)**n0(k))
                                end if
                            end do
                        end if
                    end if
                    coef2 = hmoy+dv
                    coef2 = (1.5d0*deuxmu*trans+1.d0)*coef2
                    coef2 = (1.5d0*deuxmu)+coef2
                    coef2 = 1/coef2-dp/sieleq
                    coef2 = ((1.5d0*deuxmu)**2)*coef2
                end if
            end if
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                if (mode .eq. 2) coef2 = ((1.5d0*deuxmu)**2)/(1.5d0*deuxmu+hmoy)
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
