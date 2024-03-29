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
! aslint: disable=W1504,W1306
! person_in_charge: daniele.colombo at ifpen.fr
!
subroutine xasshv_frac(ds_thm, &
                       nddls, nddlm, nnop, nnops, &
                       lact, elrefp, elrefc, elc, contac, &
                       dimuel, nface, npgf, nbspg, nptf, &
                       jcohes, jptint, igeom, jbasec, &
                       nlact, cface, rinstp, &
                       rinstm, carcri, fpg, ncompv, vect, &
                       compor, jmate, ndim, idepm, idepd, pla, &
                       algocr, rela, jheavn, ncompn, ifiss, nfiss, &
                       nfh, jheafa, ncomph, pos)
!
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/xfract.h"
#include "asterfort/xhlag4.h"
#include "asterfort/xhlan5.h"
#include "asterfort/xhmsa6.h"
#include "asterfort/xhvco4.h"
#include "asterfort/xhvco5.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xvecha.h"
#include "asterfort/xvechb.h"
#include "asterfort/xvechc.h"
#include "asterfort/xvechm.h"
#include "asterfort/xmodfc.h"
#include "asterfort/xmofhm.h"
#include "asterfort/xsautl.h"
#include "asterfort/xvinhm.h"
#include "asterfort/thmGetParaBiot.h"
#include "asterfort/thmGetBehaviourVari.h"
#include "asterfort/thmGetBehaviour.h"
#include "asterfort/thmGetParaCoupling.h"
#include "asterfort/thmGetBehaviourChck.h"
#include "asterfort/Behaviour_type.h"
!
! ======================================================================
!
! ROUTINE MODELE HM-XFEM (CAS DE LA FRACTURE)
!
! CALCUL DES MATRICES POUR LE FRACTURE (OPTION RIGI_CONT)
!
! ----------------------------------------------------------------------
! IO  ds_thm           : datastructure for THM

    type(THM_DS), intent(inout) :: ds_thm
    integer :: nddls, nnop, dimuel, i, ndim, nnops
    integer :: nddlm, contac, jmate, ncompv, nvec, pla(27), pos(16)
    integer :: nface, npgf, nbspg, isspg, nptf, cface(30, 6)
    integer :: ipgf, ifa, jcohes, lact(16), algocr, jheavn, ncompn
    integer :: jptint, igeom, jbasec, nlact(2), nfiss, ibid, ino
    integer :: ifiss, nfh, jheafa, ncomph
    integer :: idepm, idepd
    real(kind=8) :: dt, parm_theta, ta1, cohes(5), lamb(3), g(3)
    real(kind=8) :: rinstp, rinstm, carcri(*), jac, r
    real(kind=8) :: ffp(27), ffpc(27), ffc(16), dfdic(nnops, 3)
    real(kind=8) :: dfbid(27, 3), am(3), ad(3)
    real(kind=8) :: nd(3), tau1(3), tau2(3), raug, sigma(6), rho110
    real(kind=8) :: vect(560), dffc(16, 3)
    real(kind=8) :: saut(3), gradpf(3), q1, q2, dpf
    real(kind=8) :: q1m, q2m, gradpfm(3), sautm(3)
    real(kind=8) :: w11m, rho11m, alpha(5), w11
    real(kind=8) :: pf, psup, pinf, ffp2(27), t, vihydr(64), temp
    real(kind=8) :: rho11, wsaut(3), mu(3), wsautm(3)
    real(kind=8) :: dsidep(6, 6), delta(6), rela, p(3, 3)
    character(len=8) :: elrefp, elrefc, elc, fpg, job, champ
    character(len=16):: compor(*)

!   DETERMINATION DES CONSTANTES TEMPORELLES (INSTANT+THETA SCHEMA)
    dt = rinstp-rinstm
    parm_theta = carcri(PARM_THETA_THM)
    ta1 = 1.d0-parm_theta
!
! - Get parameters for behaviour
!
    call thmGetBehaviour(compor, ds_thm)
!
! - Get parameters for internal variables
!
    call thmGetBehaviourVari(ds_thm)
!
! - Some checks between behaviour and model
!
    call thmGetBehaviourChck(ds_thm)
!
! - Get parameters for coupling
!
    temp = 0.d0
    call thmGetParaCoupling(ds_thm, zi(jmate), temp)
!
! - Get Biot parameters (for porosity evolution)
!
    call thmGetParaBiot(zi(jmate), ds_thm)
!
    dfdic(:, :) = 0.d0
    dffc(:, :) = 0.d0
    ffp(:) = 0.d0
    ffpc(:) = 0.d0
    ffc(:) = 0.d0
    ffp2(:) = 0.d0
    vihydr(:) = 0.d0
!
!   BOUCLE SUR LES FACETTES DE CONTACT
    do ifa = 1, nface

!      BOUCLE SUR LES POINTS D'INTEGRATION DE LA FACETTE DE CONTACT COURANTE IFA
        do ipgf = 1, npgf

!         DECALAGE POUR ACCEDER AU POINT DE GAUSS IPGF DE LA FACETTE DE CONTACT
!         COURANTE IFA

            isspg = npgf*(ifa-1)+ipgf

!          POUR L'ACTUALISATION DES VARIABLES INTERNES
            if (contac .ne. 2) then
                do i = 1, ncompv
                    cohes(i) = zr(jcohes+ncompv*(nbspg+isspg-1)-1+i)
                end do
            end if
!
!         CALCUL DU PRODUIT DU JACOBIEN AVEC LE jac D'INTEGRATION, DES FONCTIONS
!         DE FORME POUR L'ELEMENT PARENT QUADRATIQUE ET DE LA NORMALE A LA FACETTE
            if (ndim .eq. 2) then
                call xjacf2(elrefp, elrefc, elc, ndim, fpg, &
                            jptint, ifa, cface, nptf, ipgf, &
                            nnop, nnops, igeom, jbasec, g, jac, &
                            ffp, ffpc, dfbid, nd, tau1, dfdic)
            elseif (ndim .eq. 3) then
                call xjacff(elrefp, elrefc, elc, ndim, fpg, &
                            jptint, ifa, cface, ipgf, nnop, &
                            nnops, igeom, jbasec, g, jac, ffp, &
                            ffpc, dfbid, nd, tau1, tau2, dfdic)
            end if
!
            ffp2(1:nnops) = ffpc(1:nnops)
!         CALCUL DES FONCTIONS DE FORME DE CONTACT
!
            call xmofhm(lact, nlact, nnops, ffpc, ffc)
!
!         CALCUL DU GRADIENT DES FONCTIONS DE FORME DE CONTACT
!
            call xmodfc(lact, nlact, nnops, dfdic, dffc, ndim)
!
            if (algocr .eq. 3) then
                if ((nint(rela) .eq. 3) .or. (nint(rela) .eq. 4)) then
!
                    nvec = 2
                    job = 'VECTEUR'
                    call xfract(ds_thm, nvec, nnop, nnops, nddls, nddlm, &
                                ndim, pla, zr(idepd), zr(idepm), &
                                ffp, ffc, dffc, saut, gradpf, &
                                q1, q2, dpf, q1m, q2m, sautm, &
                                gradpfm, pf, ffp2, psup, pinf, &
                                job, zi(jmate), &
                                t, dimuel, lamb, jheavn, ncompn, &
                                ifiss, nfiss, nfh, ifa, jheafa, &
                                ncomph, contac)
!
!                CALCUL DU CHANGEMENT DE BASE POUR LE SAUT DE DEPLACEMENT
!
                    call xsautl(ndim, nd, tau1, tau2, saut, sautm, p, am, ad)
!
                    do i = 1, ndim
                        am(i) = -am(i)
                    end do
!
!                CALCUL DE LA VARIABLE INTERNE (MASSE VOLUMIQUE DU LIQUIDE
!                CIRCULANT DANS LA FRACTURE)
!
                    job = 'VECTEUR'
                    call xvinhm(ds_thm, zi(jmate), ndim, &
                                cohes, dpf, saut, sautm, nd, lamb, &
                                w11m, rho11m, alpha, job, pf, &
                                rho11, w11, ipgf, rela, dsidep, &
                                delta, r, am)
!
!                CALCUL DES SECONDS MEMBRES (CF. DOC R7.02.18)
!
                    call xvechm(ds_thm, nnops, nddls, nddlm, ndim, pla, &
                                saut, sautm, nd, ffc, w11, w11m, jac, &
                                q1, dt, parm_theta, q1m, ta1, q2, q2m, dffc, &
                                rho11, gradpf, rho11m, gradpfm, ffp2, &
                                vect, &
                                ffp, nnop, delta, lamb, am, r, p, &
                                psup, pinf, pf, ncompn, jheavn, ifiss, &
                                nfiss, nfh, ifa, jheafa, ncomph)
                else if (nint(rela) .eq. 5) then
                    nvec = 2
                    job = 'VECTEUR'
                    call xfract(ds_thm, nvec, nnop, nnops, nddls, nddlm, &
                                ndim, pla, zr(idepd), zr(idepm), &
                                ffp, ffc, dffc, saut, gradpf, &
                                q1, q2, dpf, q1m, q2m, sautm, &
                                gradpfm, pf, ffp2, psup, pinf, &
                                job, zi(jmate), &
                                t, dimuel, lamb, jheavn, ncompn, &
                                ifiss, nfiss, nfh, ifa, jheafa, &
                                ncomph, contac)
!           CALCUL W AU POINT DE GAUSS
                    champ = 'W'
                    call xhlag4(ffc, idepd, idepm, lact, ndim, &
                                nnops, pla, wsaut, nvec, champ)
!
!           CALCUL MU AU POINT DE GAUSS
                    champ = 'MU'
                    call xhlag4(ffc, idepd, idepm, lact, ndim, &
                                nnops, pla, mu, nvec, champ)
! --- CALCUL DES SECONDS MEMBRES DE COHESION
!
                    call xhvco5(ndim, nnop, nnops, pla, nd, &
                                tau1, tau2, mu, nddls, jac, &
                                ffc, ffp, nddlm, wsaut, &
                                saut, vect, ifiss, nfiss, nfh, &
                                ifa, jheafa, ncomph, &
                                jheavn, ncompn, pf)
!
                    do ino = 1, nnops
                        do i = 1, ncompv
                            cohes(i) = zr(jcohes+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+i)
                        end do
                        champ = 'LAMBDA'
                        call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                                    pla, lamb, nvec, champ, job, dpf)
                        champ = 'W'
                        call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                                    pla, wsaut, nvec, champ, job, dpf)
                        champ = 'WM'
                        call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                                    pla, wsautm, nvec, champ, job, dpf)
!
                        call xhmsa6(ds_thm, ndim, ipgf, zi(jmate), lamb, wsaut, nd, &
                                    tau1, tau2, cohes, job, rela, &
                                    alpha, dsidep, sigma, p, am, raug, &
                                    wsautm, dpf, rho110)
                        call xhvco4(ino, ndim, sigma, lamb, pla, &
                                    jac, ffc, p, raug, vect)
                        vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+1) = alpha(4)
                        vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+2) = alpha(5)
                    end do
!
! --- INTERPOLATION DES VI HYDROS AU POINT DE GAUSS
!
                    rho11 = rho110
                    rho11m = rho110
                    w11 = 0.d0
                    w11m = 0.d0
                    do ino = 1, nnops
                        rho11 = rho11+ffc(ino)*vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+1)
                    rho11m = rho11m+ffc(ino)*zr(jcohes+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+4)
                        w11 = w11+ffc(ino)*vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+2)
                        w11m = w11m+ffc(ino)*zr(jcohes+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+5)
                    end do
!
! --- CALCUL DES MATRICES HYDROS
!
!
                    call xvecha(ds_thm, ndim, pla, nnops, saut, &
                                sautm, nd, ffc, w11, w11m, jac, &
                                q1, q1m, q2, q2m, dt, parm_theta, ta1, &
                                dffc, rho11, gradpf, rho11m, &
                                gradpfm, vect)
!
                    call xvechb(nnops, nddls, nddlm, ndim, &
                                ffp2, q1, dt, parm_theta, jac, q1m, ta1, &
                                q2, q2m, vect, ncompn, jheavn, ifiss, &
                                nfiss, nfh, ifa, jheafa, ncomph)
!
                    call xvechc(nnops, pla, ffc, pinf, &
                                pf, psup, jac, vect)
!
                end if
            end if
        end do
    end do
end subroutine
