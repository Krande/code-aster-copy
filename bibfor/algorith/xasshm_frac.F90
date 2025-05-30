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
subroutine xasshm_frac(ds_thm, &
                       nddls, nddlm, nnop, nnops, &
                       lact, elrefp, elrefc, elc, contac, &
                       dimuel, nptf, &
                       jptint, igeom, jbasec, &
                       jcohes, jcoheo, &
                       nlact, cface, rinstp, &
                       rinstm, carcri, fpg, ncompv, &
                       compor, jmate, ndim, idepm, idepd, &
                       pla, algocr, rela, ifa, ipgf, matri, &
                       cohes, coheo, jheavn, ncompn, ifiss, &
                       nfiss, nfh, jheafa, ncomph, pos)
!
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/xfract.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xhlan5.h"
#include "asterfort/xhmco3.h"
#include "asterfort/xhmco4.h"
#include "asterfort/xhmsa6.h"
#include "asterfort/xmathm.h"
#include "asterfort/xmmata.h"
#include "asterfort/xmmatb.h"
#include "asterfort/xmmatc.h"
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
    integer :: nddls, nnop, dimuel, i, ndim, nnops, jheavn
    integer :: nddlm, contac, jmate, ncompv, nvec, pla(27), ncompn
    integer :: nptf, nfiss, jcohes, jcoheo
    integer :: ipgf, ifa, cface(30, 6), algocr, idepd, idepm, pos(16)
    integer :: jptint, igeom, jbasec, nlact(2), lact(16), ibid, ino
    integer :: ifiss, nfh, jheafa, ncomph
    real(kind=8) :: dt, parm_theta, ta1, cohes(5), rela, g(3), vihydr(64)
    real(kind=8) :: rinstp, rinstm, carcri(*), jac, raug
    real(kind=8) :: ffp(27), ffpc(27), ffc(16), dfdic(nnops, 3)
    real(kind=8) :: dfbid(27, 3), lamb(3), wsaut(3), wsautm(3)
    real(kind=8) :: nd(3), tau1(3), tau2(3)
    real(kind=8) :: matri(560, 560), dffc(16, 3), rho11, w11
    real(kind=8) :: saut(3), gradpf(3), q1, q2, dpf, rho110
    real(kind=8) :: q1m, q2m, gradpfm(3), sautm(3)
    real(kind=8) :: w11m, rho11m, alpha(5), coheo(5)
    real(kind=8) :: pf, psup, pinf, ffp2(27), t
    real(kind=8) :: dsidep(6, 6), delta(6), p(3, 3), sigma(6)
    real(kind=8) :: am(3), ad(3), r
    real(kind=8) :: temp
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

!   INITIALISATION DE LA DIMENSION DE LA MATRICE DE TRAVAIL
!
    dfdic(:, :) = 0.d0
    dffc(:, :) = 0.d0
    ffp(:) = 0.d0
    ffpc(:) = 0.d0
    ffc(:) = 0.d0
    ffp2(:) = 0.d0
    vihydr(:) = 0.d0
!
!         CALCUL DU PRODUIT DU JACOBIEN AVEC LE jac D'INTEGRATION, DES FONCTIONS
!         DE FORME POUR L'ELEMENT PARENT QUADRATIQUE ET DE LA NORMALE A LA
!         FACETTE
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
!
!   CALCUL DES FONCTIONS DE FORME DE CONTACT
    call xmofhm(lact, nlact, nnops, ffpc, ffc)
!   CALCUL DU GRADIENT DES FONCTIONS DE FORME DE CONTACT
!
    call xmodfc(lact, nlact, nnops, dfdic, dffc, ndim)
!
    if (algocr .eq. 3) then
        if ((nint(rela) .eq. 3) .or. (nint(rela) .eq. 4)) then
!
            nvec = 2
            job = 'MATRICE'
            call xfract(ds_thm, nvec, nnop, nnops, nddls, nddlm, &
                        ndim, pla, zr(idepd), zr(idepm), &
                        ffp, ffc, dffc, saut, gradpf, &
                        q1, q2, dpf, q1m, q2m, sautm, &
                        gradpfm, pf, ffp2, psup, pinf, &
                        job, zi(jmate), &
                        t, dimuel, lamb, jheavn, ncompn, ifiss, &
                        nfiss, nfh, ifa, jheafa, ncomph, &
                        contac)
!
!          CALCUL DU CHANGEMENT DE BASE POUR LE SAUT DE DEPLACEMENT
!
            call xsautl(ndim, nd, tau1, tau2, saut, sautm, p, am, ad)
!
            do i = 1, ndim
                am(i) = -am(i)
            end do
!
!          CALCUL DE LA VARIABLE INTERNE (MASSE VOLUMIQUE DU LIQUIDE
!          CIRCULANT DANS LA FRACTURE)
!
            job = 'MATRICE'
            call xvinhm(ds_thm, zi(jmate), ndim, &
                        cohes, dpf, saut, sautm, nd, lamb, &
                        w11m, rho11m, alpha, job, pf, &
                        rho11, w11, ipgf, rela, dsidep, &
                        delta, r, am)
!
!          CALCUL DES MATRICES (CF. DOC R7.02.18)
            call xmathm(ds_thm, ndim, &
                        nnops, nnop, nddls, nddlm, ffc, &
                        pla, nd, jac, ffp, ffp2, dt, parm_theta, saut, &
                        dffc, rho11, gradpf, matri, &
                        dsidep, p, r, jheavn, ncompn, ifiss, &
                        nfiss, nfh, ifa, jheafa, ncomph)
!
!          ACTUALISATION INDICATEUR PREDICTION / CORRECTION
!
            coheo(1) = cohes(1)
            coheo(2) = cohes(2)
            coheo(3) = 2.d0
            coheo(4) = cohes(4)
            coheo(5) = cohes(5)
!
        else if (nint(rela) .eq. 5) then
!
! --- CALCUL DES MATRICES "MORTAR" QUI NE DEPENDENT
! --- PAS DE LA LOI DE COMPORTEMENT
!
            call xhmco4(ndim, nnop, nnops, pla, nd, tau1, &
                        tau2, ffc, nddls, jac, ffp, &
                        nddlm, matri, ifiss, nfiss, nfh, &
                        ifa, jheafa, ncomph, jheavn, ncompn)
!
            nvec = 2
            job = 'MATRICE'
            call xfract(ds_thm, nvec, nnop, nnops, nddls, nddlm, &
                        ndim, pla, zr(idepd), zr(idepm), &
                        ffp, ffc, dffc, saut, gradpf, &
                        q1, q2, dpf, q1m, q2m, sautm, &
                        gradpfm, pf, ffp2, psup, pinf, &
                        job, zi(jmate), &
                        t, dimuel, lamb, jheavn, ncompn, ifiss, &
                        nfiss, nfh, ifa, jheafa, ncomph, contac)
!
            do ino = 1, nnops
                champ = 'LAMBDA'
                call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                            pla, lamb, nvec, champ, job, dpf)
!
                champ = 'W'
                call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                            pla, wsaut, nvec, champ, job, dpf)
!
                champ = 'WM'
                call xhlan5(ino, idepd, idepm, ibid, lact, ndim, &
                            pla, wsautm, nvec, champ, job, dpf)
!
                do i = 1, ncompv
                    cohes(i) = zr(jcohes+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+i)
                end do
!
                call xhmsa6(ds_thm, ndim, ipgf, zi(jmate), lamb, wsaut, nd, &
                            tau1, tau2, cohes, job, rela, &
                            alpha, dsidep, sigma, p, am, raug, &
                            wsautm, dpf, rho110)
                call xhmco3(ino, ndim, dsidep, pla, p, &
                            ffc, jac, raug, matri)
!
! --- ACTUALISATION INDICATEUR PREDICTION / CORRECTION
!
                coheo(1) = cohes(1)
                coheo(2) = cohes(2)
                coheo(3) = 2.d0
                coheo(4) = cohes(4)
                coheo(5) = cohes(5)
!
                do i = 1, ncompv
                    zr(jcoheo+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+i) = coheo(i)
                end do
                vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+1) = alpha(4)
                vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+2) = alpha(5)
            end do
!
! --- INTERPOLATION DE RHO11 AU POINT DE GAUSS
!
            rho11 = rho110
            do ino = 1, nnops
                rho11 = rho11+ffc(ino)*vihydr(2*nnops*(pos(ino)-1)+2*(ino-1)+1)
            end do
!
! --- CALCUL DES MATRICES HYDROS
!
            call xmmatc(ndim, nnops, nddls, nddlm, ffc, &
                        pla, jac, ffp2, matri, &
                        jheavn, ncompn, ifiss, nfiss, &
                        nfh, ifa, jheafa, ncomph)
!
            call xmmatb(ndim, nnops, nddls, nddlm, ffc, &
                        pla, dt, parm_theta, jac, ffp2, matri, &
                        jheavn, ncompn, ifiss, nfiss, nfh, &
                        ifa, jheafa, ncomph)
!
            call xmmata(ds_thm, ndim, nnops, nnop, nddls, nddlm, saut, &
                        nd, pla, ffc, dffc, matri, rho11, &
                        gradpf, ffp, dt, parm_theta, jac, &
                        jheavn, ncompn, ifiss, nfiss, &
                        nfh, ifa, jheafa, ncomph)
!
        end if
    end if
!
end subroutine
