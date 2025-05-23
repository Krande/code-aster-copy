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
! person_in_charge: daniele.colombo at ifpen.fr
! aslint: disable=W1504,W1306
!
subroutine xassha_frac(ds_thm, &
                       nddls, nddlm, nnop, nnops, &
                       lact, elrefp, elrefc, elc, contac, &
                       dimuel, nface, npgf, nbspg, nptf, &
                       jcohes, jptint, igeom, jbasec, &
                       nlact, cface, fpg, ncompv, &
                       compor, jmate, ndim, idepm, idepd, jcoheo, incoca, &
                       pla, rela, algocr, jheavn, ncompn, ifiss, &
                       nfiss, nfh, jheafa, ncomph, pos)
!
    use THM_type
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/xfract.h"
#include "asterfort/xhlan5.h"
#include "asterfort/xhmsa6.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xmodfc.h"
#include "asterfort/xmofhm.h"
#include "asterfort/xsautl.h"
#include "asterfort/xvinhm.h"
#include "asterfort/thmGetBehaviourVari.h"
#include "asterfort/thmGetBehaviour.h"
#include "asterfort/thmGetParaCoupling.h"
#include "asterfort/thmGetBehaviourChck.h"
!
! ======================================================================
!
! ROUTINE MODELE HM-XFEM (CAS DE LA FRACTURE)
!
! ACTUALISATION DES VARIABLES INTERNES (OPTION XCVBCA)
!
! ----------------------------------------------------------------------
! IO  ds_thm           : datastructure for THM
    type(THM_DS), intent(inout) :: ds_thm
    integer :: nddls, nnop, dimuel, i, ndim, nnops, lact(16)
    integer :: nddlm, contac, jmate, ncompv, nvec, pla(27), incoca
    integer :: nface, npgf, nbspg, isspg, nptf
    integer :: ipgf, ifa, jcohes, algocr, jheavn, ncompn, nfiss
    integer :: jptint, igeom, jbasec, nlact(2), cface(30, 6), jcoheo
    integer :: ifiss, nfh, jheafa, ncomph, pos(16), ino
    integer :: idepm, idepd, ibid
    real(kind=8) :: cohes(5), jac, lamb(3), am(3), ad(3), g(3)
    real(kind=8) :: ffp(27), ffpc(27), ffc(16), dfdic(nnops, 3)
    real(kind=8) :: dfbid(27, 3), wsaut(3), wsautm(3), sigma(6)
    real(kind=8) :: nd(3), tau1(3), tau2(3)
    real(kind=8) :: dffc(16, 3), saut(3), gradpf(3)
    real(kind=8) :: q1, q2, dpf, q1m, q2m, w11m, rho11m
    real(kind=8) :: gradpfm(3), sautm(3), alpha(5)
    real(kind=8) :: pf, psup, pinf, ffp2(27), t, eps, rela, temp
    real(kind=8) :: rho11, w11, rho110, raug
    real(kind=8) :: dsidep(6, 6), delta(6), p(3, 3), r
    character(len=8) :: elrefp, elrefc, elc, fpg, job, champ
    character(len=16):: compor(*)

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
    dfdic(:, :) = 0.d0
    dffc(:, :) = 0.d0
    ffp(:) = 0.d0
    ffpc(:) = 0.d0
    ffc(:) = 0.d0
    ffp2(:) = 0.d0
    nd(:) = 0.d0
    tau1(:) = 0.d0
    tau2(:) = 0.d0
!   BOUCLE SUR LES FACETTES DE CONTACT
    do ifa = 1, nface
!      BOUCLE SUR LES POINTS D'INTEGRATION DE LA FACETTE DE CONTACT COURANTE IFA
        do ipgf = 1, npgf
!         DECALAGE POUR ACCEDER AU POINT DE GAUSS IPGF DE LA FACETTE DE CONTACT
!         COURANTE IFA
!
            isspg = npgf*(ifa-1)+ipgf
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
!
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
!          POUR L'ACTUALISATION DES VARIABLES INTERNES
!
                    do i = 1, ncompv
                        cohes(i) = zr(jcohes+ncompv*(nbspg+isspg-1)-1+i)
                    end do
                    nvec = 1
                    job = 'ACTU_VI'
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
                    job = 'ACTU_VI'
                    call xvinhm(ds_thm, zi(jmate), ndim, &
                                cohes, dpf, saut, sautm, nd, lamb, &
                                w11m, rho11m, alpha, job, pf, &
                                rho11, w11, ipgf, rela, dsidep, &
                                delta, r, am)
!
                    do i = 1, ncompv
                        zr(jcoheo+ncompv*(nbspg+isspg-1)-1+i) = alpha(i)
                    end do
                    eps = r8prem()
                    ASSERT((alpha(1)+eps) .ge. cohes(1))
                else if (nint(rela) .eq. 5) then
                    job = 'ACTU_VI'
                    nvec = 1
                    do ino = 1, nnops
                        do i = 1, ncompv
                            cohes(i) = zr(jcohes+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+i)
                        end do
!
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
                        call xhmsa6(ds_thm, ndim, ipgf, zi(jmate), lamb, wsaut, nd, &
                                    tau1, tau2, cohes, job, rela, &
                                    alpha, dsidep, sigma, p, am, raug, &
                                    wsautm, dpf, rho110)
!
                        do i = 1, ncompv
                            zr(jcoheo+ncompv*nnops*(pos(ino)-1)+ncompv*(ino-1)-1+i) = alpha(i)
                        end do
                        eps = r8prem()
                        ASSERT((alpha(1)+eps) .ge. cohes(1))
                    end do
                end if
            end if
!         CHAMPS LIES AU CONTACT INUTILES POUR LE COHESIF
            incoca = 0
        end do
    end do

end subroutine
