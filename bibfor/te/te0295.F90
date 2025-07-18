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
subroutine te0295(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/cgverho.h"
#include "asterfort/chauxi.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/gbil3d.h"
#include "asterfort/jevech.h"
#include "asterfort/nmgeom.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!
!     BUT:
!       CALCUL DES FACTEURS D'INTENSITÉ DES CONTRAINTES
!       A PARTIR DE LA FORME BILINEAIRE SYMETRIQUE G ET
!       DES DEPLACEMENTS SINGULIERS EN FOND DE FISSURE
!      POUR LES ELEMENTS ISOPARAMETRIQUES 3D
!
!       OPTION : 'CALC_K_G'    (CHARGES REELLES)
!                'CALC_K_G_F'  (CHARGES FONCTIONS)
!
!
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: icodre(4)
    integer(kind=8) :: ipoids, ivf, idfde, nno, kp, npg, compt, ier, nnos
    integer(kind=8) :: jgano, ibalo, icour, isigi
    integer(kind=8) :: igeom, ithet, ificg, irota, ipesa, idepl, iret, ncmp
    integer(kind=8) :: imate, iforc, iforf, itemps, k, i, j, kk, l, ndim, ino, ipuls
    integer(kind=8) :: jlsn, jlst, jtab(7)
!
!
    real(kind=8) :: r8bid, rac2
    real(kind=8) :: dfdi(60), f(3, 3), eps(6), fno(81)
    real(kind=8) :: dudm(3, 4), dfdm(3, 4), dtdm(3, 4), der(4)
    real(kind=8) :: u1l(3), u2l(3), u3l(3), dfvdm(3, 4), epsref(6)
    real(kind=8) :: du1dm(3, 4), du2dm(3, 4), du3dm(3, 4)
    real(kind=8) :: p(3, 3), invp(3, 3), sigin(6), dsigin(6, 3)
    real(kind=8) :: courb(3, 3, 3)
    real(kind=8) :: rhocst, rho, om, omo, rbid, e, nu, alpha, tref
    real(kind=8) :: thet, tpg(27), tno(20), tgdm(3), ttrg, la, mu, ka
    real(kind=8) :: xg, yg, zg
    real(kind=8) :: c1, c2, c3, rg, phig
    real(kind=8) :: val(1), valres(4)
    real(kind=8) :: coeff, coeff3
    real(kind=8) :: guv, guv1, guv2, guv3, k1, k2, k3, g, poids
    real(kind=8) :: k3a, ttrgv, tgvdm(3)
    real(kind=8) :: valpar(4), lsng, lstg, puls, coef, ffp(27)
!
    character(len=4) :: fami
    character(len=8) :: nompar(4)
    character(len=16) :: nomres(4)
    character(len=32) :: phenom
    character(len=16), pointer :: compor(:) => null()
!
    aster_logical :: lcour, fonc, lpesa, lrota, l_not_zero
!
! ----------------------------------------------------------------------
!
!
!
    rac2 = sqrt(2.d0)
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! --- RECUPERATION DES CHAMPS IN
!
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PBASLOR', 'L', ibalo)
    call jevech('PCOURB', 'L', icour)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
! NOMBRE DE COMPOSANTES DES TENSEURS
    ncmp = 2*ndim
!
! --- RECUPERATION DU CHAMP OUT
!
    call jevech('PGTHETA', 'E', ificg)
!
! NOMBRE DE COMPOSANTES DES TENSEURS
    ncmp = 2*ndim
    g = 0.d0
    k1 = 0.d0
    k2 = 0.d0
    k3 = 0.d0
    coeff = 1.d0
    coeff3 = 1.d0
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    nomres(4) = 'RHO'
!
! --- PAS DE CALCUL DE G POUR LES ELEMENTS OU THETA EST NULLE
!
    compt = 0
    do i = 1, nno
        thet = 0.d0
        do j = 1, ndim
            thet = thet+abs(zr(ithet+ndim*(i-1)+j-1))
        end do
        if (thet .lt. r8prem()) compt = compt+1
    end do
    if (compt .eq. nno) goto 999
!
! --- VERIFS DE COHERENCE RHO <-> PESANTEUR, ROTATION, PULSATION
!
    if (.not. cgverho(imate)) call utmess('F', 'RUPTURE1_26')
!
! --- RECUPERATION DES FORCES
!
    if (option .eq. 'CALC_K_G_XFEM_F') then
        fonc = .true.
        call jevech('PFFVOLU', 'L', iforf)
        call jevech('PINSTR', 'L', itemps)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        if (ndim .eq. 2) then
            nompar(3) = 'INST'
        else if (ndim .eq. 3) then
            nompar(3) = 'Z'
            nompar(4) = 'INST'
        end if
        valpar(ndim+1) = zr(itemps)
    else if (option .eq. 'CALC_K_G_XFEM') then
        fonc = .false.
        call jevech('PFRVOLU', 'L', iforc)
    else
        ASSERT(.false.)
    end if
!
    lpesa = .false.
    call tecach('ONO', 'PPESANR', 'L', iret, nval=7, &
                itab=jtab)
    ipesa = jtab(1)
    if (iret .eq. 0) then
        lpesa = .true.
    end if
!
    lrota = .false.
    call tecach('ONO', 'PROTATR', 'L', iret, nval=7, &
                itab=jtab)
    irota = jtab(1)
    if (iret .eq. 0) then
        lrota = .true.
    end if
!
! --- VERFICATION DU COMPORTEMENT
!
    if ((compor(1) .ne. 'ELAS') .or. (compor(4) .eq. 'COMP_INCR')) then
        if (compor(1) .ne. 'ELAS') then
            call utmess('F', 'RUPTURE1_24')
        end if
    end if
!
! --- RECUPERATION DE LA PULSATION
!
    call tecach('ONO', 'PPULPRO', 'L', iret, nval=7, &
                itab=jtab)
    ipuls = jtab(1)
    if (iret .eq. 0) then
        puls = zr(ipuls)
    else
        puls = 0.d0
    end if
!
! --- RECUPERATION DES CHARGES
!
    if (fonc) then
        do i = 1, nno
            do j = 1, ndim
                valpar(j) = zr(igeom+ndim*(i-1)+j-1)
            end do
            do j = 1, ndim
                kk = ndim*(i-1)+j
                call fointe('FM', zk8(iforf+j-1), ndim+1, nompar, valpar, &
                            fno(kk), ier)
            end do
        end do
    else
        do i = 1, nno
            do j = 1, ndim
                fno(ndim*(i-1)+j) = zr(iforc+ndim*(i-1)+j-1)
            end do
        end do
    end if
!
! --- RECUPERATION DE LA PESANTEUR ET DE LA ROTATION
!
    if (lpesa .or. lrota) then
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
        call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                    ' ', phenom, 1, ' ', [rbid], &
                    1, 'RHO', val, icodre, 1)
        rhocst = val(1)
        if (lpesa) then
            do i = 1, nno
                do j = 1, ndim
                    kk = ndim*(i-1)+j
                    fno(kk) = fno(kk)+rhocst*zr(ipesa)*zr(ipesa+j)
                end do
            end do
        end if
!
        if (lrota) then
            om = zr(irota)
            do i = 1, nno
                omo = 0.d0
                do j = 1, ndim
                    omo = omo+zr(irota+j)*zr(igeom+ndim*(i-1)+j-1)
                end do
                do j = 1, ndim
                    kk = ndim*(i-1)+j
                    fno(kk) = fno(kk)+rhocst*om*om*(zr(igeom+kk-1)-omo*zr( &
                                                    irota+j))
                end do
            end do
        end if
    end if
!
! --- RECUPERATION DE LA TEMPERATURE
!
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tref, iret)
    if (iret .ne. 0) tref = 0.d0
    do kp = 1, npg
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                    1, tpg(kp), iret)
        if (iret .ne. 0) tpg(kp) = 0.d0
    end do
!
    do ino = 1, nno
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, &
                    1, tno(ino), iret)
        if (iret .ne. 0) tno(ino) = 0.d0
    end do
! --- RECUPERATION DE LA CONTRAINTE INITIALE
    call tecach('ONO', 'PSIGINR', 'L', iret, iad=isigi)
!
! ----------------------------------------------------------------------
!
!     BOUCLE SUR LES POINTS DE GAUSS
!
! ----------------------------------------------------------------------
!
    do kp = 1, npg
!INITIALISATIONS
        l = (kp-1)*nno
        xg = 0.d0
        yg = 0.d0
        zg = 0.d0
        lsng = 0.d0
        lstg = 0.d0
        do i = 1, 3
            tgdm(i) = 0.d0
            tgvdm(i) = 0.d0
            do j = 1, 4
                dudm(i, j) = 0.d0
                du1dm(i, j) = 0.d0
                du2dm(i, j) = 0.d0
                du3dm(i, j) = 0.d0
                dtdm(i, j) = 0.d0
                dfdm(i, j) = 0.d0
                dfvdm(i, j) = 0.d0
            end do
        end do
        do i = 1, 6
            sigin(i) = 0.d0
            epsref(i) = 0.d0
            do j = 1, 3
                dsigin(i, j) = 0.d0
            end do
        end do
!
! ----- CALCUL DES ELEMENTS CINEMATIQUES (MATRICES F ET E)
!       EN UN PT DE GAUSS
!
        call nmgeom(ndim, nno, .false._1, .false._1, zr(igeom), &
                    kp, ipoids, ivf, idfde, zr(idepl), &
                    .true._1, poids, dfdi, f, eps, &
                    rbid)
!
! ----- CALCULS DES GRADIENTS DE U (DUDM),THETA (DTDM) ET FORCE(DFDM)
! ----- DU GRADIENT DE TEMPERATURE AUX POINTS DE GAUSS (TGDM)
! ----- ET LEVEL SETS
!
        do i = 1, nno
            der(1) = dfdi(i)
            der(2) = dfdi(i+nno)
            der(3) = dfdi(i+2*nno)
            der(4) = zr(ivf+l+i-1)
!
            xg = xg+zr(igeom-1+ndim*(i-1)+1)*der(4)
            yg = yg+zr(igeom-1+ndim*(i-1)+2)*der(4)
            zg = zg+zr(igeom-1+ndim*(i-1)+3)*der(4)
!
            lsng = lsng+zr(jlsn-1+i)*der(4)
            lstg = lstg+zr(jlst-1+i)*der(4)
!
            do j = 1, ndim
                tgdm(j) = tgdm(j)+tno(i)*der(j)
                do k = 1, ndim
                    dudm(j, k) = dudm(j, k)+zr(idepl+ndim*(i-1)+j-1)*der(k)
                    dtdm(j, k) = dtdm(j, k)+zr(ithet+ndim*(i-1)+j-1)*der(k)
                    dfdm(j, k) = dfdm(j, k)+fno(ndim*(i-1)+j)*der(k)
                end do
                dudm(j, 4) = dudm(j, 4)+zr(idepl+ndim*(i-1)+j-1)*der(4)
                dtdm(j, 4) = dtdm(j, 4)+zr(ithet+ndim*(i-1)+j-1)*der(4)
                dfdm(j, 4) = dfdm(j, 4)+fno(ndim*(i-1)+j)*der(4)
            end do
        end do
!
        ttrg = tpg(kp)-tref
        ttrgv = 0.d0
!
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
!
! ----- RECUPERATION DE E, NU, ALPHA ET RHO
!
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                    1, r8bid, iret)
        call rcvalb(fami, kp, 1, '+', zi(imate), &
                    ' ', phenom, 0, ' ', [0.d0], &
                    4, nomres, valres, icodre, 0)
        ASSERT(icodre(1)+icodre(2) .eq. 0)
        if (icodre(3) .ne. 0) then
            ASSERT(iret .ne. 0)
            valres(3) = 0.d0
        end if
        if (icodre(4) .ne. 0) then
            valres(4) = 0.d0
        end if
!
        e = valres(1)
        nu = valres(2)
        alpha = valres(3)
        rho = valres(4)
        k3a = alpha*e/(1.d0-2.d0*nu)
!
        la = nu*e/((1.d0+nu)*(1.d0-2.d0*nu))
        mu = e/(2.d0*(1.d0+nu))
!       EN DP
        ka = 3.d0-4.d0*nu
        coeff = e/(1.d0-nu*nu)
        coeff3 = 2.d0*mu
!       EN CP
!       KA=(3.D0-NU)/(1.D0+NU)
!       COEFF=E
!       COEFF3=2.D0 * MU
!
        c1 = la+2.d0*mu
        c2 = la
        c3 = mu
!
! ---- DETERMINATION DES GRANDEURS UTILES A LA PRISE EN COMPTE DE LA CONTRAINTE INITIALE
! ---- CONTRAINTE INITIALE, SA DERIVEE ET LA DEFORMATION ASSOCIEE EPSREF
!
        if (isigi .ne. 0) then
            do i = 1, nno
                der(1) = dfdi(i)
                der(2) = dfdi(i+nno)
                der(3) = dfdi(i+2*nno)
                der(4) = zr(ivf+l+i-1)
! CALCUL DE SIGMA INITIAL
!
                do j = 1, ncmp
                    sigin(j) = sigin(j)+zr(isigi+ncmp*(i-1)+j-1)*der(4)
                end do
!
! CALCUL DU GRADIENT DE SIGMA INITIAL
                do j = 1, ncmp
                    do k = 1, ndim
                        dsigin(j, k) = dsigin(j, k)+zr(isigi+ncmp*(i-1)+j-1)*der(k)
                    end do
                end do
            end do
!
! TRAITEMENTS PARTICULIERS DES TERMES CROISES
            do i = 4, ncmp
                sigin(i) = sigin(i)*rac2
                do j = 1, ndim
                    dsigin(i, j) = dsigin(i, j)*rac2
                end do
            end do
!
!
!
! CALCUL DE LA DEFORMATION DE REFERENCE
!
!
            epsref(1) = -(1.d0/e)*(sigin(1)-(nu*(sigin(2)+sigin(3))))
            epsref(2) = -(1.d0/e)*(sigin(2)-(nu*(sigin(3)+sigin(1))))
            epsref(3) = -(1.d0/e)*(sigin(3)-(nu*(sigin(1)+sigin(2))))
            epsref(4) = -(1.d0/mu)*sigin(4)
            epsref(5) = -(1.d0/mu)*sigin(5)
            epsref(6) = -(1.d0/mu)*sigin(6)
!
        end if
!
!
! ----- CALCUL DES CHAMPS AUXILIAIRES ET DE LEURS DERIVEES
!
        do ino = 1, nno
            ffp(ino) = zr(ivf-1+nno*(kp-1)+ino)
        end do
        p(:, :) = 0.d0
        invp(:, :) = 0.d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCUL DES COOR. CYL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call coor_cyl(ndim, nno, zr(ibalo), zr(igeom), ffp, &
                      p, invp, rg, phig, l_not_zero)
!
!       PRISE EN COMPTE DE LA COURBURE : OUI
!
        lcour = .true.
!       RECUPERATION DU TENSEUR DE COURBURE
        call jevech('PCOURB', 'L', icour)
        do i = 1, ndim
            do j = 1, ndim
                courb(i, 1, j) = zr(icour-1+ndim*(i-1)+j)
                courb(i, 2, j) = zr(icour-1+ndim*(i+3-1)+j)
                courb(i, 3, j) = zr(icour-1+ndim*(i+6-1)+j)
            end do
        end do
!
!
        call chauxi(ndim, mu, ka, rg, phig, &
                    invp, lcour, courb, du1dm, du2dm, &
                    du3dm, u1l, u2l, u3l)
!
!
!
!
!-----------------------------------------------------------------------
!       CALCUL DE G, K1, K2, K3 AU POINT DE GAUSS
!-----------------------------------------------------------------------
!
        guv = 0.d0
        coef = 2.d0
        call gbil3d(dudm, dudm, dtdm, dfdm, dfdm, &
                    tgdm, tgdm, ttrg, ttrg, poids, &
                    sigin, dsigin, epsref, c1, c2, &
                    c3, k3a, alpha, coef, rho, &
                    puls, guv)
        g = g+guv
!
        guv1 = 0.d0
        coef = 1.d0
        call gbil3d(dudm, du1dm, dtdm, dfdm, dfvdm, &
                    tgdm, tgvdm, ttrg, ttrgv, poids, &
                    sigin, dsigin, epsref, c1, c2, &
                    c3, k3a, alpha, coef, rho, &
                    puls, guv1)
        k1 = k1+guv1
!
        guv2 = 0.d0
        coef = 1.d0
        call gbil3d(dudm, du2dm, dtdm, dfdm, dfvdm, &
                    tgdm, tgvdm, ttrg, ttrgv, poids, &
                    sigin, dsigin, epsref, c1, c2, &
                    c3, k3a, alpha, coef, rho, &
                    puls, guv2)
        k2 = k2+guv2
!
        guv3 = 0.d0
        coef = 1.d0
        call gbil3d(dudm, du3dm, dtdm, dfdm, dfvdm, &
                    tgdm, tgvdm, ttrg, ttrgv, poids, &
                    sigin, dsigin, epsref, c1, c2, &
                    c3, k3a, alpha, coef, rho, &
                    puls, guv3)
        k3 = k3+guv3
!
    end do
!
    k1 = k1*coeff
    k2 = k2*coeff
    k3 = k3*coeff3
!
    zr(ificg) = g
    zr(ificg+1) = k1/sqrt(coeff)
    zr(ificg+2) = k2/sqrt(coeff)
    zr(ificg+3) = k3/sqrt(coeff3)
    zr(ificg+4) = k1
    zr(ificg+5) = k2
    zr(ificg+6) = k3
!
999 continue
!
!
end subroutine
