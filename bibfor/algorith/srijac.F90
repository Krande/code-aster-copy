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

subroutine srijac(nmat, materf, timed, timef, &
                  yf, deps, nr, nvi, vind, &
                  vinf, yd, drdy)

!

!!!
!!! MODELE LKR : CALCUL DU JACOBIEN DE LKR = DRDY(DY)
!!!

! ===================================================================================
! IN  : MOD               : TYPE DE MODELISATION
!     : NMAT              : DIMENSION MATER
!     : MATERF(NMAT,2)    : COEFFICIENTS MATERIAU A T+DT
!     : YF(NDT+3)         : VARIABLES A T + DT =  ( SIGF DLAMBDA XI_P XI_VP)
!     : DEPS(6)           : INCREMENT DE DEFORMATION
!     : TIMED             : INSTANT  T
!     : TIMEF             : INSTANT  T+DT
!     : NR                : DIMENSION DECLAREE DRDY
!     : NVI               : NOMBRE DE VARIABLES INTERNES
!     : VIND(NVI)         : VARIABLE INTERNES A T
!     : VINF(NVI)         : VARIABLE INTERNES A T+DT
!     : YD(NDT+3)         : VARIABLES A T  = ( SIGD  0 XI_P XI_VP) A T
!     : DY(NDT+3)         : SOLUTION = ( DSIG  DLAMBDA  DXI_P DXI_VP )
! OUT : DRDY(NDT+3,NDT+3) : JACOBIEN DU SYSTEME NON LINEAIRE
!     : IRET              : CODE RETOUR
! ===================================================================================

    implicit none

#include "asterc/r8prem.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcprte.h"
#include "asterfort/srbpri.h"
#include "asterfort/srcalg.h"
#include "asterfort/srcaln.h"
#include "asterfort/srcrip.h"
#include "asterfort/srcriv.h"
#include "asterfort/srdepp.h"
#include "asterfort/srdfds.h"
#include "asterfort/srdfdx.h"
#include "asterfort/srdgde.h"
#include "asterfort/srdgds.h"
#include "asterfort/srdhds.h"
#include "asterfort/srdndx.h"
#include "asterfort/srds2h.h"
#include "asterfort/srelas.h"
#include "asterfort/srfsxi.h"
#include "asterfort/srvacp.h"
#include "asterfort/srvacv.h"
#include "asterfort/srvarp.h"
#include "asterfort/srvarv.h"
#include "asterf_types.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nr, nmat, nvi
    real(kind=8) :: deps(6), drdy(nr, nr), yf(nr), yd(nr)
    real(kind=8) :: materf(nmat, 2)
    real(kind=8) :: timed, timef, vind(nvi), vinf(nvi)

    !!!
    !!! Varibales locales
    !!!

    integer(kind=8) :: i, j, varv, valv, valp, retcom, ndt, ndi

    real(kind=8) :: sigft(6), depst(6)
    real(kind=8) :: devgii, vint(nvi), dt, devsig(6), i1
    real(kind=8) :: xi5, xi1, dsdenl(6, 6), kk, mu
    real(kind=8) :: ucriv, seuilv, depsv(6), dgamv
    real(kind=8) :: seuilp, ucrip, seuivm, dhds(6), ds2hds(6)
    real(kind=8) :: paraep(3), varpl(4), dfdsp(6), bprimp
    real(kind=8) :: vecnp(6), gp(6), vetemp(6), derpar(3), dfdxip
    real(kind=8) :: depse(6), hook(6, 6), mue, ke
    real(kind=8) :: dsige(6), vident(6), dhokds(6, 6), patm, nelas
    real(kind=8) :: paravi(3), varavi(4), dfvdsi(6)
    real(kind=8) :: dgpds(6, 6), dgvds(6, 6)
    real(kind=8) :: dldgds(6, 6), dlambd, hnldgp(6, 6), dsgvds(6, 6)
    real(kind=8) :: hnldgv(6, 6), bprimv, vecnv(6), ucrim, devgiv, gv(6)
    real(kind=8) :: hnlgv(6), hnldfg(6, 6), phiv, av, nv, dphiv
    real(kind=8) :: dphvds(6), dr1dy3(6), dgpdxi(6), dfsdxp(6)
    real(kind=8) :: term1, term2, term3, coupl
    real(kind=8) :: dndxip(6), dfsdxv(6), dpadxp(3)
    real(kind=8) :: dndxiv(6), dfdxiv, dpadxv(3), dphidx
    real(kind=8) :: dphdxg(6), dgvdxi(6), phdgdx(6)
    real(kind=8) :: dr1dy4(6), dgipds(6), dgivds(6), dgipdx, dgivdx
    real(kind=8) :: mident(6, 6), kron(6), kron2(6, 6), unstro
    real(kind=8) :: dsdsig(6, 6), dgtvds(6, 6), dgtpds(6, 6), kron3(6, 6)
    real(kind=8) :: devgp(6), devgv(6), dgtpdx(6), dgtvdx(6), dxiv
    real(kind=8) :: z, r, a0, tpp, trr, dtmp, xi2, xi20, rx2
    real(kind=8) :: xi10, xi50, rx1, rx5

    aster_logical :: plas

    common/tdim/ndt, ndi

    !!!
    !!! Init
    !!!

    sigft = 0.d0
    depst = 0.d0
    devgii = 0.d0
    vint = 0.d0
    dt = 0.d0
    devsig = 0.d0
    i1 = 0.d0
    xi5 = 0.d0
    xi1 = 0.d0
    dsdenl = 0.d0
    kk = 0.d0
    mu = 0.d0
    ucriv = 0.d0
    seuilv = 0.d0
    depsv = 0.d0
    dgamv = 0.d0
    seuilp = 0.d0
    ucrip = 0.d0
    seuivm = 0.d0
    dhds = 0.d0
    ds2hds = 0.d0
    paraep = 0.d0
    varpl = 0.d0
    dfdsp = 0.d0
    bprimp = 0.d0
    vecnp = 0.d0
    gp = 0.d0
    vetemp = 0.d0
    derpar = 0.d0
    dfdxip = 0.d0
    depse = 0.d0
    hook = 0.d0
    mue = 0.d0
    ke = 0.d0
    dsige = 0.d0
    vident = 0.d0
    dhokds = 0.d0
    patm = 0.d0
    nelas = 0.d0
    paravi = 0.d0
    varavi = 0.d0
    dfvdsi = 0.d0
    dgpds = 0.d0
    dgvds = 0.d0
    dldgds = 0.d0
    dlambd = 0.d0
    hnldgp = 0.d0
    dsgvds = 0.d0
    hnldgv = 0.d0
    bprimv = 0.d0
    vecnv = 0.d0
    ucrim = 0.d0
    devgiv = 0.d0
    gv = 0.d0
    hnlgv = 0.d0
    hnldfg = 0.d0
    phiv = 0.d0
    av = 0.d0
    nv = 0.d0
    dphiv = 0.d0
    dphvds = 0.d0
    dr1dy3 = 0.d0
    dgpdxi = 0.d0
    dfsdxp = 0.d0
    term1 = 0.d0
    term2 = 0.d0
    term3 = 0.d0
    coupl = 0.d0
    dndxip = 0.d0
    dfsdxv = 0.d0
    dpadxp = 0.d0
    dndxiv = 0.d0
    dfdxiv = 0.d0
    dpadxv = 0.d0
    dphidx = 0.d0
    dphdxg = 0.d0
    dgvdxi = 0.d0
    phdgdx = 0.d0
    dr1dy4 = 0.d0
    dgipds = 0.d0
    dgivds = 0.d0
    dgipdx = 0.d0
    dgivdx = 0.d0
    mident = 0.d0
    kron = 0.d0
    kron2 = 0.d0
    unstro = 0.d0
    dsdsig = 0.d0
    dgtvds = 0.d0
    dgtpds = 0.d0
    kron3 = 0.d0
    devgp = 0.d0
    devgv = 0.d0
    dgtpdx = 0.d0
    dgtvdx = 0.d0
    dxiv = 0.d0
    z = 0.d0
    r = 0.d0
    a0 = 0.d0
    tpp = 0.d0
    trr = 0.d0
    dtmp = 0.d0
    xi2 = 0.d0
    xi20 = 0.d0
    rx2 = 0.d0
    xi10 = 0.d0
    xi50 = 0.d0
    rx1 = 0.d0
    rx5 = 0.d0
    plas = .false.

    !!!
    !!! Recup. des temperatures
    !!!

    tpp = materf(7, 1)
    trr = materf(8, 1)

    if ((tpp .ge. trr) .and. (trr .gt. 0.d0)) then
        dtmp = tpp-trr
    else
        dtmp = 0.d0
    end if

    !!!
    !!! Passage en convention meca. des sols
    !!!

    do i = 1, ndt
        sigft(i) = -yf(i)
        depst(i) = -deps(i)
    end do

    !!!
    !!! Variables locales tmp
    !!!

    mident(:, :) = 0.d0

    do i = 1, ndt
        mident(i, i) = 1.d0
    end do

    varv = 0
    devgii = 0.d0
    dlambd = yf(ndt+1)

    !!!
    !!! Vecteur variables internes tmp
    !!!

    vint(1:nvi) = vind(1:nvi)

    if (yf(ndt+2) .ge. vind(1)) then
        vint(1) = yf(ndt+2)
    else
        vint(1) = vind(1)
    end if

    if (yf(ndt+3) .ge. vind(3)) then
        vint(3) = yf(ndt+3)
    else
        vint(3) = vind(3)
    end if

    !!! Increment de temps
    dt = timef-timed

    !!!
    !!! Construction deviateur des contraintes et 1er invariant
    !!!

    call lcdevi(sigft, devsig)

    i1 = sigft(1)+sigft(2)+sigft(3)

    !!!
    !!! Recup. des para. materiau
    !!!

    xi10 = materf(12, 2)
    xi20 = materf(13, 2)
    xi50 = materf(14, 2)
    rx1 = materf(24, 2)
    rx2 = materf(25, 2)
    rx5 = materf(26, 2)
    xi5 = xi50*exp(rx5*dtmp)
    xi1 = xi10*exp(rx1*dtmp)
    xi2 = xi20*exp(rx2*dtmp)

    !!!
    !!! Construction tenseur elastique non lineaire
    !!!

    call srelas(ndi, ndt, nmat, materf, sigft, dsdenl, kk, mu)

    !!!
    !!! 1) Calcul de la def. visco. et du para. d'ecrouissage
    !!!

    !!! 1-1) Indicateur angle de dilatance visco.
    valv = 0

    !!! 1-2) Calcul seuil visco. par rapport a yf
    seuilv = 0.d0

    call srcriv(vint(3), i1, devsig, nmat, materf, tpp, ucriv, seuilv)

    if (seuilv .ge. 0.d0) then

        call srdgde(valv, vint(3), dt, seuilv, ucriv, &
                    i1, devsig, vint, nvi, nmat, materf, &
                    tpp, depsv, dgamv, retcom)

    else

        dgamv = 0.d0

        do i = 1, ndt
            depsv(i) = 0.d0
        end do

        seuilv = 0.d0
        ucriv = 0.d0

    end if

    !!!
    !!! 2) Calcul de depsp et dgamp
    !!!

    call srdhds(nmat, materf, devsig, dhds, retcom)
    call srds2h(nmat, materf, devsig, dhds, ds2hds, retcom)

    !!! 2-1) Calcul fonction seuil en yf
    seuilp = 0.d0

    call srcrip(i1, devsig, vint, nvi, nmat, materf, tpp, ucrip, seuilp)

    !!! 2-2) Indicateur contractance ou dilatance
    seuivm = 0.d0

    call srcriv(xi5, i1, devsig, nmat, materf, tpp, ucrim, seuivm)

    if (seuivm .lt. 0.d0) then
        varv = 0
    else
        varv = 1
    end if

    !!! 2-3) Si seuilp >= 0, alors plasticite
    if ((seuilp .ge. 0.d0) .or. (vinf(7) .gt. 0.d0)) then

        !!! 2-3-1) Indicateur angle de dilatance
        if (yf(ndt+2) .lt. xi1) then
            valp = 0
        else
            valp = 1
        end if

        !!! 2-3-2) Calcul de df/dsig
        call srvarp(vint, nvi, nmat, materf, tpp, paraep)
        call srvacp(nmat, materf, paraep, varpl)
        call srdepp(vint, nvi, nmat, materf, paraep, derpar)
        call srdfds(nmat, materf, paraep, varpl, ds2hds, ucrip, dfdsp)

        !!! 2-3-2) Calcul de gp
        bprimp = srbpri(valp, vint, nvi, nmat, materf, paraep, i1, devsig, tpp)

        vecnp(:) = 0.d0
        call srcaln(devsig, bprimp, vecnp, retcom)
        call srcalg(dfdsp, vecnp, gp, devgii)

        !!! 2-3-4) Calcul def. elastique
        do i = 1, ndt
            depse(i) = depst(i)-yf(ndt+1)*gp(i)-depsv(i)
        end do

        !!! 2-3-5) calcul de dgp/dsigma
        call srdgds(nmat, materf, paraep, varpl, devsig, &
                    i1, valp, ds2hds, vecnp, dfdsp, &
                    bprimp, nvi, vint, dhds, tpp, dgpds)

        !!! 2-3-6) Produit matriciel hook*dlambda*d(gp)/d(sig)
        dldgds(1:ndt, 1:ndt) = dlambd*dgpds(1:ndt, 1:ndt)
        hnldgp(1:ndt, 1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), dldgds(1:ndt, 1:ndt))

        !!! 2-3-7) Calcul de d2(fp)/d(sig)d(xi)
        plas = .true.
        call srfsxi(nmat, materf, i1, devsig, ds2hds, &
                    plas, vint(1), paraep, varpl, tpp, dfsdxp, dpadxp)

        !!! 2-3-8) Calcul de d(n)/d(xi)
        call srdndx(nmat, materf, i1, devsig, bprimp, &
                    valp, paraep, vint(1), tpp, derpar, dndxip)

    !!! 2-4) Pas de plasticite
    else

        do i = 1, ndt
            depse(i) = depst(i)-depsv(i)
        end do

        hnldgp(:, :) = 0.d0
        dgpds(:, :) = 0.d0
        dfdsp(:) = 0.d0
        gp(:) = 0.d0
        vecnp(:) = 0.d0
        dfsdxp(:) = 0.d0
        dndxip(:) = 0.d0

        devgii = 0.d0

    end if

    !!!
    !!! Calcul de d(r1)/d(y)
    !!!

    !!! Calcul de d(r1)/d(y1) : y1 == sigma

    !!! construction du tenseur elasticite lineaire

    mue = materf(4, 1)
    ke = materf(5, 1)

    hook(:, :) = 0.d0

    do i = 1, ndi
        do j = 1, ndi
            hook(i, j) = ke-2.d0*mue/3.d0
        end do
    end do

    do i = 1, ndt
        hook(i, i) = hook(i, i)+2.d0*mue
    end do

    !!! contrainte elastique
    dsige(1:ndt) = matmul(hook(1:ndt, 1:ndt), depse(1:ndt))

    !!! produit tensoriel d(sige) x vident
    patm = materf(1, 2)
    nelas = materf(2, 2)

    vident(:) = 0.d0

    do i = 1, ndi
        vident(i) = nelas/3.d0/patm*(i1/(3.d0*patm))**(nelas-1.d0)
    end do

    call lcprte(dsige, vident, dhokds)

    !!! calcul de d(fv)/d(sig)
    call srvarv(vint(3), nmat, materf, tpp, paravi)
    call srvacv(nmat, materf, paravi, varavi)

    bprimv = srbpri(valv, vint, nvi, nmat, materf, paravi, i1, devsig, tpp)

    call srcaln(devsig, bprimv, vecnv, retcom)
    call srdfds(nmat, materf, paravi, varavi, ds2hds, ucriv, dfvdsi)

    !!! calcul de gvp
    call srcalg(dfvdsi, vecnv, gv, devgiv)

    !!! calcul de d(gvp)/d(sig)
    call srdgds(nmat, materf, paravi, varavi, devsig, &
                i1, valv, ds2hds, vecnv, dfvdsi, &
                bprimv, nvi, vint, dhds, tpp, dgvds)

    !!! produit matriciel de hook * phi_v * d(gv)/d(sig)
    r = 8.3144621d0
    a0 = materf(16, 2)
    nv = materf(17, 2)
    z = materf(27, 2)

    if ((tpp .ge. trr) .and. (trr .gt. 0.d0)) then
        av = a0*exp(-z/r/tpp*(1.d0-tpp/trr))
    else
        av = a0
    end if

    phiv = av*(seuilv/patm)**nv

    dsgvds(1:ndt, 1:ndt) = phiv*dgvds(1:ndt, 1:ndt)
    hnldgv(1:ndt, 1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), dsgvds(1:ndt, 1:ndt))

    !!! produit matriciel hook*d(phiv)/d(sig)*gv
    dphiv = av*nv/patm*(seuilv/patm)**(nv-1.d0)

    dphvds(1:ndt) = dphiv*dfvdsi(1:ndt)
    hnlgv(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), gv(1:ndt))
    call lcprte(hnlgv, dphvds, hnldfg)

    !!! assemblage
    do i = 1, ndt
        do j = 1, ndt
            drdy(i, j) = -(mident(i, j)-dhokds(i, j)+hnldgp(i, j)+hnldgv(i, j)*dt+&
                    & hnldfg(i, j)*dt)/mu
        end do
    end do

    !!! Calcul de d(r1)/d(y2) - y2 == dlambda

    if (vinf(7) .le. 0.d0) then
        do i = 1, ndt
            drdy(i, ndt+1) = 0.d0
        end do
    else
        vetemp(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), gp(1:ndt))
        do i = 1, ndt
            drdy(i, ndt+1) = vetemp(i)/mu
        end do
    end if

    !!! Calcul de d(r1)/d(y3) - y3 == xip

    !!! assemblage d(gp)/d(xi)
    term1 = dot_product(dfsdxp(1:ndt), vecnp(1:ndt))
    term2 = dot_product(dfdsp(1:ndt), dndxip(1:ndt))
    term3 = dot_product(dfdsp(1:ndt), vecnp(1:ndt))

    do i = 1, ndt
        dgpdxi(i) = dfsdxp(i)-term1*vecnp(i)-term2*vecnp(i)-term3*dndxip(i)
    end do

    !!! assemblage final
    dr1dy3(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), dgpdxi(1:ndt))
    dr1dy4(1:ndt) = dlambd*dr1dy3(1:ndt)

    do i = 1, ndt
        drdy(i, ndt+2) = dr1dy4(i)/mu
    end do

    !!! Calcul de d(r1)/d(y4) - y4 == xivp

    dxiv = min(dgamv, xi5-yd(ndt+3))
    if (abs(dxiv-dgamv) .lt. r8prem()) then
        !!! calcul de d(fvp)/d(sig)d(xiv)
        plas = .false.

        call srfsxi(nmat, materf, i1, devsig, ds2hds, &
                    plas, vint(3), paravi, varavi, tpp, dfsdxv, dpadxv)

        !!! calcull de d(n)/d(xi)
        call srdndx(nmat, materf, i1, devsig, bprimv, &
                    valv, paravi, vint(3), tpp, dpadxv, dndxiv)

        !!! assemblage de d(gv)/d(xiv)
        term1 = dot_product(dfsdxv(1:ndt), vecnv(1:ndt))
        term2 = dot_product(dfvdsi(1:ndt), dndxiv(1:ndt))
        term3 = dot_product(dfvdsi(1:ndt), vecnv(1:ndt))

        do i = 1, ndt
            dgvdxi(i) = dfsdxv(i)-term1*vecnv(i)-term2*vecnv(i)-term3*dndxiv(i)
        end do

        !!! calcul de d(phiv)/d(xiv)
        call srdfdx(nmat, materf, ucriv, i1, devsig, &
                    paravi, varavi, dpadxv, dfdxiv)

        !!! assemblage d(r1)/d(y4)
        dphidx = dphiv*dfdxiv

        dphdxg(1:ndt) = dphidx*gv(1:ndt)
        phdgdx(1:ndt) = phiv*dgvdxi(1:ndt)
        vetemp(1:ndt) = dphdxg(1:ndt)+phdgdx(1:ndt)
        dr1dy4(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), vetemp(1:ndt))

        do i = 1, ndt
            drdy(i, ndt+3) = dr1dy4(i)/mu*dt
        end do

    else

        do i = 1, ndt
            drdy(i, ndt+3) = 0.d0
        end do

    end if

    !!!
    !!! Calcul de d(r2)/d(y)
    !!!

    !!! Application de la condition de kt sur r(ndt+1)

    if (vinf(7) .le. 0.d0) then

        !!! calcul de d(r2)/d(y1) - y1 == sigma
        do i = 1, ndt
            drdy(ndt+1, i) = 0.d0
        end do

        !!! calcul de d(r2)/d(y2) - y2 == dlambda
        drdy(ndt+1, ndt+1) = 1.d0

        !!! calcul de d(r2)/d(y3) - y3 == xip
        drdy(ndt+1, ndt+2) = 0.d0

    else

        !!! calcul de d(r2)/d(y1) - y1 == sigma
        do i = 1, ndt
            drdy(ndt+1, i) = -dfdsp(i)/mu
        end do

        !!! calcul de d(r2)/d(y2) - y2 == dlambda
        drdy(ndt+1, ndt+1) = 0.d0

        !!! calcul de d(r2)/d(y3) - y3 = xip
        call srdfdx(nmat, materf, ucrip, i1, devsig, paraep, varpl, derpar, dfdxip)

        drdy(ndt+1, ndt+2) = dfdxip/mu

    end if

    !!! Calcul de d(r2)/d(y4) - y4 = xivp
    drdy(ndt+1, ndt+3) = 0.d0

    !!!
    !!! Calcul de d(r3)/d(y)
    !!!

    !!! Calcul de d(r3)/d(y1) - y1 == sigma

    kron(:) = 0.d0

    do i = 1, ndi
        kron(i) = 1.d0
    end do

    !!! construction de d(s)/d(sig)
    unstro = 1.d0/3.d0

    call lcprte(kron, kron, kron2)
    kron3(1:ndt, 1:ndt) = unstro*kron2(1:ndt, 1:ndt)
    dsdsig(1:ndt, 1:ndt) = mident(1:ndt, 1:ndt)-kron3(1:ndt, 1:ndt)

    !!! construction de dev(g)
    call lcdevi(gv, devgv)
    call lcdevi(gp, devgp)

    !!! construction de d(devgii)/d(sig)
    dgtvds(1:ndt, 1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), dgvds(1:ndt, 1:ndt))
    dgtpds(1:ndt, 1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), dgpds(1:ndt, 1:ndt))
    dgivds(:) = 0.d0
    dgipds(:) = 0.d0

    if ((seuilp .ge. 0.d0) .or. (vinf(7) .gt. 0.d0)) then
        do i = 1, ndt
            do j = 1, ndt
                dgivds(i) = dgivds(i)+devgv(j)/devgiv*dgtvds(j, i)
                dgipds(i) = dgipds(i)+devgp(j)/devgii*dgtpds(j, i)
            end do
        end do
    else
        do i = 1, ndt
            do j = 1, ndt
                dgivds(i) = dgivds(i)+devgv(j)/devgiv*dgtvds(j, i)
            end do
        end do
    end if

    coupl = materf(28, 2)

    if ((varv .eq. 1) .and. (coupl .ge. 1.d0/2.d0)) then
        do i = 1, ndt
            drdy(ndt+2, i) = sqrt(2.d0/3.d0)*(dlambd*dgipds(i)+ &
                        & (dphvds(i)*devgiv+phiv*dgivds(i))*dt)
        end do
    else
        do i = 1, ndt
            drdy(ndt+2, i) = dlambd*sqrt(2.d0/3.d0)*dgipds(i)
        end do
    end if

    !!! Calcul de d(r3)/d(y2) - y2 == dlambda

    !!! application de la condition de kt
    if (vinf(7) .le. 0.d0) then
        drdy(ndt+2, ndt+1) = 0.d0
    else
        drdy(ndt+2, ndt+1) = -devgii*sqrt(2.d0/3.d0)
    end if

    !!! Calcul de d(r3)/d(y3) - y3 == xip

    dgtpdx(1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), dgpdxi(1:ndt))
    dgipdx = dot_product(devgp(1:ndt), dgtpdx(1:ndt))

    if (vinf(7) .gt. 0.d0) then
        drdy(ndt+2, ndt+2) = 1.d0-dlambd*sqrt(2.d0/3.d0)*dgipdx/devgii
    else
        drdy(ndt+2, ndt+2) = 1.d0
    end if

    !!! Calcul de d(r3)/d(y4) - y4 == xivp

    dgtvdx(1:ndt) = matmul(dsdsig(1:ndt, 1:ndt), dgvdxi(1:ndt))
    dgivdx = dot_product(devgv(1:ndt), dgtvdx(1:ndt))

    if (abs(dxiv-dgamv) .lt. r8prem()) then
        drdy(ndt+2, ndt+3) = -(dphidx*devgiv+phiv*dgivdx/devgiv)*sqrt(2.d0/3.d0)*dt
    else
        drdy(ndt+2, ndt+3) = 0.d0
    end if

    !!!
    !!! Calcul de d(r4)/d(y)
    !!!

    if (abs(dxiv-dgamv) .lt. r8prem()) then

        !!! Calcul de d(r4)/d(y1) - y1 == sigma
        do i = 1, ndt
            drdy(ndt+3, i) = (dphvds(i)*devgiv+phiv*dgivds(i))*sqrt(2.d0/3.d0)*dt
        end do

        !!! Calcul de d(r4)/d(y2) - y2 == dlambda
        drdy(ndt+3, ndt+1) = 0.d0

        !!! Calcul de d(r4)/d(y3) - y3 == xip
        drdy(ndt+3, ndt+2) = 0.d0

        !!! Calcul de d(r4)/d(y4) - y4 == xivp
        drdy(ndt+3, ndt+3) = 1.d0-sqrt(2.d0/3.d0)*dt*(dphidx*devgiv+phiv*dgivdx/devgiv)

    else

        !!! Calcul de d(r4)/d(y1) - y1 == sigma
        do i = 1, ndt
            drdy(ndt+3, i) = 0.d0
        end do

        !!! Calcul de d(r4)/d(y2) - y2 == dlambda
        drdy(ndt+3, ndt+1) = 0.d0

        !!! Calcul de d(r4)/d(y3) - y3 == xip
        drdy(ndt+3, ndt+2) = 0.d0

        !!! Calcul de d(r4)/d(y4) - y4 == xivp
        drdy(ndt+3, ndt+3) = 1.d0

    end if

end subroutine
