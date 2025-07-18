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

subroutine srresi(nmat, materf, timed, timef, &
                  nvi, vind, vinf, yd, yf, deps, nr, r)

!

!!!
!!! MODELE LKR : CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = -R(DY)
!!!

! ===================================================================================
! IN  : TYPMOD         : TYPE DE MODELISATION
!     : NMAT           : DIMENSION MATER
!     : MATERF(NMAT,2) : COEFFICIENTS MATERIAU A T+DT
!     : TIMED          : INSTANT  T
!     : TIMEF          : INSTANT  T+DT
!     : NVI            : NOMBRE DE VARIABLES INTERNES
!     : DEPS(6)        : INCREMENT DE DEFORMATION
!     : VIND(NVI)      : VARIABLES INTERNES A T
!     : VINF(NVI)      : VARIABLES INTERNES A T+DT
!     : YD(NDT+3)      : VARIABLES A T    = ( SIGD 0    XIPD XIVPD (EPSD3))
!     : YF(NDT+3)      : VARIABLES A T+DT = ( SIGF DLAM XIPF XIVPF (EPS3F))
!     : DEPS(6)        : INCREMENT DE DEFORMATIONS
!     : DY(NDT+3)      : SOLUTION         = ( DSIG DLAM DXIP DXIVP (DEPS3))
!     : NR(NDT+3)      : DIMENSION DU VECTEUR INCONNUES
! OUT : R              : SYSTEME NL A T+DT
! ===================================================================================

    implicit none

#include "asterfort/lcdevi.h"
#include "asterfort/srbpri.h"
#include "asterfort/srcalg.h"
#include "asterfort/srcaln.h"
#include "asterfort/srcrip.h"
#include "asterfort/srcriv.h"
#include "asterfort/srdfds.h"
#include "asterfort/srdgde.h"
#include "asterfort/srdhds.h"
#include "asterfort/srds2h.h"
#include "asterfort/srelas.h"
#include "asterfort/srvacp.h"
#include "asterfort/srvarp.h"

    !!!
    !!! Variables globales
    !!!

    integer(kind=8) :: nmat, nr, nvi, ndi, ndt
    real(kind=8) :: deps(6), vind(nvi), vinf(nvi), r(nr), yd(nr), yf(nr), materf(nmat, 2)
    real(kind=8) :: timed, timef

    !!!
    !!! Variables locales
    !!!

    integer(kind=8) :: i, retcom, val, varv

    real(kind=8) :: vint(nvi), devsig(6), i1, ucrip, seuilp
    real(kind=8) :: dt, seuilv, depsv(6), dgamv
    real(kind=8) :: dxiv, xi5, xi1, xi2, seuivm, ucriv, coupl
    real(kind=8) :: dsdenl(6, 6), kk, mu, dhds(6), ds2hds(6)
    real(kind=8) :: paraep(3), varpl(4), dfdsp(6), bprimp
    real(kind=8) :: vecnp(6), gp(6), devgii, depse(6)
    real(kind=8) :: dsige(6), sigdt(6), sigft(6), depst(6), lamgd2
    real(kind=8) :: xi10, xi20, xi50, rx1, rx2, rx5, trr, tpp, dtmp

    common/tdim/ndt, ndi

    !!!
    !!! Redcup. des temp.
    !!!

    tpp = materf(7, 1)
    trr = materf(8, 1)

    if ((tpp .gt. trr) .and. (trr .gt. 0.d0)) then
        dtmp = tpp-trr
    else
        dtmp = 0.d0
    end if

    !!!
    !!! Passage en convention mecanique des sols
    !!!

    do i = 1, ndt
        sigft(i) = -yf(i)
        sigdt(i) = -yd(i)
        depst(i) = -deps(i)
    end do

    !!!
    !!! Variable tmp.
    !!!

    varv = 0
    val = 0
    devgii = 0.d0

    !!!
    !!! Variables itnernes tmp
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

    !!!! Increment de temps
    dt = timef-timed

    !!!
    !!! Tenseur deviatoire et premier invariant
    !!!

    call lcdevi(sigft, devsig)
    i1 = sigft(1)+sigft(2)+sigft(3)

    !!!
    !!! Parametres materiaux
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
    !!! 1) Calcul de la def. visco. depsv et de dgamv
    !!!

    !!! 1-1) Indicateur sur angle de dilatance visqueux psi - val = 0 pour visco.
    val = 0

    !!! 1-2) Calcul du seuil visco. par rapport a yf
    call srcriv(vint(3), i1, devsig, nmat, materf, tpp, ucriv, seuilv)

    !!! 1-3) Si seuil >= 0, appel a srdgde pour calcul de depsv et dgamv
    if (seuilv .ge. 0.d0) then

        call srdgde(val, vint(3), dt, seuilv, ucriv, i1, devsig, vint, nvi, nmat, &
                    materf, tpp, depsv, dgamv, retcom)

    else

        dgamv = 0.d0

        do i = 1, ndt
            depsv(i) = 0.d0
        end do

    end if

    !!!
    !!! 2) Calcul de la def. plastique depsp et de dgamp
    !!!

    !!! 2-1) Calcul de la position de sigma / seuil cara. pour varv
    call srcriv(xi5, i1, devsig, nmat, materf, tpp, ucriv, seuivm)

    if (seuivm .le. 0.d0) then
        varv = 0
    else
        varv = 1
    end if

    !!! 2-2) Calcul de la fonction seuil plastique en yf
    seuilp = 0.d0

    call srcrip(i1, devsig, vint, nvi, nmat, materf, tpp, ucrip, seuilp)

    !!! 2-2-1) Si seuilp >= 0, plasticite
    if ((seuilp .ge. 0.d0) .or. (vinf(7) .gt. 0.d0)) then
        !!! indicateur angle de dil.
        if (yf(ndt+2) .lt. xi1) then
            val = 0
        else
            val = 1
        end if

        !!! calcul de df/dsig
        call srdhds(nmat, materf, devsig, dhds, retcom)
        call srds2h(nmat, materf, devsig, dhds, ds2hds, retcom)
        call srvarp(vint, nvi, nmat, materf, tpp, paraep)
        call srvacp(nmat, materf, paraep, varpl)
        call srdfds(nmat, materf, paraep, varpl, ds2hds, ucrip, dfdsp)

        !!! calcl de g
        bprimp = srbpri(val, vint, nvi, nmat, materf, paraep, i1, devsig, tpp)

        call srcaln(devsig, bprimp, vecnp, retcom)
        call srcalg(dfdsp, vecnp, gp, devgii)
    end if

    !!!
    !!! 3) eq. d'equilibre : conv. meca sol sigdt - sigft + dsde:(depst-depsp-depsvp) = 0
    !!!

    if ((seuilp .ge. 0.d0) .or. (vinf(7) .gt. 0.d0)) then
        do i = 1, ndt
            depse(i) = depst(i)-depsv(i)-yf(ndt+1)*gp(i)
        end do

        dsige(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), depse(1:ndt))

        do i = 1, ndt
            r(i) = dsige(i)+sigdt(i)-sigft(i)
        end do

    else

        do i = 1, ndt
            depse(i) = depst(i)-depsv(i)
        end do

        dsige(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), depse(1:ndt))

        do i = 1, ndt
            r(i) = dsige(i)+sigdt(i)-sigft(i)
        end do

    end if

    !!! Mise a l'echelle
    do i = 1, ndt
        r(i) = r(i)/mu
    end do

    !!!
    !!! 4) condition de kt : fp=0 ou dlam=0
    !!!

    if (vinf(7) .le. 0.d0) then
        r(ndt+1) = -yf(ndt+1)
    else
        r(ndt+1) = -seuilp/mu
    end if

    !!!
    !!! 5) evolution de xip
    !!!

    coupl = materf(28, 2)
    lamgd2 = max(0.d0, yf(ndt+1)*devgii*sqrt(2.d0/3.d0))

    if ((varv .eq. 1) .and. (coupl .ge. 1.d0/2.d0)) then
        r(ndt+2) = yd(ndt+2)-yf(ndt+2)+lamgd2+dgamv
    else
        r(ndt+2) = yd(ndt+2)-yf(ndt+2)+lamgd2
    end if

    !!!
    !!! 6) evolution de xivp
    !!!

    !!! Attention si algo_inte = 'NEWTON_PERT', et si xi5 depend de T,
    !!! le calcul de drdyb dans lcjacp peut ne pas aboutir car eps2 peut etre tres faible
    dxiv = min(dgamv, xi5-yd(ndt+3))

    ! if (dxiv.lt.1.d-10) dxiv = 0.d0
    r(ndt+3) = yd(ndt+3)-yf(ndt+3)+dxiv

end subroutine
