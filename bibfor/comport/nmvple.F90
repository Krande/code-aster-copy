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

subroutine nmvple(fami, kpg, ksp, ndim, imate, &
                  compor, crit, typmod, instam, instap, &
                  deps, sigm, vim, option, sigp, &
                  vip, dsidep, iret)
! aslint: disable=
    implicit none
#include "asterfort/ggplem.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/vpalem.h"
#include "asterfort/zerofr.h"
    integer(kind=8) :: ndim, imate, iret, kpg, ksp
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: compor(*), option
    real(kind=8) :: crit(4), instam, instap
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6), vim(1), sigp(6), vip(1), dsidep(6, 6)
! ----------------------------------------------------------------------
!     REALISE LA LOI DE VISCOPLASTICITE DE LEMAITRE
!  POUR LES ELEMENTS
!     ISOPARAMETRIQUES EN PETITES DEFORMATIONS
!
!
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE
! OUT IRET    : CODE RETOUR DE LA RECHERCHE DE ZERO DE F(X)=0
!                   IRET=0 => PAS DE PROBLEME
!                   IRET=1 => ECHEC
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX YY ZZ XY XZ YZ
!
! ----------------------------------------------------------------------
!
!     COMMON POUR LES PARAMETRES DES LOIS VISCOPLASTIQUES
    common/nmpavp/dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
    real(kind=8) :: dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
!     COMMON POUR LES PARAMETRES DE LA LOI DE LEMAITRE (NON IRRADIEE)
    common/nmpale/unsurk, unsurm, valden
    real(kind=8) :: unsurk, unsurm, valden
!
    real(kind=8) :: depsth(6), valres(5), epsthe
    real(kind=8) :: depsdv(6), sigdv(6), sigel(6), epsmo, sigmo, e, nu
    real(kind=8) :: troisk, valpar(2), rac2, t1, t2
    real(kind=8) :: em, num, troikm, deumum, sigmp(6)
    real(kind=8) :: deltkl, deltp2
    real(kind=8) :: degran(6)
    integer(kind=8) :: k, l, iret1, ibid
    integer(kind=8) :: ndimsi, iret2
    real(kind=8) :: a0, xap, x, tm, tp
    real(kind=8) :: fg, fdgdst, fdgdev, defam(6), defap(6)
    real(kind=8) :: coef1, coef2, deltev
    integer(kind=8) :: icodre(5)
    character(len=8) :: nompar(2)
    character(len=16) :: nomres(5)
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    character(len=6), parameter :: epsa(6) = (/'EPSAXX', 'EPSAYY', 'EPSAZZ', &
                                               'EPSAXY', 'EPSAXZ', 'EPSAYZ'/)
! DEB ------------------------------------------------------------------
!
    call verift(fami, kpg, ksp, 'T', imate, &
                epsth_=epsthe)
!
    iret = 0
    theta = crit(4)
    t1 = abs(theta-0.5d0)
    t2 = abs(theta-1.d0)
    prec = 0.01d0
    if ((t1 .gt. prec) .and. (t2 .gt. prec)) then
        call utmess('F', 'ALGORITH6_55')
    end if
!
    if (typmod(1) .eq. 'C_PLAN') then
        call utmess('F', 'ALGORITH6_92')
    end if
!
    degran(1:6) = 0.d0
    rac2 = sqrt(2.d0)
    call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                ksp, tm, iret1)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, tp, iret2)
    if ((iret1+iret2) .eq. 0) then
        tschem = tm*(1.d0-theta)+tp*theta
    else
        tschem = 0.d0
    end if
    dpc = vim(1)
    deltat = instap-instam
!
    dsidep(:, :) = 0.d0
!
    if (ndim .eq. 2) then
        ndimsi = 4
    else
        ndimsi = 6
    end if
! VARIABLE DE COMMANDE ANELASTIQUE
!
    do k = 1, ndimsi
        call rcvarc(' ', epsa(k), '-', fami, kpg, &
                    ksp, defam(k), iret2)
        if (iret2 .eq. 1) defam(k) = 0.d0
!
        call rcvarc(' ', epsa(k), '+', fami, kpg, &
                    ksp, defap(k), iret2)
        if (iret2 .eq. 1) defap(k) = 0.d0
    end do
!
!
! MISE AU FORMAT DES TERMES NON DIAGONAUX
!
    do k = 4, ndimsi
        defam(k) = defam(k)*rac2
        defap(k) = defap(k)*rac2
    end do
!
    nompar(1) = 'INST'
    valpar(1) = instam
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 1, nompar, [valpar], &
                2, nomres, valres, icodre, 2)
    em = valres(1)
    num = valres(2)
    deumum = em/(1.d0+num)
    troikm = em/(1.d0-2.d0*num)
    valpar(1) = instap
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 1, nompar, [valpar], &
                2, nomres, valres, icodre, 2)
    e = valres(1)
    nu = valres(2)
    deuxmu = e/(1.d0+nu)
    troisk = e/(1.d0-2.d0*nu)
!
!
    nomres(1) = 'N'
    nomres(2) = 'UN_SUR_K'
    nomres(3) = 'UN_SUR_M'
    nompar(1) = 'TEMP'
    valpar(1) = tschem
    call rcvalb(fami, 1, 1, '+', imate, &
                ' ', 'LEMAITRE', 1, nompar, [valpar], &
                3, nomres, valres, icodre, 2)
    valden = valres(1)
    unsurk = valres(2)
    unsurm = valres(3)
!
    epsmo = 0.d0
    do k = 1, 3
        depsth(k) = deps(k)-epsthe-(defap(k)-defam(k))
        depsth(k) = depsth(k)-degran(k)
        depsth(k) = depsth(k)*theta
        if ((k .eq. 1) .or. (ndimsi .eq. 6)) then
            depsth(k+3) = deps(k+3)-(defap(k+3)-defam(k+3))
            depsth(k+3) = depsth(k+3)-degran(k+3)
            depsth(k+3) = depsth(k+3)*theta
        end if
        epsmo = epsmo+depsth(k)
    end do
!
    epsmo = epsmo/3.d0
    do k = 1, ndimsi
        depsdv(k) = depsth(k)-epsmo*kron(k)
    end do
!
    sigmo = 0.d0
    do k = 1, 3
        sigmo = sigmo+sigm(k)
    end do
    sigmo = sigmo/3.d0
!
    do k = 1, ndimsi
        sigmp(k) = (theta*deuxmu+(1.d0-theta)*deumum)/deumum*(sigm(k)- &
                                                              sigmo*kron(k))+ &
                   (theta*troisk+(1.d0-theta)*troikm)/troikm* &
                   sigmo*kron(k)
    end do
!
    sigmo = 0.d0
    do k = 1, 3
        sigmo = sigmo+sigmp(k)
    end do
    sigmo = sigmo/3.d0
    sieleq = 0.d0
    do k = 1, ndimsi
        sigdv(k) = sigmp(k)-sigmo*kron(k)
        sigel(k) = sigdv(k)+deuxmu*depsdv(k)
        sieleq = sieleq+sigel(k)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
!
!----RESOLUTION DE L'EQUATION SCALAIRE----
!
    prec = crit(3)
    niter = nint(crit(1))
!
    a0 = -sieleq
!
    xap = sieleq
    xap = xap-sieleq*1.d-12
    if (abs(a0) .le. prec) then
        x = 0.d0
    else
        call zerofr(0, 'DEKKER2', vpalem, 0.d0, xap, &
                    prec, int(niter), x, iret, ibid)
        if (iret .eq. 1) then
            if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
                goto 999
            else
                x = 0.d0
                iret = 0
            end if
        end if
    end if
    if (x .ne. 0.d0) call ggplem(x, dpc+(sieleq-x)/(1.5d0*deuxmu), valden, unsurk, unsurm, &
                                 theta, deuxmu, fg, fdgdst, fdgdev)
!
!-----------------------------------------
    if (x .ne. 0.d0) then
        coef1 = 1.d0/(1.d0+1.5d0*deuxmu*deltat*fg/x)
    else
!       COEF1 = 1.d0/(1.d0+1.5d0*DEUXMU*DELTAT*FDGDST)
        coef1 = 1.d0
    end if
!
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        deltp2 = 0.d0
        do k = 1, ndimsi
            sigdv(k) = sigel(k)*coef1
            sigp(k) = sigdv(k)+(sigmo+troisk*epsmo)*kron(k)
            sigp(k) = (sigp(k)-sigm(k))/theta+sigm(k)
            deltev = (sigel(k)-sigdv(k))/(deuxmu*theta)
            deltp2 = deltp2+deltev**2
        end do
        vip(1) = vim(1)+sqrt(2.d0*deltp2/3.d0)
    end if
!
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:14) .eq. 'RIGI_MECA_TANG') then
        if (x .ne. 0.d0) then
            coef2 = sieleq*(1.d0-deltat*fdgdev)
            coef2 = coef2/(1.d0+1.5d0*deuxmu*deltat*fdgdst)
            coef2 = coef2-x
            coef2 = coef2*1.5d0/(sieleq**3)
        else
            coef2 = 0.d0
        end if
        do k = 1, ndimsi
            do l = 1, ndimsi
                deltkl = 0.d0
                if (k .eq. l) deltkl = 1.d0
                dsidep(k, l) = coef1*(deltkl-kron(k)*kron(l)/3.d0)
                dsidep(k, l) = deuxmu*(dsidep(k, l)+coef2*sigel(k)*sigel(l))
                dsidep(k, l) = dsidep(k, l)+troisk*kron(k)*kron(l)/3.d0
            end do
        end do
    end if
!
999 continue
!
! FIN ------------------------------------------------------------------
end subroutine
