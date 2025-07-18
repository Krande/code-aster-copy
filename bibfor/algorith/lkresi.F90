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

subroutine lkresi(typmod, nmat, materf, timed, timef, &
                  nvi, vind, vinf, yd, yf, &
                  deps, nr, r)
! person_in_charge: alexandre.foucault at edf.fr
    implicit none
!       ----------------------------------------------------------------
!       CALCUL DES TERMES DU SYSTEME NL A RESOUDRE = -R(DY) POUR LETK
!       IN  TYPMOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           TIMED  :  INSTANT  T
!           TIMEF  :  INSTANT  T+DT
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           DEPS   :  INCREMENT DE DEFORMATION
!           VIND   :  VARIABLES INTERNES A T
!           VINF   :  VARIABLES INTERNES A T+DT
!           YD     :  VARIABLES A T    = ( SIGD 0    XIPD XIVPD (EPSD3))
!           YF     :  VARIABLES A T+DT = ( SIGF DLAM XIPF XIVPF (EPS3F))
!           DEPS   :  INCREMENT DE DEFORMATIONS
!           DY     :  SOLUTION         = ( DSIG DLAM DXIP DXIVP (DEPS3))
!           NR     :  DIMENSION DU VECTEUR INCONNUES
!       OUT R      :  SYSTEME NL A T+DT
!       ----------------------------------------------------------------
#include "asterfort/lcdevi.h"
#include "asterfort/lkbpri.h"
#include "asterfort/lkcalg.h"
#include "asterfort/lkcaln.h"
#include "asterfort/lkcrip.h"
#include "asterfort/lkcriv.h"
#include "asterfort/lkdfds.h"
#include "asterfort/lkdgde.h"
#include "asterfort/lkdhds.h"
#include "asterfort/lkds2h.h"
#include "asterfort/lkelas.h"
#include "asterfort/lkvacp.h"
#include "asterfort/lkvarp.h"
    integer(kind=8) :: nmat, nr, nvi, ndi, ndt
    real(kind=8) :: deps(6), vind(*), vinf(*)
    real(kind=8) :: r(nr), yd(nr), yf(nr), materf(nmat, 2)
    real(kind=8) :: timed, timef
    character(len=8) :: typmod
!
    integer(kind=8) :: i, retcom, val, varv
    real(kind=8) :: zero, vint(7), devsig(6), i1, ucrip, seuilp
    real(kind=8) :: dt, seuilv, depsv(6), dgamv
    real(kind=8) :: dxiv, xivmax, xippic, seuivm, ucriv
    real(kind=8) :: dsdenl(6, 6), kk, mu, dhds(6), ds2hds(6)
    real(kind=8) :: paraep(3), varpl(4), dfdsp(6), bprimp
    real(kind=8) :: vecnp(6), gp(6), devgii, deux, trois, depse(6)
    real(kind=8) :: dsige(6), sigdt(6), sigft(6), depst(6), lamgd2
    parameter(zero=0.d0)
    parameter(deux=2.d0)
    parameter(trois=3.d0)
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
!       --------------------------------------------------------------
!
! --------------------------------------------------------------------
! --- PASSAGE EN CONVENTION MECANIQUE DES SOLS
! --------------------------------------------------------------------
    sigft = 0.d0
    sigdt = 0.d0
    depst = 0.d0
    do i = 1, ndt
        sigft(i) = -yf(i)
        sigdt(i) = -yd(i)
        depst(i) = -deps(i)
    end do
!
! ----------------------------------------------------------------------
! --- VARIABLES LOCALES TEMPORAIRES
! ----------------------------------------------------------------------
    varv = 0
    val = 0
    devgii = zero
!
! --- VECTEUR VARIABLES INTERNES TEMPORAIRES
    vint(1:nvi) = vind(1:nvi)
!
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
!
! --- INCREMENT DE TEMPS
    dt = timef-timed
!
! --- CONSTRUCTION TENSEUR DEVIATOIRE DES CONTRAINTES ET 1ER INVARIANT
    call lcdevi(sigft, devsig)
    i1 = sigft(1)+sigft(2)+sigft(3)
!
! --- DONNEES MATERIAU : VALEUR MAX DE XIV; XI_PIC
    xivmax = materf(20, 2)
    xippic = materf(18, 2)
!
! --- CONSTRUCTION TENSEUR ELASTIQUE NON LINEAIRE DSDENL
    call lkelas(ndi, ndt, nmat, materf, depst, &
                sigft, dsdenl, kk, mu)
!
! ----------------------------------------------------------------------
! --- I) - BUT : CALCUL DE LA DEFORMATION VISQUEUSE -DEPSV- ET DU
! ---      PARAMETRE D ECROUISSAGE VISQUEUX -DGAMV-
! ----------------------------------------------------------------------
!
! --- I-1) INDICATEUR SUR ANGLE DE DILATANCE VISQUEUX PSI -> VAL = 0
    val = 0
!
! --- I-2) VARIABLE D'ECROUISSAGE VISQUEUSE VINTR = YF(NDT+3)
! --- I-3) CALCUL SEUIL VISQUEUX PAR RAPPORT A YF(1:6)=SIGF -> SEUILV
! --- I-3-1)  XIT   = YF(NDT+3)
    call lkcriv(vint(3), i1, devsig, vint, nmat, &
                materf, ucriv, seuilv)
!
    if (seuilv .ge. zero) then
        call lkdgde(val, vint(3), dt, seuilv, ucriv, &
                    i1, devsig, vint, nmat, materf, &
                    depsv, dgamv, retcom)
    else
        dgamv = zero
        do i = 1, ndt
            depsv(i) = zero
        end do
    end if
! ----------------------------------------------------------------------
! --- II) - BUT : CALCUL DE LA DEFORMATION PLASTIQUE -DEPSP- ET DU
! ---       PARAMETRE D ECROUISSAGE PLASTIQUE -DGAMP-
! ----------------------------------------------------------------------
! --- II-2-B-2) INDICATEUR CONTRACTANCE OU DILATANCE -> VARV = 0 OU 1
! --- II-2-B-2)-1) CALCUL POSITION YF PAR RAPPORT SEUIL VISQUEUX MAX
    call lkcriv(xivmax, i1, devsig, vint, nmat, &
                materf, ucriv, seuivm)
!
! --- II-2-B-2)-2) TEST SUR SEUIL >0 OU <0 POUR DEFINIR VARV
    if (seuivm .le. zero) then
        varv = 0
    else
        varv = 1
    end if
!
! --- II-1) CALCUL FONCTION SEUIL PLASTIQUE EN YF
    seuilp = zero
    call lkcrip(i1, devsig, vint, nmat, materf, &
                ucrip, seuilp)
!
! --- II-2)SI SEUILP >= 0 ALORS PLASTICITE A PRENDRE EN COMPTE
    if ((seuilp .ge. zero) .or. (vinf(7) .gt. zero)) then
! --- II-2-B-1) INDICATEUR ANGLE DE DILATANCE PLASTIQUE PSI -> 0 OU 1
        if (yf(ndt+2) .le. xippic) then
            val = 0
        else
            val = 1
        end if
!
! --- II-2-B-3) CALCUL DE DF/DSIG
        call lkdhds(nmat, materf, i1, devsig, dhds, &
                    retcom)
        call lkds2h(nmat, materf, i1, devsig, dhds, &
                    ds2hds, retcom)
        call lkvarp(vint, nmat, materf, paraep)
        call lkvacp(nmat, materf, paraep, varpl)
        call lkdfds(nmat, materf, devsig, paraep, varpl, &
                    ds2hds, ucrip, dfdsp)
!
! --- II-2-B-4) CALCUL DE G
        bprimp = lkbpri(val, vint, nmat, materf, paraep, i1, devsig)
        call lkcaln(devsig, bprimp, vecnp, retcom)
        call lkcalg(dfdsp, vecnp, gp, devgii)
    end if
!
! ----------------------------------------------------------------------
! --- III) EQUATION D'EQUILIBRE : (CONVENTION MECANIQUE DES SOLS)
! ---      SIGDT - SIGFT + DSDE:(DEPST-DEPSP-DEPSVP) = 0
! ----------------------------------------------------------------------
    if ((seuilp .ge. zero) .or. (vinf(7) .gt. zero)) then
        do i = 1, ndt
            depse(i) = depst(i)-depsv(i)-yf(ndt+1)*gp(i)
        end do
!
        dsige(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), depse(1:ndt))
!
        do i = 1, ndt
            r(i) = dsige(i)+sigdt(i)-sigft(i)
        end do
    else
        do i = 1, ndt
            depse(i) = depst(i)-depsv(i)
        end do
!
        dsige(1:ndt) = matmul(dsdenl(1:ndt, 1:ndt), depse(1:ndt))
!
        do i = 1, ndt
            r(i) = dsige(i)+sigdt(i)-sigft(i)
        end do
    end if
! === =================================================================
! --- MISE A L'ECHELLE DE DEFORMATIONS -> R(I)/MODULE_CISAILLEMENT
! === =================================================================
    do i = 1, ndt
        r(i) = r(i)/mu
    end do
!
! ----------------------------------------------------------------------
! --- IV) CONDITION DE KUHN-TUCKER : -FP = 0 OU -DLAM = 0
! ----------------------------------------------------------------------
! --- APPLICATION DE LA CONDITION DE KHUN-TUCKER SUR R(NDT+1)
    if (vinf(7) .eq. zero) then
        r(ndt+1) = -yf(ndt+1)
    else
        r(ndt+1) = -seuilp/mu
    end if
! ----------------------------------------------------------------------
! --- V) EVOLUTION DE XIP :
! ---    XIPD - XIPF + DLAM*G_II*SQRT(2/3)(+ DGAMVP) = 0
! ----------------------------------------------------------------------
    lamgd2 = max(zero, yf(ndt+1)*devgii*sqrt(deux/trois))
!
    if (varv .eq. 0) then
        r(ndt+2) = yd(ndt+2)-yf(ndt+2)+lamgd2
    else
        r(ndt+2) = yd(ndt+2)-yf(ndt+2)+lamgd2+dgamv
    end if
! ----------------------------------------------------------------------
! --- VI) EVOLUTION DE XIVP :
! ---     XIVPD - XIVPF + MIN(DGAM_VP,XIV_MAX-XIVPD) = 0
! ----------------------------------------------------------------------
! --- TEST POUR DEFINIR MIN(DGAMV,XIV_MAX-XIV)
    dxiv = min(dgamv, xivmax-yd(ndt+3))
!
    r(ndt+3) = yd(ndt+3)-yf(ndt+3)+dxiv
!
end subroutine
