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
subroutine cjsmid(mod, crit, mater, nvi, epsd, &
                  deps, sigd, sigf, vind, vinf, &
                  noconv, aredec, stopnc, niter, epscon)
    implicit none
!     INTEGRATION PLASTIQUE (MECANISMES ISOTROPE ET DEVIATOIRE) DE CJS
!     IN   MOD      :  MODELISATION
!          CRIT     :  CRITERES DE CONVERGENCE
!          MATER    :  COEFFICIENTS MATERIAU A T+DT
!          NVI      :  NB DE VARIABLES INTERNES
!          EPSD     :  DEFORMATIONS A T
!          DEPS     :  INCREMENTS DE DEFORMATION
!          SIGD     :  CONTRAINTE  A T
!          VIND     :  VARIABLES INTERNES  A T
!          AREDEC   :  ARRET DES DECOUPAGES
!          STOPNC   :  ARRET EN CAS DE NON CONVERGENCE
!     VAR  SIGF     :  CONTRAINTE  A T+DT
!          VINF     :  VARIABLES INTERNES  A T+DT
!          NOCONV   :  PAS DE CONVERGENCE
!          NITER    :  NOMBRE D ITERATIONS A CONVERGENCE
!          EPSCON   :  VALEUR ERR FINALE
!       ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/cjsiid.h"
#include "asterfort/cjsjid.h"
#include "asterfort/cjsncn.h"
#include "asterfort/cjsncv.h"
#include "asterfort/cjsnor.h"
#include "asterfort/iunifi.h"
#include "asterfort/mgauss.h"
    integer(kind=8) :: ndt, ndi, nvi, nr, nmod, niter, iret
    integer(kind=8) :: nitimp
    parameter(nmod=16)
    parameter(nitimp=200)
!
    integer(kind=8) :: iter
    aster_logical :: noconv, aredec, stopnc
!
    real(kind=8) :: epsd(6), deps(6)
    real(kind=8) :: sigd(6), sigf(6), gd(6)
    real(kind=8) :: vind(*), vinf(*), epscon
    real(kind=8) :: crit(*), mater(14, 2)
    real(kind=8) :: r(nmod), drdy(nmod, nmod)
    real(kind=8) :: ddy(nmod), dy(nmod), yd(nmod), yf(nmod)
    real(kind=8) :: err, err1, err2, signe
    real(kind=8) :: det, pa, qinit
    integer(kind=8) :: umess
!
    integer(kind=8) :: essai, essmax
    parameter(essmax=10)
!
!    SI ABS(COS_NORMALES) < TOLROT RELAX = RELAX*DECREL
!
    real(kind=8) :: tolrot, decrel
    parameter(tolrot=0.8d0)
    parameter(decrel=0.5d0)
!
    real(kind=8) :: relax(essmax+1), rotagd(essmax+1), xf(6), nor1(7), nor2(7)
    real(kind=8) :: erimp(nitimp, 4)
!
    aster_logical :: devnu1, devnu2, tra1, tra2
    integer(kind=8) :: i, j
!
    character(len=8) :: mod
!
    common/tdim/ndt, ndi
!
!
!     ------------------------------------------------------------------
!
!
! -> DIMENSION DU PROBLEME:
!    NR = NDT(SIG)+ 1(QISO) + 1(R)+ NDT(X) + 1(DLAMBI) + 1(DLAMBD)
!
!
    umess = iunifi('MESSAGE')
    noconv = .false.
    pa = mater(12, 2)
    qinit = mater(13, 2)
!
!
!
!
    nr = 2*ndt+4
!
! -> MISE A ZERO DES DATAS
!
    do i = 1, nr
        ddy(i) = 0.d0
        dy(i) = 0.d0
        yd(i) = 0.d0
        yf(i) = 0.d0
    end do
!
!
    do i = 1, ndt
        gd(i) = 0.d0
    end do
!
!
! -> INITIALISATION DE YD PAR LES CHAMPS (SIGD, VIND, ZERO)
!
    yd(1:ndt) = sigd(1:ndt)
    yd(ndt+1) = vind(1)
    yd(ndt+2) = vind(2)
!
    do i = 1, ndt
        yd(ndt+2+i) = vind(i+2)
    end do
    yd(2*ndt+3) = 0.d0
    yd(2*ndt+4) = 0.d0
!
!
! -> INITIALISATION : DY : CALCUL DE LA SOLUTION D ESSAI INITIALE EN DY
!    (SOLUTION EXPLICITE)
!
    call cjsiid(mod, mater, epsd, deps, yd, &
                gd, dy)
!
!
!
!
! -> INCREMENTATION DE YF = YD + DY
!
    yf(1:nr) = yd(1:nr)+dy(1:nr)
!
!
!---------------------------------------
! -> BOUCLE SUR LES ITERATIONS DE NEWTON
!---------------------------------------
!
    iter = 0
100 continue
!
    iter = iter+1
!
!
!--->   EN CAS D'ENTREE EN TRACTION, ON A
!          1.  LES CONTRAINTES SONT RAMENEES SUR L'AXE HYDROSTATIQUE
!              A DES VALEURS FAIBLES ( EGALES A PA/100.0 SOIT -1 KPA )
!          2.  LES VARIABLES INTERNES N'EVOLUENT PAS
!          3.  ON SORT DIRECTEMENT DE CJSMID, CAR ON SUPPOSE ALORS
!              ETRE REVENU DANS LE DOMAINE ELASTIQUE
!
!       SINON IL APPARAIT ENSUITE UNE ERREUR DANS LA ROUTINE MGAUSS
!
!
    if ((yf(1)+yf(2)+yf(3)) .ge. -qinit) then
        do i = 1, ndi
            sigf(i) = -qinit/3.d0+pa/100.0d0
        end do
        do i = ndi+1, ndt
            sigf(i) = 0.d0
        end do
        vinf(1:nvi-1) = vind(1:nvi-1)
        vinf(nvi) = 0.d0
        goto 999
    end if
!
!
! -> CALCUL DU SECOND MEMBRE A T+DT :  -R(DY)
!    CALCUL DE SIGNE(S:DEPSDP)
! ET CALCUL DU JACOBIEN DU SYSTEME A T+DT :  DRDY(DY)
!
    do i = 1, nr
        r(i) = 0.d0
        do j = 1, nr
            drdy(i, j) = 0.d0
        end do
    end do
!
    call cjsjid(mod, mater, epsd, deps, yd, &
                yf, gd, r, signe, drdy)
!
!
! -> RESOLUTION DU SYSTEME LINEAIRE : DRDY(DY).DDY = -R(DY)
!
    ddy(1:nr) = r(1:nr)
    call mgauss('NFVP', drdy, ddy, nmod, nr, &
                1, det, iret)
!
    relax(1) = 1.d0
!
!   ESTIMATION NORMALE AU POINT YD + DY
!
    do i = 1, ndt
        sigf(i) = yd(i)+dy(i)
        xf(i) = yd(ndt+2+i)+dy(ndt+2+i)
    end do
    call cjsnor(mater, sigf, xf, nor1, devnu1, &
                tra1)
!
    essai = 0
!
40  continue
    essai = essai+1
    if ((.not. devnu1) .and. (.not. tra1)) then
        if (essai .gt. essmax) then
            if (aredec .and. stopnc) then
                call cjsncn('CJSMID  ', essmax, ndt, nvi, umess, &
                            relax, rotagd, epsd, deps, sigd, &
                            vind)
            else
                noconv = .true.
                goto 200
            end if
        end if
!
!   ESTIMATION NORMALE AU POINT YD + RELAX*DY
!
        do i = 1, ndt
            sigf(i) = yd(i)+dy(i)+relax(essai)*ddy(i)
            xf(i) = yd(ndt+2+i)+dy(ndt+2+i)+relax(essai)*ddy(ndt+2+i)
        end do
        call cjsnor(mater, sigf, xf, nor2, devnu2, &
                    tra2)
!
        rotagd(essai) = 0.d0
        do i = 1, ndt
            rotagd(essai) = rotagd(essai)+nor1(i)*nor2(i)
        end do
        rotagd(essai) = rotagd(essai)/(nor1(ndt+1)*nor2(ndt+1))
!
        if (abs(rotagd(essai)) .lt. tolrot .and. (.not. devnu2) .and. (.not. tra2)) then
            relax(essai+1) = relax(essai)*decrel
            goto 40
        end if
    end if
!
!  DANS LES CAS OU DEVNU1 OU DEVNU2 (ON A DETCTE DES DEVIATEURS NULS)
!  OU TRA1 OU TRA2 ( ON A DETECTE DES TRACTIONS) ON ABANDONNE
!  LA RELAXATION SUR LES NORMALES
!
!
    do i = 1, nr
        dy(i) = dy(i)+relax(essai)*ddy(i)
        yf(i) = yd(i)+dy(i)
    end do
!
!
!
!
! -> VERIFICATION DE LA CONVERGENCE : ERREUR = !!DDY!!/!!DY!! < TOLER
!
!
    err1 = norm2(ddy(1:nr))
    err2 = norm2(dy(1:nr))
    if (err2 .eq. 0.d0) then
        err = err1
    else
        err = err1/err2
    end if
    if (iter .le. nitimp) then
        erimp(iter, 1) = err1
        erimp(iter, 2) = err2
        erimp(iter, 3) = err
        erimp(iter, 4) = relax(essai)
    end if
!
    if (iter .le. int(abs(crit(1)))) then
!
!          --   CONVERVENCE   --
        if (err .lt. crit(3)) then
            goto 200
!
!          --  NON CONVERVENCE : ITERATION SUIVANTE  --
        else
            goto 100
        end if
!
    else
!
!          --  NON CONVERVENCE : ITERATION MAXI ATTEINTE  --
!
        if (aredec .and. stopnc) then
            call cjsncv('CJSMID', nitimp, iter, ndt, nvi, &
                        umess, erimp, epsd, deps, sigd, &
                        vind)
        else
            noconv = .true.
        end if
    end if
!
200 continue
    niter = iter
    epscon = err
!
!
! -> MISE A JOUR DES CONTRAINTES ET VARIABLES INTERNES
!
    sigf(1:ndt) = yf(1:ndt)
    vinf(1) = yf(ndt+1)
    vinf(2) = yf(ndt+2)
    do i = 1, ndt
        vinf(2+i) = yf(ndt+2+i)
    end do
    vinf(nvi-1) = signe
!
!
999 continue
!
end subroutine
