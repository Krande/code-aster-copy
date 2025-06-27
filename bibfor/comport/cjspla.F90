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
subroutine cjspla(mod, crit, mater, seuili, seuild, &
                  nvi, epsd, deps, sigd, vind, &
                  sigf, vinf, mecani, nivcjs, niter, &
                  ndec, epscon, iret, trac)
    implicit none
!       INTEGRATION PLASTIQUE DE LA LOI CJS
!       IN  MOD    :  MODELISATION
!           CRIT   :  CRITERES DE CONVERGENCE
!           MATER  :  COEFFICIENTS MATERIAU A T+DT
!           SEUILI :  FONCTION DE CHARGE ISO. CALCULEE AVEC PREDICT ELAS
!           SEUILD :  FONCTION DE CHARGE DEV. CALCULEE AVEC PREDICT ELAS
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           EPSD   :  DEFORMATIONS A T
!           DEPS   :  INCREMENT DE DEFORMATION
!           SIGD   :  CONTRAINTE  A T
!           VIND   :  VARIABLES INTERNES  A T
!       VAR SIGF   :  CONTRAINTE A T+DT  (IN -> ELAS, OUT -> PLASTI )
!       OUT VINF   :  VARIABLES INTERNES A T+DT
!           MECANI :  MECANISME(S) ACTIVE(S)
!           NITER  :  NOMBRE D ITERATIONS POUR PLASTICITE
!                          (CUMUL DECOUPAGE)
!           NDEC   :  NOMBRE DE DECOUPAGE
!           EPSCON :  EPSILON A CONVERGENCE
!           IRET   :  CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!       ----------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cjsmde.h"
#include "asterfort/cjsmid.h"
#include "asterfort/cjsmis.h"
#include "asterfort/cjssmd.h"
#include "asterfort/cjssmi.h"
    integer(kind=8) :: ndt, ndi, nvi, niter, ndec, iret
    integer(kind=8) :: nvimax
    parameter(nvimax=16)
    real(kind=8) :: epsd(6), deps(6)
    real(kind=8) :: sigd(6), sigf(6), predic(6)
    real(kind=8) :: sigd0(6), deps0(6), predi0(6)
    real(kind=8) :: vind(*), vinf(*), vind0(nvimax), epscon
    real(kind=8) :: mater(14, 2), crit(*)
    real(kind=8) :: i1f
    real(kind=8) :: seuili, seuild, pa, pref, qinit
    real(kind=8) :: zero
    aster_logical :: chgmec, noconv, aredec, stopnc, trac
    character(len=6) :: mecani
    character(len=4) :: nivcjs
!
    character(len=8) :: mod
    parameter(zero=0.d0)
    integer(kind=8) :: idec
    integer(kind=8) :: i, niter0
!
    common/tdim/ndt, ndi
!
!
    ASSERT(nvi .le. nvimax)
!
    pa = mater(12, 2)
    qinit = mater(13, 2)
    trac = .false.
!
!  SAUVEGARDE DES GRANDEURS D ENTREE INITIALES
    predi0(1:ndt) = sigf(1:ndt)
    sigd0(1:ndt) = sigd(1:ndt)
    deps0(1:ndt) = deps(1:ndt)
    vind0(1:nvi) = vind(1:nvi)
!
!
!  ARRET OU NON EN NON CONVERGENCE INTERNE
!  -------------------------------------------
!
!
!
    if (int(crit(1)) .lt. 0) then
        stopnc = .true.
    else
        stopnc = .false.
    end if
!
!
!  INITIALISATION DES VARIABLES DE REDECOUPAGE
!  -------------------------------------------
!
! INT(CRIT(5)) = 0  1 OU -1 -> PAS DE REDECOUPAGE DU PAS DE TEMPS
!
!
    if ((int(crit(5)) .eq. 0) .or. (int(crit(5)) .eq. -1) .or. (int(crit(5)) .eq. 1)) then
        ndec = 1
        aredec = .true.
        noconv = .false.
!
!
! INT(CRIT(5)) < -1 -> REDECOUPAGE DU PAS DE TEMPS SI NON CONVERGENCE
!
    else if (int(crit(5)) .lt. -1) then
        ndec = 1
        aredec = .false.
        noconv = .false.
!
!
! INT(CRIT(5)) > 1 -> REDECOUPAGE IMPOSE DU PAS DE TEMPS
!
    else if (int(crit(5)) .gt. 1) then
        ndec = int(crit(5))
        aredec = .true.
        noconv = .false.
    end if
!
!
!
!  POINT DE RETOUR EN CAS DE DECOUPAGE
!  APRES UNE NON CONVERGENCE, POUR  INT(CRIT(5)) < -1
!
500 continue
    if (noconv) then
        ndec = -int(crit(5))
        aredec = .true.
    end if
!
!
!   RESTAURATION DE SIGD VIND DEPS ET PREDIC ELAS SIGF
!   EN TENANT COMPTE DU DECOUPAGE EVENTUEL
!
    sigd(1:ndt) = sigd0(1:ndt)
    vind(1:nvi) = vind0(1:nvi)
!
    do i = 1, ndt
        deps(i) = deps0(i)/ndec
        sigf(i) = sigd0(i)+(predi0(i)-sigd(i))/ndec
    end do
!
!
!  BOUCLE SUR LES DECOUPAGES
!  -------------------------
    do idec = 1, ndec
!
!
! SAUVEGARDE PREDIC ELASTIQUE POUR EVENTUEL CHANGEMENT
! DE MECANISME
!
        predic(1:ndt) = sigf(1:ndt)
!
        i1f = zero
        do i = 1, ndi
            i1f = i1f+sigf(i)
        end do
!
        if ((i1f+qinit) .eq. 0.d0) then
            i1f = -qinit+1.d-12*pa
            pref = abs(pa)
        else
            pref = abs(i1f+qinit)
        end if
!
        call cjssmi(mater, sigf, vind, seuili)
        call cjssmd(mater, sigf, vind, seuild)
        seuili = seuili/pref
        seuild = seuild/pref
        chgmec = .false.
        if (seuili .gt. zero .and. seuild .le. zero) then
            mecani = 'ISOTRO'
        end if
        if (seuili .le. zero .and. seuild .gt. zero) then
            mecani = 'DEVIAT'
        end if
        if (seuili .gt. zero .and. seuild .gt. zero) then
            mecani = 'ISODEV'
        end if
!
!
        do i = 1, nvi-1
            vinf(i) = vind(i)
        end do
!
!
100     continue
!
!--->   RESOLUTION EN FONCTION DES MECANISMES ACTIVES
!
!       MECANISME ISOTROPE SEUL
!       -----------------------
!
        if (mecani .eq. 'ISOTRO') then
            call cjsmis(mod, crit, mater, nvi, epsd, &
                        deps, sigd, sigf, vind, vinf, &
                        noconv, aredec, stopnc, niter0, epscon)
            niter = niter+niter0
            if (noconv .and. (.not. aredec)) goto 500
            if (noconv) then
                iret = 1
                goto 999
            end if
        end if
!
!
!       MECANISME DEVIATOIRE SEUL
!       -------------------------
!
        if (mecani .eq. 'DEVIAT') then
            call cjsmde(mod, crit, mater, nvi, epsd, &
                        deps, sigd, sigf, vind, vinf, &
                        noconv, aredec, stopnc, niter0, epscon, &
                        trac)
            niter = niter+niter0
            if (trac) goto 999
!
            if (noconv .and. (.not. aredec)) goto 500
            if (noconv) then
                iret = 1
                goto 999
            end if
        end if
!
!       MECANISMES ISOTROPE ET DEVIATOIRE
!       ---------------------------------
!
        if (mecani .eq. 'ISODEV') then
            call cjsmid(mod, crit, mater, nvi, epsd, &
                        deps, sigd, sigf, vind, vinf, &
                        noconv, aredec, stopnc, niter0, epscon)
            niter = niter+niter0
            if (noconv .and. (.not. aredec)) goto 500
            if (noconv) then
                iret = 1
                goto 999
            end if
        end if
!
!
!
!--->   CALCUL DES FONCTIONS DE CHARGES SUR ETAT FINAL
!
        call cjssmi(mater, sigf, vinf, seuili)
        call cjssmd(mater, sigf, vinf, seuild)
!
        i1f = zero
        do i = 1, ndi
            i1f = i1f+sigf(i)
        end do
        if ((i1f+qinit) .eq. 0.d0) then
            i1f = -qinit+1.d-12*pa
            pref = abs(pa)
        else
            pref = abs(i1f+qinit)
        end if
!--->   VERIFICATION DES MECANISMES ACTIVES
!
!
        if ((mecani .eq. 'ISOTRO') .and. (seuild .gt. zero)) then
            mecani = 'ISODEV'
            chgmec = .true.
        end if
        if ((mecani .eq. 'DEVIAT') .and. (seuili .gt. zero)) then
            mecani = 'ISODEV'
            chgmec = .true.
        end if
!
!
! - SI ON ACTIVE EN FAIT LES DEUX MECANISMES AU LIEU D'UN SEUL : RETOUR
!   ET SI ON AVAIT CONVERGE
!
        if (chgmec .and. (.not. noconv)) then
            chgmec = .false.
            sigf(1:ndt) = predic(1:ndt)
            goto 100
        else
            if (idec .lt. ndec) then
                sigd(1:ndt) = sigf(1:ndt)
                do i = 1, nvi-1
                    vind(i) = vinf(i)
                end do
                do i = 1, ndt
                    sigf(i) = sigd(i)+(predi0(i)-sigd(i))/ndec
                end do
            end if
        end if
    end do
!
!
!--->   CALCUL DE LA VARIABLE INTERNE CORRESPONDANT AU MECANISME PLASTIC
!
!
    if (mecani .eq. 'ISOTRO') vinf(nvi) = 1.d0
    if (mecani .eq. 'DEVIAT') vinf(nvi) = 2.d0
    if (mecani .eq. 'ISODEV') vinf(nvi) = 3.d0
!
999 continue
end subroutine
