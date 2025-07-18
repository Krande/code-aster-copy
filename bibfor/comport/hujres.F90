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

subroutine hujres(fami, kpg, ksp, mod, crit, &
                  mater, imat, nvi, deps, sigd, &
                  vind, sigf, vinf, iret, etatf)
    implicit none
!   INTEGRATION PLASTIQUE DE LA LOI DE HUJEUX
!   IN  MOD    :  MODELISATION
!       CRIT   :  CRITERES DE CONVERGENCE
!       MATER  :  COEFFICIENTS MATERIAU A T+DT
!       NVI    :  NOMBRE DE VARIABLES INTERNES
!       DEPS   :  INCREMENT DE DEFORMATION
!       SIGD   :  CONTRAINTE  A T
!       VIND   :  VARIABLES INTERNES  A T
!   VAR SIGF   :  CONTRAINTE A T+DT  (IN -> ELAS, OUT -> PLASTI)
!   OUT VINF   :  VARIABLES INTERNES A T+DT
!                 (CUMUL DECOUPAGE)
!       NDEC   :  NOMBRE DE DECOUPAGE
!       IRET   :  CODE RETOUR DE  L'INTEGRATION DE LA LOI DE HUJEUX
!                    IRET=0 => PAS DE PROBLEME
!                    IRET=1 => ECHEC
!       ETATF  :  ETAT PLASTIQUE OU ELASTIQUE DU POINT DE GAUSS
!   ------------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/hujact.h"
#include "asterfort/hujmid.h"
#include "asterfort/hujpot.h"
#include "asterfort/hujpre.h"
#include "asterfort/hujprj.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndt, ndi, nvi, ndec, iret, kpg, ksp
    integer(kind=8) :: i, k, ifm, niv, nsubd, maj, niter
    integer(kind=8) :: nvimax, idec, imat, indi(7)
    parameter(nvimax=50)
    real(kind=8) :: deps(6), vins(nvimax)
    real(kind=8) :: sigd(6), sigf(6), predic(6), ptrac
    real(kind=8) :: sigd0(6), deps0(6), predi0(6)
    real(kind=8) :: vind(*), vinf(*), vind0(nvimax)
    real(kind=8) :: mater(22, 2), crit(*)
    real(kind=8) :: pf, qf, vec(3), zero, pref
    real(kind=8) :: deux
    real(kind=8) :: tol, asig, adsig, tole1, rtrac
    aster_logical :: chgmec, noconv, aredec, stopnc, negmul(8)
    aster_logical :: subd, rdctps, loop, impose
    character(len=8) :: mod
    character(len=7) :: etatf
    character(len=*) :: fami
    aster_logical :: debug, plas, try, nodec, tract
!
    common/tdim/ndt, ndi
    common/meshuj/debug
!
    data zero, deux/0.d0, 2.d0/
    data tol, tole1/.5d0, 1.d-7/
!
    if (nvi .gt. nvimax) then
        call utmess('F', 'COMPOR1_1')
    end if
!
    try = .true.
    loop = .false.
    nodec = .false.
    niter = int(crit(1))
!
    do i = 1, 7
        indi(i) = 0
    end do
!
    ptrac = mater(21, 2)
    pref = mater(8, 2)
    rtrac = abs(1.d-6*pref)
!
    call infniv(ifm, niv)
!
! ----  SAUVEGARDE DES GRANDEURS D ENTREE INITIALES
    predi0(1:ndt) = sigf(1:ndt)
    sigd0(1:ndt) = sigd(1:ndt)
    deps0(1:ndt) = deps(1:ndt)
    vind0(1:nvi) = vind(1:nvi)
!
!  ARRET OU NON EN NON CONVERGENCE INTERNE
!  ---------------------------------------
    stopnc = .false.
!     if (int(crit(1)) .lt. 0) then
!         stopnc = .true.
!     else
!         stopnc = .false.
!     endif
!
!
!  INITIALISATION DES VARIABLES DE REDECOUPAGE
!  -------------------------------------------
!  INT( CRIT(1) ) = 0  1 OU -1 ==> PAS DE REDECOUPAGE DU PAS DE TEMPS
    if ((niter .eq. 0) .or. (niter .eq. -1) .or. (niter .eq. 1)) then
        ndec = 1
        nodec = .true.
        aredec = .true.
        noconv = .false.
!
!
! ---- INT( CRIT(1) ) < -1 ==> REDECOUPAGE DU PAS DE TEMPS
!                            SI NON CONVERGENCE
    else if (niter .lt. -1) then
        ndec = 1
        aredec = .false.
        noconv = .false.
!
!
! ---- INT( CRIT(5) ) > 1 ==> REDECOUPAGE IMPOSE DU PAS DE TEMPS
    else if (niter .gt. 1) then
        ndec = niter
        aredec = .true.
        noconv = .false.
    end if
!
!
!  POINT DE RETOUR EN CAS DE DECOUPAGE POUR DR/R > TOLE OU
!  APRES UNE NON CONVERGENCE, POUR  INT(CRIT(1)) < -1
!
    subd = .false.
    rdctps = .false.
!
! ============================
500 continue
! ============================
!
    if (noconv) then
        ndec = 3
        iret = 0
        aredec = .true.
        noconv = .false.
    end if
!
    if (subd) then
        ndec = nsubd
        iret = 0
        aredec = .true.
        try = .false.
    end if
!
    if (rdctps) then
        ndec = 3
        iret = 0
        aredec = .true.
    end if
!
!   RESTAURATION DE SIGD VIND DEPS ET PREDIC ELAS SIGF
!   EN TENANT COMPTE DU DECOUPAGE EVENTUEL
    sigd(1:ndt) = sigd0(1:ndt)
    vind(1:nvi) = vind0(1:nvi)
!
    do i = 1, ndt
        deps(i) = deps0(i)/ndec
    end do
!
    call hujpre(fami, kpg, ksp, 'ELASTIC', mod, &
                imat, mater, deps, sigd, &
                sigf, vind, iret)
!
!KH ON IMPOSE UNE CONDITION SUR LA VARIATION DE DSIGMA
!   < 50% de SIGMA_INIT
    asig = 0.d0
    adsig = 0.d0
!
    do i = 1, ndt
        asig = asig+sigd(i)**2.d0
        adsig = adsig+(sigf(i)-sigd(i))**2.d0
    end do
    asig = sqrt(asig)
    adsig = sqrt(adsig)
!
    if ((-asig/pref) .gt. tole1) then
        nsubd = nint(adsig/asig/tol)
        if (nsubd .gt. 1) then
            ndec = min(ndec*nsubd, 100)
            do i = 1, ndt
                deps(i) = deps0(i)/ndec
            end do
        end if
    end if
!
!
    if (debug) then
!
        write (ifm, *)
        write (ifm, '(A)') '================== HUJRES =================='
        write (ifm, *)
        write (ifm, 1001) 'NDEC=', ndec
!
    end if
!
!               INIT BOUCLE SUR LES DECOUPAGES
!  =============================================================
    do idec = 1, ndec
!
! --- MISE A JOUR DU COMPTEUR D'ITERATIONS LOCALES
        loop = .false.
        vind(35) = vind0(35)+idec
!
        if (debug) then
!
            write (ifm, *)
            write (ifm, 1001) '%% IDEC=', idec
!
        end if
!
        maj = 0
        predic(1:ndt) = sigf(1:ndt)
!
        do k = 1, 8
            negmul(k) = .false.
        end do
!
! ---> SAUVEGARDE DES SURFACES DE CHARGE AVANT MODIFICATION
        vins(1:nvi) = vind(1:nvi)
! ---> DEFINITION DU DOMAINE POTENTIEL DE MECANISMES ACTIFS
!
        call hujpot(mod, mater, vind, deps, sigd, &
                    sigf, etatf, rdctps, iret, aredec)
        if (iret .eq. 1) goto 9999
!
        if (rdctps .and. (.not. aredec)) then
            goto 500
        else if (rdctps .and. (aredec)) then
            iret = 1
            goto 9999
        end if
!
!
! ---> SI ELASTICITE PASSAGE A L'ITERATION SUIVANTE
        plas = .false.
        if (etatf .eq. 'ELASTIC') then
            do i = 1, 3
                call hujprj(i, sigf, vec, pf, qf)
                if (((pf+deux*rtrac-ptrac)/abs(pref)) .ge. zero) then
                    plas = .true.
                    etatf = 'PLASTIC'
                end if
            end do
!
            if (.not. plas) then
                do i = 1, nvi
                    vinf(i) = vind(i)
                end do
                chgmec = .false.
                goto 40
            end if
        end if
!
!
! ---> SINON RESOLUTION VIA L'ALGORITHME DE NEWTON
        chgmec = .false.
!
        vinf(1:nvi) = vind(1:nvi)
        if (debug) write (6, *) 'HUJRES - VINF =', (vinf(i), i=24, 31)
!
100     continue
!--->   RESOLUTION EN FONCTION DES MECANISMES ACTIVES
!       MECANISMES ISOTROPE ET DEVIATOIRE
!-----------------------------------------------------
!
        call hujmid(mod, crit, mater, nvi, deps, &
                    sigd, sigf, vind, vinf, noconv, &
                    aredec, stopnc, negmul, iret, subd, &
                    loop, nsubd, indi, tract)
!
! --- ON CONTROLE QUE LES PRESSIONS ISOTROPES PAR PLAN
!     NE PRESENTENT PAS DE TRACTION
        do i = 1, ndi
            call hujprj(i, sigf, vec, pf, qf)
            if (((pf+rtrac-ptrac)/abs(pref)) .gt. zero) then
                noconv = .true.
                if (debug) write (6, '(A)') 'HUJRES :: SOL EN TRACTION'
            end if
        end do
! --- SI TRACTION DETECTEE ET NON CONVERGENCE, ON IMPOSE
! --- ETAT DE CONTRAINTES ISOTROPE
        if ((noconv) .and. (tract)) then
            noconv = .false.
            do i = 1, 3
                sigf(i) = -deux*rtrac
                sigf(3+i) = zero
                vind(23+i) = zero
                vind(27+i) = zero
                vind(5+4*i) = zero
                vind(6+4*i) = zero
                vind(7+4*i) = zero
                vind(8+4*i) = zero
            end do
            vinf(1:nvi) = vind(1:nvi)
            iret = 0
            if (debug) write (6, '(A)') 'HUJRES :: CORRECTION SIGMA'
        end if
! --- SI PROCHE TRACTION ET NON CONVERGENCE, ON IMPOSE
! --- ETAT DE CONTRAINTES ISOTROPE
        if (noconv) then
            impose = .false.
            do i = 1, ndi
                call hujprj(i, sigf, vec, pf, qf)
                if ((abs(pf-ptrac)/abs(pref)) .lt. 1.d-5) then
                    impose = .true.
                    if (debug) write (6, '(A)') 'HUJRES :: SOL LIQUEFIE'
                end if
            end do
            if (impose) then
                noconv = .false.
                do i = 1, 3
                    sigf(i) = -deux*rtrac
                    sigf(3+i) = zero
                    vind(23+i) = zero
                    vind(27+i) = zero
                    vind(5+4*i) = zero
                    vind(6+4*i) = zero
                    vind(7+4*i) = zero
                    vind(8+4*i) = zero
                end do
                vinf(1:nvi) = vind(1:nvi)
                iret = 0
            end if
        end if
!
! --- PRISE DE DECISION POUR LA SUITE DE L'ALGORITHME
! --- ECHEC DE L'INTEGRATION
        if (iret .eq. 1) goto 9999
!
! --- REDECOUPAGE LOCAL SI POSSIBLE
        if ((noconv .or. subd) .and. (.not. aredec)) goto 500
!
! --- REDECOUPAGE LOCAL SI ENCORE POSSIBLE
        if (subd .and. try .and. nodec) goto 500
!
! --- NON CONVERGENCE D'OU ECHEC D'INTEGRATION LOCAL
        if (noconv) then
            if (debug) write (6, '(A)') 'HUJRES :: PROBLEME AVEC NOCONV '
            iret = 1
            goto 9999
        end if
!
! --- REPRISE A CE NIVEAU SI MECANISME SUPPOSE ELASTIQUE
40      continue
!
! --- VERIFICATION DES MECANISMES SUPPOSES ACTIFS
!                 DES MULTIPLICATEURS PLASTIQUES
!                 DES SEUILS EN DECHARGE
!
        call hujact(mater, vind, vinf, vins, sigd, &
                    sigf, negmul, chgmec, indi)
        if (iret .eq. 1) goto 9999
!
! - SI ON MODIFIE LE DOMAINE POTENTIEL DE MECANISMES ACTIVES : RETOUR
!   AVEC UNE CONDITION DE CONVERGENCE DE LA METHODE DE NEWTON
        if (chgmec) then
            if (.not. aredec) then
! --- REDECOUPAGE ACTIVE S'IL NE L'ETAIT PAS ENCORE
                chgmec = .false.
                noconv = .true.
                goto 500
            else
! --- REPRISE DE L'ITERATION EN TENANT COMPTE DES MODIFICATIONS
                if (debug) write (ifm, '(A)') 'HUJRES :: CHANGEMENT DE MECANISME'
                chgmec = .false.
!
! --- REINITIALISATION DE SIGF A LA PREDICTION ELASTIQUE PREDIC
!            CALL LCEQVE (PREDIC, SIGF)
! --- PREDICTEUR MIS A JOUR TENANT COMPTE DE L'ETAT PRECEDEMMENT OBTENU
                maj = maj+1
                if (maj .lt. 5) then
                    loop = .true.
                    vinf(1:nvi) = vind(1:nvi)
                    goto 100
                else
                    if (debug) write (6, '(A)') 'HUJRES :: SOLUTION EXPLICITE'
                end if
            end if
        end if
! End - CHGMEC
! --- S'IL N'Y A PAS DE CHANGEMENT DE MECANISME, ON POURSUIT
!
        maj = 0
        if (idec .lt. ndec) then
            sigd(1:ndt) = sigf(1:ndt)
            do i = 1, nvi
                vind(i) = vinf(i)
            end do
            do i = 1, ndt
                deps(i) = deps0(i)/ndec
            end do
! --- APPLICATION DE L'INCREMENT DE DÉFORMATIONS, SUPPOSE ELASTIQUE
            call hujpre(fami, kpg, ksp, 'ELASTIC', mod, &
                        imat, mater, deps, sigd, &
                        sigf, vind, iret)
        end if
!
    end do
! End - Boucle sur les redecoupages
!
!
9999 continue
!
1001 format(a, i3)
    if (debug) write (6, *) 'IRET - HUJRES =', iret
end subroutine
