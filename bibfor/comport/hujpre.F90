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

subroutine hujpre(fami, kpg, ksp, etat, mod, &
                  imat, mater, deps, sigd, &
                  sigf, vind, iret)
    implicit none
!                   CALCUL DE LA PREDICTION EN CONTRAINTE
!       ================================================================
!       IN      ETAT    COMPORTEMENT DU POINT DU POINT DE CALCUL
!                               'ELASTIC'     > ELASTIQUE
!                               'PLASTIC'     > PLASTIQUE
!               MOD     TYPE DE MODELISATION
!               DEPS    INCREMENT DE DEFORMATION TOTALE
!               SIGD    CONTRAINTE A T
!               VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!               OPT     OPTION DE CALCUL A FAIRE
!                               'RIGI_MECA_TANG'> DSDE(T)
!                               'FULL_MECA'     > DSDE(T+DT) , SIG(T+DT)
!                               'RAPH_MECA'     > SIG(T+DT)
!       OUT     SIGF    CONTRAINTE A T+DT
!               IRET    CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/hujela.h"
#include "asterfort/hujprj.h"
#include "asterfort/hujtid.h"
#include "asterfort/tecael.h"
#include "asterfort/trace.h"
    integer(kind=8) :: ndt, ndi, imat, iret, iadzi, iazk24, i, kpg, ksp
    real(kind=8) :: vind(*)
    real(kind=8) :: deps(6), dev(3), pf(3), q, pd(3), dp(3)
    real(kind=8) :: sigd(6), sigf(6), dsig(6), dsde(6, 6), rtrac
    real(kind=8) :: mater(22, 2), i1, d13, tole1, un, zero
    real(kind=8) :: ptrac, pref, maxi, cohes, factor
    character(len=7) :: etat
    character(len=8) :: mod, nomail
    character(len=*) :: fami
    logical :: debug
!
    common/tdim/ndt, ndi
    common/meshuj/debug
!
    data un, zero/1.d0, 0.d0/
    data d13, tole1/0.33333333334d0, 1.0d-7/
!
!
    pref = mater(8, 2)
    ptrac = mater(21, 2)
    rtrac = abs(pref*1.d-6)
!
    if (etat .eq. 'ELASTIC') then
!
        call hujela(mod, mater, deps, sigd, sigf, iret)
!
    else if (etat .eq. 'PLASTIC') then
!
        call hujtid(fami, kpg, ksp, mod, imat, &
                    sigd, vind, dsde, iret)
        if (iret .eq. 0) then
            dsig(1:ndt) = matmul(dsde(1:ndt, 1:ndt), deps(1:ndt))
            sigf(1:ndt) = sigd(1:ndt)+dsig(1:ndt)
            i1 = d13*trace(ndi, sigf)
        else
            iret = 0
            i1 = -un
            if (debug) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                write (6, '(10(A))')&
     &     'HUJPRE :: ECHEC DANS LA PSEUDO-PREDICTION ELASTIQUE DANS ',&
     &     'LA MAILLE ', nomail
            end if
        end if
!
        if ((i1+un)/abs(pref) .ge. tole1) then
            if (debug) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                write (6, '(10(A))')&
     &      'HUJPRE :: TRACTION DANS LA PSEUDO-PREDICTION ELASTIQUE ',&
     &      'DANS LA MAILLE ', nomail
            end if
            call hujela(mod, mater, deps, sigd, sigf, iret)
        end if
!
    end if
!
!
! ---> CONTROLE QU'AUCUNE COMPOSANTE DU VECTEUR SIGF NE SOIT POSITIVE
    do i = 1, ndt
        dsig(i) = sigf(i)-sigd(i)
    end do
!
    maxi = un
    cohes = -rtrac+ptrac
    factor = un
!
    do i = 1, ndi
        call hujprj(i, sigf, dev, pf(i), q)
        call hujprj(i, sigd, dev, pd(i), q)
        call hujprj(i, dsig, dev, dp(i), q)
        if (pf(i) .gt. cohes .and. dp(i) .gt. tole1) then
            factor = (-pd(i)+cohes)/dp(i)
            if ((factor .gt. zero) .and. (factor .lt. maxi)) then
                maxi = factor
            end if
        end if
    end do
!
!
! ---> SI IL EXISTE PF(I)>0, ALORS MODIFICATION DE LA PREDICTION
    if (maxi .lt. un) then
        do i = 1, ndt
            dsig(i) = maxi*dsig(i)
        end do
        sigf(1:ndt) = sigd(1:ndt)+dsig(1:ndt)
        if (debug) then
            write (6, '(A,A,E12.5)')&
     &    'HUJPRE :: APPLICATION DE FACTOR POUR MODIFIER ',&
     &    'LA PREDICTION -> FACTOR =', maxi
            write (6, '(A,6(1X,E12.5))') 'SIGF =', (sigf(i), i=1, ndt)
        end if
    end if
!
end subroutine
