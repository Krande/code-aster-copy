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

subroutine mxwell(fami, kpg, ksp, ndim, typmod, &
                  imate, instam, instap, deps, sigm, option, &
                  sigp, vip, dsidep, iret)
!
!     REALISE LA LOI DE MAXWELL ISOTROPE : VISC_MAXWELL
!
! IN  FAMI    : FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
! IN  KPG,KSP : NUMERO DU (SOUS)POINT DE GAUSS
! IN  NDIM    : DIMENSION DE L'ESPACE (3D=3,2D=2,1D=1)
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  INSTAM  : INSTANT T
! IN  INSTAP  : INSTANT T+DT
! IN  DEPS    : INCREMENT DE DEFORMATION TOTALE
!               SI C_PLAN DEPS(3) EST EN FAIT INCONNU (ICI:0)
!                 =>  ATTENTION LA PLACE DE DEPS(3) EST ALORS UTILISEE.
! IN  SIGM    : CONTRAINTES A T
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
! OUT IRET    : CODE RETOUR DE L'INTEGRATION DE LA LOI DE VOM MISES
!               = 1  => PAS DE PROBLEME
!               = 0  => ECHEC DANS L'INTEGRATION DE LA LOI
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verift.h"

    aster_logical :: cplan
    integer(kind=8) :: ndim, imate, kpg, ksp, iret, ndimsi
    integer(kind=8) :: k, l, icodre(3)
!
    real(kind=8) :: instam, instap, dt
    real(kind=8) :: rac2, valres(3)
    real(kind=8) :: em, num, e, nu, deuxGm, deuxG, bulkmodulusm, bulkmodulus
    real(kind=8) :: etadm, etad, etavm, etav
    real(kind=8) :: alpha, depsth(6), epsthe
    real(kind=8) :: depsmo, deps(6), depsdv(6), kron(6)
    real(kind=8) :: sigmmo, sigm(6), sigmdv(6), sigpmo, sigpdv(6), sigp(6), vip(*)
    real(kind=8) :: dsidep(6, 6)
!
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: nomres(4), option
!-----------------------------------------------------------------------
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!   data epsa/'EPSAXX','EPSAYY','EPSAZZ','EPSAXY','EPSAXZ','EPSAYZ'/
!------------------------------------------------------------------
!
!     -- 1 INITIALISATIONS :
!     ----------------------

    iret = 0
!
    cplan = typmod(1) .eq. 'C_PLAN'
!
    ndimsi = 2*ndim
!
    dt = instap-instam
    rac2 = sqrt(2.d0)
!
!
!     -- 2 RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Remplissage du vecteur nomres avec les nom des parametres
!
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!

! Recuperation des parametres elastiques E et NU aux instants - et +
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                3, nomres(1), valres(1), icodre(1), 2)
    em = valres(1)
    num = valres(2)
    alpha = valres(3)
!
    deuxGm = em/(1.d0+num)
    bulkmodulusm = em/(3.d0*(1.d0-2.d0*num))
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres(1), valres(1), icodre(1), 2)
    e = valres(1)
    nu = valres(2)
    deuxG = e/(1.d0+nu)
    bulkmodulus = e/(3.d0*(1.d0-2.d0*nu))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recuperation des parametres de fluage aux instants - et +
!
    nomres(1) = 'ETA_D'
    nomres(2) = 'ETA_V'
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'VISC_MAXWELL', 0, ' ', [0.d0], &
                2, nomres(1), valres(1), icodre(1), 2)
    etadm = valres(1)
    etavm = valres(2)
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'VISC_MAXWELL', 0, ' ', [0.d0], &
                2, nomres(1), valres(1), icodre(1), 2)
    etad = valres(1)
    etav = valres(2)

!    CALCUL DE DEPSMO ET DEPSDV (deformation moyenne et deviateur des deformations)
!     --------------------------------
!!
! - Get temperatures
!
    call verift(fami, kpg, ksp, '-', imate, &
                epsth_=epsthe)
!
!
    depsmo = 0.d0
    depsdv(:) = 0.d0
    depsth(:) = 0.d0
    do k = 1, ndimsi
        depsth(k) = deps(k)
    end do
    do k = 1, 3
        depsth(k) = depsth(k)-epsthe
        depsmo = depsmo+depsth(k)
    end do
    depsmo = depsmo/3.
    do k = 1, ndimsi
        depsdv(k) = depsth(k)-depsmo*kron(k)
    end do
!
!     CALCUL DE SIGMMO, SIGMDV :
!     -------------------------------------------------------
    sigmmo = 0.d0
    sigmdv(:) = 0.d0
    do k = 1, 3
        sigmmo = sigmmo+sigm(k)
    end do
    sigmmo = sigmmo/3.d0
    do k = 1, ndimsi
        sigmdv(k) = sigm(k)-sigmmo*kron(k)
    end do
!
!    CALCUL DE SIGPMO, SIGPDV, SIGP .... ):
!     -------------------------------------
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA' .or. option(1:16) &
        .eq. 'RIGI_MECA_IMPLEX') then
!
        sigpmo = 0.d0
        sigpdv(:) = 0.d0
        sigp(:) = 0.d0
!
    !! VISCO-ELASTIQUE
    !!
!        ! SIGPMO
        sigpmo = (depsmo*3.d0*bulkmodulus+sigmmo)/(1.d0+bulkmodulus*dt/etav)
!        ! SIGPDV
        do k = 1, ndimsi
            sigpdv(k) = (depsdv(k)*deuxG+sigmdv(k))/(1.d0+deuxG*dt/etad)
        end do
!        ! SIGP
        do k = 1, ndimsi
            sigp(k) = sigpdv(k)+sigpmo*kron(k)
        end do

        vip(1) = 0.d0
    end if
!
!      CALCUL DE DSIDEP(6,6) :
!     ----------------------------
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!
        dsidep(:, :) = 0.d0
!        ! Les 9 premiers termes
        do k = 1, 3
            do l = 1, 3
                dsidep(k, l) = dsidep(k, l)+bulkmodulus/(1.d0+bulkmodulus*dt/etav) &
                               -(1.d0/3.d0)*deuxG/(1.d0+deuxG*dt/etad)
            end do
        end do
        ! Les termes diagonaux
        do k = 1, ndimsi
            dsidep(k, k) = dsidep(k, k)+deuxG/(1.d0+deuxG*dt/etad)
        end do
!       --  CORRECTION POUR LES CONTRAINTES PLANES
        if (cplan) then
            do k = 1, ndimsi
                if (k .ne. 3) then
                    do l = 1, ndimsi
                        if (l .ne. 3) then
                            dsidep(k, l) = dsidep(k, l)-1.d0/dsidep(3, 3)*dsidep(k, 3)*dsidep(3, l)
                        end if
                    end do
                end if
            end do
        end if
!
    end if
!
end subroutine
