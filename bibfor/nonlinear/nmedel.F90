! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine nmedel(ndim, typmod, imate, deps, sigm,&
                  option, sigp, dsidep)
!
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
    integer :: ndim, imate
    character(len=8) :: typmod(*)
    character(len=16) :: option
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6), sigp(6), dsidep(6, 6)
!
! ----------------------------------------------------------------------
!     LOI ELASTIQUE POUR L'ELEMENT A DISCONTINUITE
! ----------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  DEPS    : INCREMENT DE DEFORMATION
!               SI C_PLAN DEPS(3) EST EN FAIT INCONNU (ICI:0)
!                 =>  ATTENTION LA PLACE DE DEPS(3) EST ALORS UTILISEE.
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE (INUTILISE POUR RAPH_MECA)
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!-----------------------------------------------------------------------
!
    aster_logical :: cplan
    real(kind=8) :: deuxmu
    real(kind=8) :: valres(3)
    real(kind=8) :: depsmo, e, nu, troisk
    real(kind=8) :: kron(6), depsdv(6)
    integer :: ndimsi
    integer :: k, j, kpg, spt
    integer :: icodre(3)
    character(len=8) :: fami, poum
    character(len=16) :: nomres(3)
    data        kron/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
!
!     INITIALISATIONS
!     ---------------
!
    cplan = typmod(1) .eq. 'C_PLAN'
    ndimsi = 2*ndim
!
!     RECUPERATION DES CARACTERISTIQUES
!     ---------------------------------
    nomres(1)='E'
    nomres(2)='NU'
    fami='FPG1'
    kpg=1
    spt=1
    poum='+'
!
    call rcvalb(fami, kpg, spt, poum, imate,&
                ' ', 'ELAS', 0, ' ', [0.d0],&
                2, nomres(1), valres(1), icodre(1), 2)
    e = valres(1)
    nu = valres(2)
!
    deuxmu = e/(1.d0+nu)
    troisk = e/(1.d0-2.d0*nu)
!
!
    if (cplan) deps(3)=-nu/(1.d0-nu)*(deps(1)+deps(2)) +(1.d0+nu)/(1.d0-nu)
    depsmo = (deps(1)+deps(2)+deps(3))/3.d0
    do k = 1, ndimsi
        depsdv(k) = deps(k) - depsmo * kron(k)
    end do
!
!
    do k = 1, ndimsi
        sigp(k) = sigm(k)+deuxmu*depsdv(k)+troisk*depsmo*kron(k)
    end do
!
!
!      CALCUL DE DSIDEP(6,6) :
!     ------------------------
!
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
        call r8inir(36, 0.d0, dsidep, 1)
!
        do k = 1, 3
            do j = 1, 3
                dsidep(k,j) = troisk/3.d0-deuxmu/(3.d0)
            end do
        end do
        do k = 1, ndimsi
            dsidep(k,k) = dsidep(k,k) + deuxmu
        end do
!
    endif
!
end subroutine
