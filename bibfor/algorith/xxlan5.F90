! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine xxlan5(ino, idepl, idepm, idep0, lact, &
                  ndim, pla, lamb, nvec, champ)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
    integer :: ino
!
!
! Routine utilitaire de calcul d'un champ
!
! IN INO    : NUMERO DU NOEUD
! IN IDEPL  : ADRESSE INCREMENT PRECEDANT
! IN IDEPM  : ADRESSE DEPLACEMENT INSTANT -
! IN IDEP0  : ADRESSE INCREMENT ACTUEL
! IN LACT   : DDL DE LAGRANGE ACTIF OU NON
! IN NDIM   : DIMENSION DU MODELE
! IN PLA    : PLACE DES DDLS DE LAGRANGE
! OUT LAMB  : CHAMP DEMANDE
! IN NVEC   : NOMBRE D ADRESSES DEPLACEMENT: 1 OU 2
! IN CHAMP  : NOM DU CHAMP: LAMBDA, MU OU W
    integer :: idepl, idepm, idep0
    integer :: j, lact(8), ndim, nli
    integer :: pla(27), pli, nvec, indcha
    real(kind=8) :: lamb(3)
    character(len=8) :: champ
!
! --- RÉACTION CONTACT = SOMME DES FF(I).LAMBDA(I) POUR I=1,NNOL
! --- RÉACTION FROTT = SOMME DES FF(I).(LAMB1(I).TAU1+LAMB2(I).TAU2)
! --- (DEPDEL+DEPMOI)
    if (champ .eq. 'LAMBDA') then
        indcha = 0
    else if (champ .eq. 'W') then
        indcha = 1
    else if (champ .eq. 'MU') then
        indcha = 2
    else
        ASSERT(.false.)
    end if
    lamb(:) = 0.d0
    pli = pla(ino)
    nli = lact(ino)
    if (nli .eq. 0) goto 1
    do j = 1, ndim
        lamb(j) = lamb(j)+zr(idepl-1+pli+indcha*ndim-1+j)
        if (nvec .ge. 2) then
            lamb(j) = lamb(j)+zr(idepm-1+pli+indcha*ndim-1+j)
        end if
        if (nvec .ge. 3) then
            lamb(j) = lamb(j)+zr(idep0-1+pli+indcha*ndim-1+j)
        end if
    end do
1   continue
end subroutine
