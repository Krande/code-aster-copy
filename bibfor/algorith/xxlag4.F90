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
subroutine xxlag4(ffc, idepl, idepm, lact, ndim, &
                  nnol, pla, lamb, nvec, champ)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
!
!
! Routine utilitaire de calcul d'un champ
!
! IN FFC    : FONCTIONS DE FORME DE CONTACT
! IN IDEPL  : ADRESSE DEPLACEMENT COURANT
! IN IDEPM  : ADRESSE DEPLACEMENT INSTANT -
! IN LACT   : DDL DE LAGRANGE ACTIF OU NON
! IN NDIM   : DIMENSION DU MODELE
! IN NNOL   : NOMBRE DE NOEUDS EL PARENT PORTEURS DE DDL LAGRANGE
! IN PLA    : PLACE DES DDLS DE LAGRANGE
! OUT LAMB  : CHAMP DEMANDE AU POINT DE GAUSS
! IN NVEC   : NOMBRE D ADRESSES DEPLACEMENT: 1 OU 2
! IN CHAMP  : NOM DU CHAMP: LAMBDA, MU OU W
    integer :: i, idepl, idepm
    integer :: j, lact(8), ndim, nli, nnol
    integer :: pla(27), pli, nvec, indcha
    real(kind=8) :: ffc(8), ffi, lamb(3)
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
    do i = 1, nnol
        pli = pla(i)
        ffi = ffc(i)
        nli = lact(i)
        if (nli .eq. 0) goto 1
        do j = 1, ndim
            lamb(j) = lamb(j)+ffi*zr(idepl-1+indcha*ndim+pli-1+j)
            if (nvec .eq. 2) then
                lamb(j) = lamb(j)+ffi*zr(idepm-1+indcha*ndim+pli-1+j)
            end if
        end do
1       continue
    end do
end subroutine
